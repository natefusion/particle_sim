package particle_sim;
import rl "vendor:raylib"
import "core:fmt"
import "core:math"
import "core:math/linalg"
import "core:math/rand"
import "core:slice"
dynamic_viscosity_water :: 0.0010016; // 68 deg F kg/(ms)
dynamic_viscosity_mystery :: 1;
dynamic_viscosity_air ::    0.00001822;
density_water : f32 : 1000; // 22 deg C
density_iron : f32 : 7874; // kg/m^3
density_mystery :: 7874;
gravity : f32 : 9.81
line_width :: 2
bounds :: rl.Rectangle {10, 40, 800, 800};
restitution :: 0.5

// particle aggregation
// gamma = ((mu_0 · (71Am^2)^2)÷(2pi · (107.4nm)^3 · k_B · 296 K)) = 1.99 * 10^38
// diameter = 107.4 nm
// 107.4nm * 200,000 = 21.48mm = 2.148cm
// 200,000^3 = 8 * 10^15 particles in a (2.148cm)^3 cube

Particle :: struct {
    force        : [2]f32,
    velocity     : [2]f32,
    velocity_partial : [2]f32,
    position     : [2]f32,
    position_partial : [2]f32,
    mass         :    f32,
    radius       :    f32,
    density      :    f32,
}

Velocity_Field :: struct {
    position : [2]f32,
    width : f32,
    height : f32,
    field : [12][12][2]f32,
}

is_nan_vec2 :: proc(x: [2]f32) -> bool {
    return math.is_nan_f32(x[0]) || math.is_nan_f32(x[1])
}

sphere_volume :: proc(radius: f32) -> f32 {
    return 4 * math.PI / 3 * math.pow(radius, 3);
}

draw_velocity_field :: proc(field : ^Velocity_Field) {
    for row, i in field.field {
        rowlen := cast(f32)len(field.field);
        for col, j in row {
            collen := cast(f32)len(row);
            start : [2]f32 = {field.position.x + cast(f32)j*field.width/collen, field.position.y + cast(f32)i*field.height/rowlen};
            rl.DrawLineEx(start,
                          start + linalg.normalize(col)*20,
                          1,
                          rl.GREEN);

            rl.DrawCircleV(start, 2, rl.GREEN);
        }
    }
}

// Stokes Law stuff (https://en.wikipedia.org/wiki/Stokes%27_law)
// Assumumptions:
// Laminar flow
// no inertial effects (reynold's number is zero)
// homogenous material
// smooth surfaces
// particles no do interfere with each other (problem?)

F_d :: proc(dynamic_viscosity, radius: f32, velocity: [2]f32) -> [2]f32 {
    return -6 * math.PI * dynamic_viscosity * radius * velocity;
}

F_g :: proc(mass: f32) -> [2]f32 {
    return {0, mass * gravity};
}

x_axis_sort :: proc(a,b : Particle) -> bool {
    return a.position.x < b.position.x
}

intersection :: proc(a,b: [2]f32) -> bool {
    inside :: proc(a: f32, b: [2]f32) -> bool {
        return a >= b[0] && a <= b[1]
    }
    c1 := inside(a[0], b)
    c2 := inside(a[1], b)
    c3 := inside(b[0], a)
    c4 := inside(b[1], a)

    return c1 || c2 || c3 || c4
}

main :: proc() {
    rl.InitWindow(1000, 1000, "particle sim");

    rl.SetTargetFPS(rl.GetMonitorRefreshRate(rl.GetCurrentMonitor()));
    FPS := cast(f32)rl.GetMonitorRefreshRate(rl.GetCurrentMonitor());

    max_particles :: 16384
    ps : [dynamic]Particle;

    particle_template : Particle = {
        velocity = {1, 0},
        position = {500, 450},
        density  = density_iron,
        radius = math.pow_f32(10, -4),
    }
    
    particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density


    for i in 0..<max_particles {
        particle_template.position = {bounds.x+2 + rand.float32() * bounds.width, bounds.y+2 + rand.float32() * bounds.height}
        particle_template.velocity = {2 * rand.float32() - 1, 2 * rand.float32() - 1} * 10
        append(&ps, particle_template)
    }

    density_fluid : f32 = density_water

    damping_factor : f32 = 0.2;

    field_template : Velocity_Field;
    field_template.width = 100;
    field_template.height = 100;
    field_template.position = {10, 10};
    
    fields : [dynamic]Velocity_Field;
    selected_field := -1;
    // append(&fields, field_template);

    held_mouse_pos : [2]f32;

    for !rl.WindowShouldClose() {
        dt := rl.GetFrameTime()
        delete_field := -1;
        field_chosen := false;
        alt_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_ALT);
        ctrl_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_CONTROL);
        mousePos := rl.GetMousePosition();

        if alt_key_down && rl.IsKeyPressed(rl.KeyboardKey.D) {
            delete_field = selected_field;
        }

        if !alt_key_down do selected_field = -1;

        if rl.IsKeyPressed(rl.KeyboardKey.A) {
            field_template.position.x = mousePos.x - field_template.width/2 + 10;
            field_template.position.y = mousePos.y - field_template.height/2 + 10; 
            append(&fields, field_template);
        }

        if !ctrl_key_down && selected_field >= 0 {
            fields[selected_field].position = mousePos - held_mouse_pos;
        }

        mouse_particle_idx := -1

        for &p in ps {
            p.velocity_partial = 0
            p.position_partial = 0
            p.force = 0
        }

        slice.sort_by(ps[:], x_axis_sort)

        for i in 1..<len(ps) {
            p := &ps[i-1]
            p1 := &ps[i]

            l : [2]f32 = {-1,1}*p.radius + p.position.x
            l1 : [2]f32 = {-1,1}*p1.radius + p1.position.x
            if intersection(l, l1) {
                if rl.CheckCollisionCircles(p1.position, p.radius, p.position, p.radius) {
                    dposition1 := p.position-p1.position
                    dposition2 := p1.position-p.position
                    distance := linalg.length(dposition1)
                    distance_squared := math.pow(distance, 2)
                    v1 := (2 * p1.mass) / (p.mass + p1.mass) * linalg.dot(p.velocity-p1.velocity, dposition1) / distance_squared * dposition1
                    v2 := (2 * p.mass) / (p.mass + p1.mass) * linalg.dot(p1.velocity-p.velocity, dposition2) / distance_squared * dposition2
                    p.velocity_partial -= v1
                    p1.velocity_partial -= v2
                    
                    overlap := p.radius + p1.radius - linalg.distance(p.position,p1.position)
                    norm := linalg.normalize(p.position-p1.position)
                    p.position_partial += overlap * norm
                    p1.position_partial -= overlap * norm
                }
            }
        }

        for &p, p_idx in ps {

            if rl.CheckCollisionPointCircle(mousePos, p.position, p.radius) {
                mouse_particle_idx = p_idx
            }
            avg : [2]f32= 0
            num := 0

            for &field, field_idx in fields {
                if alt_key_down && selected_field < 0 && rl.CheckCollisionPointRec(mousePos, rl.Rectangle {field.position.x, field.position.y, field.width, field.height}) {
                    selected_field = field_idx;
                    held_mouse_pos = mousePos - fields[selected_field].position;
                }

                for &row, i in field.field {
                    rowlen := cast(f32)len(field.field);

                    for &col, j in row {
                        collen := cast(f32)len(row);
                        start : [2]f32 = {field.position.x + cast(f32)j*field.width/collen, field.position.y + cast(f32)i*field.height/rowlen};
                        if ctrl_key_down && field_idx == selected_field {
                            arrow := start - mousePos;
                            col = -100 * linalg.vector_normalize(arrow);
                        }

                        diff := linalg.abs(p.position - start);
                        if diff.x < field.width/collen && diff.y < field.height/rowlen {
                            avg += col;
                            num += 1;
                        }
                    }
                }
            }

            ybounds1 := (p.position.y > (bounds.height + bounds.y - p.radius - line_width));
            ybounds2 :=  p.position.y < p.radius + line_width + bounds.y;
            xbounds1 := (p.position.x > (bounds.width + bounds.x - p.radius - line_width));
            xbounds2 := p.position.x < line_width + p.radius + bounds.x;
            // F_d : [2]f32 = 0//F_d(dynamic_viscosity_mystery, p.radius, p.velocity) if num > 0 else F_d(dynamic_viscosity_air, p.radius, p.velocity);
            F_d : [2]f32 = F_d(dynamic_viscosity_water, p.radius, p.velocity)
            F_g : [2]f32 = F_g(p.mass);

            p.force += F_d
            p.force += F_g

            if ybounds1 || ybounds2 {
                p.velocity_partial.y += -(1 + restitution) * p.velocity.y
                if ybounds1 do p.position.y = bounds.height + bounds.y - p.radius - line_width
                if ybounds2 do p.position.y = bounds.y + line_width + p.radius
            }

            if xbounds1 || xbounds2 {
                p.velocity_partial.x += -(1 + restitution) * p.velocity.x
                if xbounds1 do p.position.x = bounds.width + bounds.x - p.radius - line_width
                if xbounds2 do p.position.x = bounds.x + line_width + p.radius
            }

            p.force += ((avg / cast(f32)num) if num > 0 else 0) * p.mass
        }

        for &p in ps {
            p.velocity += p.velocity_partial + p.force / p.mass * dt
            p.position += p.velocity * dt + p.position_partial
        }


        rl.BeginDrawing();
        rl.ClearBackground(rl.RAYWHITE);

        // if (mouse_particle_idx >= 0) {
        //     p := ps[mouse_particle_idx]
        //     // rl.DrawText(rl.TextFormat("Acceleration: (%f, %f)", p.acceleration.x / pixels_per_meter, p.acceleration.y / pixels_per_meter), 14, 14, 20, rl.BLACK);
        //     // rl.DrawText(rl.TextFormat("Velocity: (%f, %f)", p.velocity.x / pixels_per_meter, p.velocity.y / pixels_per_meter), 14, 14+20+2, 20, rl.BLACK);
        //     // rl.DrawText(rl.TextFormat("Velocity: %f", linalg.length(p.velocity) / pixels_per_meter), 14, 14+40+2, 20, rl.BLACK);
        //     // rl.DrawText(rl.TextFormat("Position: (%f, %f)", p.position.x, p.position.y), 14, 14+60+2, 20, rl.BLACK);
        //     // rl.DrawText(rl.TextFormat("mass: %f", p.mass), 14, 14+80+2, 20, rl.BLACK);
        // }

        rl.DrawFPS(14, 14+120+2)

        for p in ps {
            rl.DrawCircleV(p.position, 1 if p.radius < 1 else p.radius, rl.BLACK);
        }
        rl.DrawRectangleLinesEx(bounds, line_width, rl.RED);

        for &field, i in fields {
            if i == selected_field {
                rl.DrawRectangleLinesEx(rl.Rectangle{field.position.x-10, field.position.y-10, field.width, field.height}, line_width, rl.RED)
            }
            draw_velocity_field(&field);
        }
        
        rl.EndDrawing();

        if delete_field >= 0 {
            selected_field = -1; ordered_remove(&fields, delete_field);
        }
    }

    rl.CloseWindow();
}
