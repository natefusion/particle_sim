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
gravity : [2]f32 : {0, 9.81}
line_width :: 2
bounds :: rl.Rectangle {10, 160, 800, 800};
restitution :: 0.5
radius_visual :: 10.0

// particle aggregation
// gamma = ((mu_0 · (71Am^2)^2)÷(2pi · (107.4nm)^3 · k_B · 296 K)) = 1.99 * 10^38
// diameter = 107.4 nm
// 107.4nm * 200,000 = 21.48mm = 2.148cm
// 200,000^3 = 8 * 10^15 particles in a (2.148cm)^3 cube

Particle :: struct {
    force        : [2]f32,
    position_old : [2]f32,
    position     : [2]f32,
    mass         :    f32,
    radius       :    f32,
    density      :    f32,
}

Velocity_Field :: struct {
    position : [2]f32,
    width : f32,
    height : f32,
    field : [4][4][2]f32,
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

add_forces :: proc(p: ^Particle) {
    F_d : [2]f32 = F_d(dynamic_viscosity_water, p.radius, p.position-p.position_old)
    p.force += (F_d + gravity*p.mass)
}

update :: proc(p: ^Particle, dt: f32) {
    temp := p.position
    p.position = 2*p.position - p.position_old + p.force / p.mass * dt * dt
    p.position_old = temp
}

resolve_collision :: proc(p, p1: ^Particle) {
    if rl.CheckCollisionCircles(p.position, radius_visual, p1.position, radius_visual) {
        diff := p.position - p1.position
        dist := linalg.length(diff)
        n := diff / dist
        delta := 0.5 * (dist - radius_visual - radius_visual)
        p.position -= n * delta
        p1.position += n * delta
    }
}

constrain :: proc(p: ^Particle) {
    ybounds1 := (p.position.y > (bounds.height + bounds.y - radius_visual - line_width));
    ybounds2 :=  p.position.y < radius_visual + line_width + bounds.y;
    xbounds1 := (p.position.x > (bounds.width + bounds.x - radius_visual - line_width));
    xbounds2 := p.position.x < line_width + radius_visual + bounds.x;

    if ybounds1 || ybounds2 {
        if ybounds1 do p.position.y = bounds.height + bounds.y - radius_visual - line_width
        if ybounds2 do p.position.y = bounds.y + line_width + radius_visual
        p.position_old.y = 2 * p.position.y - p.position_old.y
    }

    if xbounds1 || xbounds2 {
        if xbounds1 do p.position.x = bounds.width + bounds.x - radius_visual - line_width
        if xbounds2 do p.position.x = bounds.x + line_width + radius_visual
        p.position_old.x = 2 * p.position.x - p.position_old.x
    }
}

main :: proc() {
    rl.InitWindow(1000, 1000, "particle sim");

    rl.SetTargetFPS(rl.GetMonitorRefreshRate(rl.GetCurrentMonitor()));
    FPS := cast(f32)rl.GetMonitorRefreshRate(rl.GetCurrentMonitor());

    max_particles :: 100
    ps : [dynamic]Particle;

    particle_template : Particle = {
        density = density_iron,
        radius = math.pow_f32(10, -5),
    }
    particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density

    for i in 0..<max_particles {
        particle_template.position_old = {bounds.x+2 + rand.float32() * bounds.width, bounds.y+2 + rand.float32() * bounds.height}
        particle_template.position = particle_template.position_old
        // particle_template.velocity = {2 * rand.float32() - 1, 2 * rand.float32() - 1} * 50
        append(&ps, particle_template)
    }

    field_template : Velocity_Field;
    field_template.width = 100;
    field_template.height = 100;
    field_template.position = {900, 900};
    
    fields : [dynamic]Velocity_Field;
    selected_field := -1;
    append(&fields, field_template);

    held_mouse_pos : [2]f32;

    mouse_particle_idx := -1

    for !rl.WindowShouldClose() {
        dt := rl.GetFrameTime()
        delete_field := -1;
        field_chosen := false;
        alt_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_ALT);
        ctrl_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_CONTROL);
        left_mouse := rl.IsMouseButtonPressed(rl.MouseButton.LEFT);
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


        
        for &p, p_idx in ps {
            if rl.CheckCollisionPointCircle(mousePos, p.position, radius_visual) {
                if left_mouse do mouse_particle_idx = p_idx
            }
            
            p.force = 0
            avg : [2]f32 = 0
            num := 0
            for &field, field_idx in fields {
                if alt_key_down && selected_field < 0 && rl.CheckCollisionPointRec(mousePos, rl.Rectangle {field.position.x, field.position.y, field.width, field.height}) {
                    selected_field = field_idx
                    held_mouse_pos = mousePos - field.position
                }

                for &row, i in field.field {
                    rowlen := cast(f32)len(field.field);

                    for &col, j in row {
                        collen := cast(f32)len(row);
                        start : [2]f32 = {field.position.x + cast(f32)j*field.width/collen, field.position.y + cast(f32)i*field.height/rowlen};
                        if ctrl_key_down && field_idx == selected_field {
                            arrow := start - mousePos;
                            col = -1000 * linalg.vector_normalize(arrow);
                        }

                        diff := linalg.abs(p.position - start);
                        if diff.x < field.width/collen && diff.y < field.height/rowlen {
                            avg += col;
                            num += 1;
                        }
                    }
                }
            }
            p.force += ((avg / cast(f32)num) if num > 0 else 0) * p.mass
            
            add_forces(&p)
            update(&p, dt)
        }
        
        for i in 0..<len(ps) {
            constrain(&ps[i])
            for j in i+1..<len(ps) {
                resolve_collision(&ps[i], &ps[j])
            }
        }
        
        rl.BeginDrawing();
        rl.ClearBackground(rl.RAYWHITE);

        if (mouse_particle_idx >= 0) {
            p := ps[mouse_particle_idx]
            v := (p.position - p.position_old) / dt
            a : = p.force / p.mass
            rl.DrawText(rl.TextFormat("Acceleration: (%f, %f) m/s^2", a.x, a.y), 14, 14, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Acceleration: %f m/s^2", linalg.length(a)), 14, 14+20+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Velocity: (%f, %f) m/s", v.x, v.y), 14, 14+40+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Velocity: %f m/s", linalg.length(v)), 14, 14+60+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Position: (%f, %f) m", p.position.x, p.position.y), 14, 14+80+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("mass: %.10f kg", p.mass), 14, 14+100+2, 20, rl.BLACK);
        }

        rl.DrawFPS(14, 14+120+2)

        for p, p_idx in ps {
            rl.DrawCircleV(p.position, radius_visual, rl.BLACK);
            if p_idx == mouse_particle_idx do rl.DrawCircleLinesV(p.position, radius_visual+1, rl.RED);
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
