package particle_sim;
import rl "vendor:raylib";
import "core:fmt";
import "core:math";
import "core:math/linalg";
dynamic_viscosity_water :: 0.0010016; // 68 deg F kg/(ms)
dynamic_viscosity_mystery :: 1;
dynamic_viscosity_air ::    0.00001822;
density_water :: 1000; // 22 deg C
density_iron :: 7874; // kg/m^3
density_mystery :: 7874;
pixels_per_meter : f32 : 30;
gravity : f32 : 9.81 * pixels_per_meter;

// particle aggregation
// gamma = ((mu_0 · (71Am^2)^2)÷(2pi · (107.4nm)^3 · k_B · 296 K)) = 1.99 * 10^38
// diameter = 107.4 nm
// 107.4nm * 200,000 = 21.48mm = 2.148cm
// 200,000^3 = 8 * 10^15 particles in a (2.148cm)^3 cube

Particle :: struct {
    force        : [2]f32,
    acceleration : [2]f32,
    velocity     : [2]f32,
    position     : [2]f32,
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

sphere_volume :: proc(radius: f32) -> f32 {
    return 4. * math.PI / 3. * math.pow(radius, 3.);
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
    return {0., mass * gravity};
}

main :: proc() {
    rl.InitWindow(1000, 1000, "particle sim");

    rl.SetTargetFPS(rl.GetMonitorRefreshRate(rl.GetCurrentMonitor()));
    FPS := cast(f32)rl.GetMonitorRefreshRate(rl.GetCurrentMonitor());
    dt := 1.0/FPS

    particle_template : Particle = {
        acceleration = {0, 0},
        velocity = {1*pixels_per_meter, 0},
        position = {500, 450},
        density  = density_mystery / math.pow(pixels_per_meter, 3),
        radius = .2 * pixels_per_meter,
    }
    particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density;

    ps : [dynamic]Particle;
    append(&ps, particle_template)

    particle_template.position.x += 50
    particle_template.velocity = {0, 0}
    append(&ps, particle_template)

    density_fluid : f32 = density_water / math.pow(pixels_per_meter, 3);

    bounds := rl.Rectangle {10, 400, 980, 100};
    line_width :: 2;
    damping_factor : f32 = 0.2;
    restitution :: 0.2

    field_template : Velocity_Field;
    field_template.width = 100;
    field_template.height = 100;
    field_template.position = {10, 10};
    
    fields : [dynamic]Velocity_Field;
    selected_field := -1;
    append(&fields, field_template);

    held_mouse_pos : [2]f32;

    for !rl.WindowShouldClose() {
        delete_field := -1;
        field_chosen := false;
        alt_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_ALT);
        ctrl_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_CONTROL);
        mousePos := rl.GetMousePosition();

        // if rl.IsKeyPressed(rl.KeyboardKey.H) {
        //     p.acceleration = 0;
        //     p.velocity = 0;
        // }

        if alt_key_down && rl.IsKeyPressed(rl.KeyboardKey.D) {

            delete_field = selected_field;
        }

        if !alt_key_down do selected_field = -1;

        // if rl.IsMouseButtonPressed(rl.MouseButton.LEFT) {
        //     p.position = mousePos;
        // }

        if rl.IsKeyPressed(rl.KeyboardKey.A) {
            field_template.position.x = mousePos.x - field_template.width/2 + 10;
            field_template.position.y = mousePos.y - field_template.height/2 + 10; 
            append(&fields, field_template);
        }

        if rl.IsKeyPressed(rl.KeyboardKey.N) {
            particle_template.position = mousePos
            append(&ps, particle_template)
        }

        if !ctrl_key_down && selected_field >= 0 {
            fields[selected_field].position = mousePos - held_mouse_pos;
        }

        mouse_particle_idx := -1

        for &p in ps do p.force = 0

        for &p, p_idx in ps {
            if rl.CheckCollisionPointCircle(mousePos, p.position, p.radius) {
                mouse_particle_idx = p_idx
            }
            avg : [2]f32= 0
            num := 0

            for &p1, p1_idx in ps[p_idx+1:] {
                if rl.CheckCollisionCircles(p1.position, p1.radius, p.position, p.radius) {
                    normal := linalg.normalize(p1.position - p.position)
                    distance := linalg.distance(p.position,p1.position)
                    minDistance := p.radius + p1.radius
                    repulsion := normal * (minDistance - distance)
                    
                    momentum := p.velocity * p.mass
                    momentum1 := p1.velocity * p1.mass
                    force := (momentum + momentum1) / (p.mass + p1.mass) * FPS
                    
                    p.force -= force + repulsion * p.mass * FPS * FPS
                    p1.force += force + repulsion * p1.mass * FPS * FPS
                }
            }
            
            
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
                            col = -5000 * linalg.vector_normalize(arrow);
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

            F_d : [2]f32 = 0//F_d(dynamic_viscosity_mystery, p.radius, p.velocity) if num > 0 else F_d(dynamic_viscosity_air, p.radius, p.velocity);
            F_g : [2]f32 = 0//F_g(p.mass);



            F_wall :[2]f32 = 0
            if ybounds1 || xbounds1 || ybounds2 || xbounds2 {
                e :: 1
                F_wall = -p.velocity * FPS * p.mass * (1 + restitution) - (F_d + F_g);
                if ybounds1 do p.position.y = bounds.height + bounds.y - p.radius - line_width
                if xbounds1 do p.position.x = bounds.width + bounds.x - p.radius - line_width
                if ybounds2 do p.position.y = bounds.y + line_width + p.radius
                if xbounds2 do p.position.x = bounds.x + line_width + p.radius
            }

            F_particle : [2]f32 = 0
            
            F_net : [2]f32 =  F_d + F_g + ((avg / cast(f32)num) if num > 0 else 0) + F_wall + p.force
            p.acceleration = F_net / p.mass;
            p.velocity += p.acceleration * dt;            
            p.position += p.velocity * dt
        }


        rl.BeginDrawing();
        rl.ClearBackground(rl.RAYWHITE);

        if (mouse_particle_idx >= 0) {
            p := ps[mouse_particle_idx]
            rl.DrawText(rl.TextFormat("Acceleration: (%f, %f)", p.acceleration.x / pixels_per_meter, p.acceleration.y / pixels_per_meter), 14, 14, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Velocity: (%f, %f)", p.velocity.x / pixels_per_meter, p.velocity.y / pixels_per_meter), 14, 14+20+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Velocity: %f", linalg.length(p.velocity) / pixels_per_meter), 14, 14+40+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Position: (%f, %f)", p.position.x, p.position.y), 14, 14+60+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("mass: %f", p.mass), 14, 14+80+2, 20, rl.BLACK);
        }
        


        // rl.DrawLineEx(p.position, p.position + F_d, line_width, rl.BLUE);
        for p in ps {
            rl.DrawCircleV(p.position, p.radius, rl.BLACK);
            rl.DrawLineEx(p.position, p.position + p.acceleration * p.mass / pixels_per_meter, line_width, rl.PURPLE);
        }
        rl.DrawRectangleLinesEx(bounds, line_width, rl.RED);

        for &field, i in fields {
            if i == selected_field {
                rl.DrawRectangleLinesEx(rl.Rectangle{field.position.x-10, field.position.y-10, field.width, field.height}, line_width, rl.RED)
            }
            draw_velocity_field(&field);
        }
        
        rl.EndDrawing();

        if rl.IsKeyPressed(rl.KeyboardKey.D) && mouse_particle_idx >= 0 {
            ordered_remove(&ps, mouse_particle_idx)
        }

        if delete_field >= 0 {
            selected_field = -1; ordered_remove(&fields, delete_field);
        }
    }

    rl.CloseWindow();
}
