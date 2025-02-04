package particle_sim;
import rl "vendor:raylib";
import "core:fmt";
import "core:math";
import "core:math/linalg";
dynamic_viscosity_water :: 0.0010016; // 68 deg F kg/(ms)
dynamic_viscosity_air ::    0.00001822;
density_water :: 1000; // 22 deg C
density_iron :: 7784;
pixels_per_meter : f32 : 30;
gravity : f32 : 9.81 * pixels_per_meter;

Particle :: struct {
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
    barrier : bool,
}

draw_velocity_field :: proc(field : ^Velocity_Field) {
    for row, i in field.field {
        rowlen := cast(f32)len(field.field);
        for col, j in row {
            collen := cast(f32)len(row);
            start : [2]f32 = {field.position.x + cast(f32)j*field.width/collen, field.position.y + cast(f32)i*field.height/rowlen};
            if !field.barrier {
                rl.DrawLineEx(start,
                              start + linalg.normalize(col)*20,
                              1,
                              rl.GREEN);
            }

            rl.DrawCircleV(start, 2, rl.RED if field.barrier else rl.GREEN);
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

    p : Particle = {
        acceleration = {0, 0},
        position = {500, 500},
        mass     = 1,
        radius   = .33 * pixels_per_meter,
        density  = density_iron / math.pow(pixels_per_meter, 3),
    };

    density_fluid : f32 = density_water / math.pow(pixels_per_meter, 3);

    bounds := rl.Rectangle {10, 10, 980, 980};
    line_width :: 2;
    damping_factor : f32 = 0.2;


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

        if rl.IsKeyPressed(rl.KeyboardKey.H) {
            p.acceleration = 0;
            p.velocity = 0;
        }

        if alt_key_down && rl.IsKeyPressed(rl.KeyboardKey.D) {
            delete_field = selected_field;
        }

        if !alt_key_down do selected_field = -1;

        if rl.IsMouseButtonPressed(rl.MouseButton.LEFT) {
            p.position = mousePos;
        }

        if rl.IsKeyPressed(rl.KeyboardKey.A) {
            field_template.position.x = mousePos.x - field_template.width/2 + 10;
            field_template.position.y = mousePos.y - field_template.height/2 + 10; 
            append(&fields, field_template);
        }

        if !ctrl_key_down && selected_field >= 0 {
            fields[selected_field].position = mousePos - held_mouse_pos;
            if rl.IsKeyPressed(rl.KeyboardKey.B) {
                fields[selected_field].barrier = !fields[selected_field].barrier;
            }
        }


        accel_to_add : [2]f32;
        for &field, field_idx in fields {
            if alt_key_down && selected_field < 0 && rl.CheckCollisionPointRec(mousePos, rl.Rectangle {field.position.x, field.position.y, field.width, field.height}) {
                selected_field = field_idx;
                held_mouse_pos = mousePos - fields[selected_field].position;
            }

            avg : [2]f32= 0
            num := 0
            for &row, i in field.field {
                rowlen := cast(f32)len(field.field);

                for &col, j in row {
                    collen := cast(f32)len(row);
                    start : [2]f32 = {field.position.x + cast(f32)j*field.width/collen, field.position.y + cast(f32)i*field.height/rowlen};
                    if ctrl_key_down && field_idx == selected_field {
                        arrow := start - mousePos;
                        col = -20 * linalg.vector_normalize(arrow);
                    }

                    diff := linalg.abs(p.position - start);
                    if diff.x < field.width/collen && diff.y < field.height/rowlen {
                        if !field.barrier {
                            avg += col;
                            num += 1;
                            p.acceleration = 0;
                        } else {
                            avg += linalg.normalize(p.velocity);
                            num += 1;
                        }
                    }
                }
            }

            if !field.barrier && num > 0 do p.velocity = avg / cast(f32)num;
            if field.barrier && num > 0 do accel_to_add = -10000 * linalg.normalize(p.velocity) * avg / cast(f32)num;
        }

        F_d : [2]f32 = F_d(dynamic_viscosity_air, p.radius, p.velocity);
        F_g : [2]f32 = F_g(p.mass);
        
        F_net : [2]f32 =  F_d + F_g;
        p.acceleration = F_net / p.mass + accel_to_add;
        p.velocity += p.acceleration * 1./FPS;
        p.position += p.velocity * 1./FPS;
        
        rl.BeginDrawing();

        rl.DrawText(rl.TextFormat("Acceleration: (%f, %f)", p.acceleration.x / pixels_per_meter, p.acceleration.y / pixels_per_meter), 14, 14, 20, rl.BLACK);
        rl.DrawText(rl.TextFormat("Velocity: (%f, %f)", p.velocity.x / pixels_per_meter, p.velocity.y / pixels_per_meter), 14, 14+20+2, 20, rl.BLACK);
        rl.DrawText(rl.TextFormat("Position: (%f, %f)", p.position.x, p.position.y), 14, 14+40+2, 20, rl.BLACK);
        
        rl.ClearBackground(rl.RAYWHITE);
        rl.DrawCircleV(p.position, p.radius, rl.BLACK);
        rl.DrawLineEx(p.position, p.position + F_d, line_width, rl.BLUE);
        rl.DrawRectangleLinesEx(bounds, line_width, rl.RED);

        ybounds := (p.position.y >= (bounds.height + bounds.y - p.radius - line_width)) || p.position.y <= p.radius + line_width;
        xbounds := (p.position.x >= (bounds.width + bounds.x - p.radius - line_width)) || p.position.x <= p.radius + line_width;
        if ybounds || xbounds {
            p.velocity *= -1 * damping_factor;
            if ybounds do p.position.y = bounds.height + bounds.y - p.radius - line_width;
            if xbounds do p.position.x = bounds.width + bounds.x - p.radius - line_width;
        }

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
