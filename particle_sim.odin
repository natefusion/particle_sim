package particle_sim;
import rl "vendor:raylib"
import "core:fmt"
import "core:math"
import "core:math/linalg"
import "core:math/rand"
import "core:slice"
NM_PER_M : f64 : 10e9
NG_PER_KG : f64 : 10e12
dynamic_viscosity_water :: 0.0010016 * NG_PER_KG / NM_PER_M; // 68 deg Fahrenheit ng/(nm*s)
dynamic_viscosity_mystery :: 1 * NG_PER_KG / NM_PER_M; // ng/(nm*s)
dynamic_viscosity_air ::    0.00001822 * NG_PER_KG / NM_PER_M; // ng/(nm*s)
density_water : f64 : 1000 * NG_PER_KG / (NM_PER_M * NM_PER_M * NM_PER_M); // 22 deg C, ng/nm^3
density_iron : f64 : 7874 * NG_PER_KG / (NM_PER_M * NM_PER_M * NM_PER_M); // ng/nm^3
gravity : [2]f64 : {0, 9.81 * NM_PER_M} // nm/s^2
line_width :: 2
bounds :: rl.Rectangle {10, 160, 800, 800}
bounds_ugh :: Big_Rect {x=cast(f64)bounds.x, y=cast(f64)bounds.y, width=cast(f64)bounds.width, height=cast(f64)bounds.height}
restitution :: 0.5
radius_visual :: 10.0
NM_PER_PX :f64: 100_000_000.0


Big_Rect :: struct {
    width, height, x, y: f64
}

// particle aggregation
// gamma = ((mu_0 · (71Am^2)^2)÷(2pi · (107.4nm)^3 · k_B · 296 K)) = 1.99 * 10^38
// diameter = 107.4 nm
// 107.4nm * 200,000 = 21.48mm = 2.148cm
// 200,000^3 = 8 * 10^15 particles in a (2.148cm)^3 cube
Empty :: struct {}
Set_int :: map[int]Empty

Particle :: struct {
    force        : [2]f64,
    position_old : [2]f64,
    position     : [2]f64,
    kissing      : Set_int,
    mass         :    f64,
    radius       :    f64,
    density      :    f64,
}

Velocity_Field :: struct {
    position : [2]f64,
    width : f64,
    height : f64,
    field : [4][4][2]f64,
}

CheckCollisionPointCircle :: proc(point, center: [2]f64, radius: f64) -> (collision: bool) {
    distanceSquared := (point.x - center.x)*(point.x - center.x) + (point.y - center.y)*(point.y - center.y);
    collision = distanceSquared <= radius*radius
    return
}

sphere_volume :: proc(radius: f64) -> f64 {
    return 4 * math.PI / 3 * math.pow(radius, 3);
}

draw_velocity_field :: proc(field : ^Velocity_Field) {
    for row, i in field.field {
        rowlen := cast(f32)len(field.field);
        for col, j in row {
            col : [2]f32 = {cast(f32)col.x, cast(f32)col.y}
            collen := cast(f32)len(row);
            start : [2]f32 = {cast(f32)field.position.x + cast(f32)j*cast(f32)field.width/collen, cast(f32)field.position.y + cast(f32)i*cast(f32)field.height/rowlen};
            rl.DrawLineEx(start,
                          start + linalg.normalize(col)*20,
                          1,
                          rl.GREEN);

            rl.DrawCircleV(start, 2, rl.GREEN);
        }
    }
}

draw_magnetic_field :: proc(magnet: [2]f64) {
    for i:f32=bounds.x; i < bounds.x+bounds.width; i+=50 {
        for j:f32=bounds.y; j < bounds.y+bounds.height; j+=50 {
            x :[2]f64= {cast(f64)i,cast(f64)j}
            F := F_m(magnet, x)
            d := linalg.length([2]f32{cast(f32)F.x,cast(f32)F.y}) / 1
            rl.DrawCircleV({i,j}, d, rl.GREEN)
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

F_d :: proc(dynamic_viscosity, radius: f64, velocity: [2]f64) -> [2]f64 {
    return -6 * math.PI * dynamic_viscosity * radius * velocity;
}

F_m :: proc(magnet: [2]f64, particle: [2]f64) -> (out: [2]f64) {
    x := (magnet - particle) * NM_PER_PX
    // north facing sensor
    a :: 0.06941
    b :: -0.05959
    c :: 3.86
    d :: 2.568

    mag := a / math.pow_f64(linalg.length(x) - b, c) + d 
    dir := linalg.normalize(x)
    out = mag * dir * 100000 // Gauss, not force, who cares
    return
}

x_axis_sort :: proc(a,b : Particle) -> bool {
    return a.position.x < b.position.x
}

intersection :: proc(a,b: [2]f64) -> bool {
    inside :: proc(a: f64, b: [2]f64) -> bool {
        return a >= b[0] && a <= b[1]
    }
    c1 := inside(a[0], b)
    c2 := inside(a[1], b)
    c3 := inside(b[0], a)
    c4 := inside(b[1], a)

    return c1 || c2 || c3 || c4
}

add_forces :: proc(p: ^Particle, magnet: [2]f64) {
    F_d : [2]f64 = F_d(dynamic_viscosity_water, p.radius, p.position-p.position_old)
    F_m : [2]f64 = F_m(magnet, p.position)
    p.force += (F_d + F_m + gravity*p.mass)
}

update :: proc(p: ^Particle, dt: f64) {
    temp := p.position
    a := p.force / p.mass
    x := a * dt * dt / NM_PER_PX
    p.position = 2*p.position - p.position_old + x
    p.position_old = temp
}

kiss_probability :: proc() -> f64 {
    return math.max(rand.float64(), rand.float64())
}

CheckCollisionCircles :: proc(center1: [2]f64, radius1: f64, center2: [2]f64, radius2: f64) -> (collision: bool) {
    dx := center2.x - center1.x      // X distance between centers
    dy := center2.y - center1.y      // Y distance between centers

    distanceSquared := dx*dx + dy*dy // Distance between centers squared
    radiusSum := radius1 + radius2

    collision = (distanceSquared <= (radiusSum*radiusSum))

    return
}

resolve_collision :: proc(ps: []Particle, p_idx, p1_idx: int) {
    p := &ps[p_idx]
    p1 := &ps[p1_idx]

    if CheckCollisionCircles(p.position, radius_visual, p1.position, radius_visual) {
        diff := p.position - p1.position
        dist := linalg.length(diff)
        n := diff / dist
        delta := 0.5 * (dist - radius_visual - radius_visual)
        p.position -= n * delta
        p1.position += n * delta

        if !(p1_idx in p.kissing && p_idx in p1.kissing) && kiss_probability() > (1 - 0.9) {
            p.kissing[p1_idx] = {}
            p1.kissing[p_idx] = {}
            pv := p.position - p.position_old
            pv1 := p1.position - p1.position_old

            v :[2]f64= (pv*p.mass + pv1*p1.mass)/(p.mass+p1.mass)
        
            p.position_old = p.position - v
            p1.position_old = p1.position - v
        }
    } else {
        if kiss_probability() > (1 - 0.1) {
            delete_key(&p.kissing, p1_idx)
            delete_key(&p1.kissing, p_idx)
        }
    }
}

constrain :: proc(p: ^Particle) {
    ybounds1 := (p.position.y > (bounds_ugh.height + bounds_ugh.y - radius_visual - line_width));
    ybounds2 :=  p.position.y < radius_visual + line_width + bounds_ugh.y;
    xbounds1 := (p.position.x > (bounds_ugh.width + bounds_ugh.x - radius_visual - line_width));
    xbounds2 := p.position.x < line_width + radius_visual + bounds_ugh.x;

    if ybounds1 || ybounds2 {
        if ybounds1 do p.position.y = bounds_ugh.height + bounds_ugh.y - radius_visual - line_width
        if ybounds2 do p.position.y = bounds_ugh.y + line_width + radius_visual
        p.position_old.y = 2 * p.position.y - p.position_old.y
    }

    if xbounds1 || xbounds2 {
        if xbounds1 do p.position.x = bounds_ugh.width + bounds_ugh.x - radius_visual - line_width
        if xbounds2 do p.position.x = bounds_ugh.x + line_width + radius_visual
        p.position_old.x = 2 * p.position.x - p.position_old.x
    }
}

main :: proc() {
    rl.InitWindow(1000, 1000, "particle sim");

    rl.SetTargetFPS(rl.GetMonitorRefreshRate(rl.GetCurrentMonitor()));
    FPS := cast(f64)rl.GetMonitorRefreshRate(rl.GetCurrentMonitor());

    max_particles :: 100
    ps : [dynamic]Particle;

    particle_template : Particle = {
        density = density_iron,
        radius = 100, // nm
    }
    particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density
    fmt.println(particle_template.mass)

    for i in 0..<max_particles {
        particle_template.position_old = {bounds_ugh.x+2 + rand.float64() * bounds_ugh.width, bounds_ugh.y+2 + rand.float64() * bounds_ugh.height}
        particle_template.position = particle_template.position_old
        // particle_template.velocity = {2 * rand.float32() - 1, 2 * rand.float32() - 1} * 50
        particle_template.kissing = make(type_of(particle_template.kissing))
        append(&ps, particle_template)
    }

    field_template : Velocity_Field;
    field_template.width = 100;
    field_template.height = 100;
    field_template.position = {900, 900};
    
    fields : [dynamic]Velocity_Field;
    selected_field := -1;
    append(&fields, field_template);

    held_mouse_pos : [2]f64;

    mouse_particle_idx := -1

    magnet: [2]f64 = {100, 100};
    radius_magnet :: 20
    magnet_selected := false

    for !rl.WindowShouldClose() {
        dt := cast(f64)rl.GetFrameTime()
        delete_field := -1;
        field_chosen := false;
        alt_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_ALT);
        ctrl_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_CONTROL);
        left_mouse := rl.IsMouseButtonPressed(rl.MouseButton.LEFT);
        mousePos32 := rl.GetMousePosition();
        mousePos : [2]f64 = {cast(f64)mousePos32.x, cast(f64)mousePos32.y}
        

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



        if left_mouse && CheckCollisionPointCircle(mousePos, magnet, radius_magnet) {
            magnet_selected = !magnet_selected
        }

        if magnet_selected {
            magnet = mousePos
        }
        
        for &p, p_idx in ps {
            if CheckCollisionPointCircle(mousePos, p.position, radius_visual) {
                if left_mouse do mouse_particle_idx = p_idx
            }
            
            p.force = 0
            avg : [2]f64 = 0
            num := 0
            for &field, field_idx in fields {
                if alt_key_down && selected_field < 0 && rl.CheckCollisionPointRec(mousePos32, rl.Rectangle {cast(f32)field.position.x, cast(f32)field.position.y, cast(f32)field.width, cast(f32)field.height}) {
                    selected_field = field_idx
                    held_mouse_pos = mousePos - field.position
                }

                for &row, i in field.field {
                    rowlen := cast(f64)len(field.field);

                    for &col, j in row {
                        collen := cast(f64)len(row);
                        start : [2]f64 = {field.position.x + cast(f64)j*field.width/collen, field.position.y + cast(f64)i*field.height/rowlen};
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
            p.force += ((avg / cast(f64)num) if num > 0 else 0) * p.mass
            

        }
        
        for i in 0..<len(ps) {
            constrain(&ps[i])
            for j in i+1..<len(ps) {
                resolve_collision(ps[:], i, j)
            }
        }

        for &p in ps {
            add_forces(&p, magnet)
            update(&p, dt)
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
            rl.DrawText(rl.TextFormat("mass: %.10f kg,\ttouching: %d", p.mass, len(p.kissing)), 14, 14+100+2, 20, rl.BLACK);
        }

        rl.DrawFPS(14, 14+120+2)

        for p, p_idx in ps {
            pos : [2]f32 = {cast(f32)p.position.x, cast(f32)p.position.y}
            rl.DrawCircleV(pos, radius_visual, rl.BLACK);
            // if p_idx == mouse_particle_idx do rl.DrawCircleLinesV(pos, radius_visual+1, rl.RED);
        }


        for &field, i in fields {
            if i == selected_field {
                rl.DrawRectangleLinesEx(rl.Rectangle{cast(f32)field.position.x-10, cast(f32)field.position.y-10, cast(f32)field.width, cast(f32)field.height}, line_width, rl.RED)
            }
            draw_velocity_field(&field);
        }

        // draw_magnetic_field(magnet)

        mag := [2]f32{cast(f32)magnet.x, cast(f32)magnet.y}
        rl.DrawCircleV(mag, radius_magnet, rl.BLUE);
        if magnet_selected do rl.DrawCircleLinesV(mag, radius_magnet+1, rl.RED);

        rl.DrawRectangleLinesEx(bounds, line_width, rl.RED);
        
        rl.EndDrawing();

        if delete_field >= 0 {
            selected_field = -1; ordered_remove(&fields, delete_field);
        }
    }

    rl.CloseWindow();
}
