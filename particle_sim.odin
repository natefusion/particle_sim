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
bounds :: rl.Rectangle {10, 10, 700, 700}
bounds_ugh :: Big_Rect {x=cast(f64)bounds.x, y=cast(f64)bounds.y, width=cast(f64)bounds.width, height=cast(f64)bounds.height}
restitution :: 0.5
radius_visual :: 5.0
NM_PER_PX :f64: 100_000_000.0

left_outlet_count := 0
right_outlet_count := 0

magnet_on := true


// if true, particles don't kiss
incel :: true

Wall :: struct {
    start: [2]f64,
    end: [2]f64,
    callback : proc(ps: ^Particle)
}


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
    radius_visual : f64,
    density      :    f64,
    disabled     : bool,
}

Force_Point :: struct {
    position : [2]f64,
    strength : [2]f64,
}

GetSide :: proc(p, a, b: [2]f64) -> int {
    crossProduct := (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);

    if crossProduct > 0 {
        return 1; // Point is on the "left" side of the line
    } else if (crossProduct < 0) {
        return -1; // Point is on the "right" side of the line
    } else {
        return 0; // Point is on the line
    }
}


draw_wall :: proc(wall : Wall) {
    start := [2]f32{cast(f32)wall.start.x, cast(f32)wall.start.y}
    end := [2]f32{cast(f32)wall.end.x, cast(f32)wall.end.y}
    rl.DrawLineEx(start, end, line_width, rl.PURPLE)
}

CheckCollisionPointCircle :: proc(point, center: [2]f64, radius: f64) -> (collision: bool) {
    distanceSquared := (point.x - center.x)*(point.x - center.x) + (point.y - center.y)*(point.y - center.y);
    collision = distanceSquared <= radius*radius
    return
}

sphere_volume :: proc(radius: f64) -> f64 {
    return 4 * math.PI / 3 * math.pow(radius, 3);
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
    out = mag * dir * 10000 // Gauss, not force, who cares
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
    F_m : [2]f64 = F_m(magnet, p.position) if magnet_on else 0
    F_g : [2]f64 = 0//gravity * p.mass
    p.force += F_d + F_m + F_g
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

    if CheckCollisionCircles(p.position, p.radius_visual, p1.position, p1.radius_visual) {
        diff := p.position - p1.position
        dist := linalg.length(diff)
        n := diff / dist
        delta := 0.5 * (dist - p.radius_visual - p1.radius_visual)
        p.position -= n * delta
        p1.position += n * delta

        if incel || !(p1_idx in p.kissing && p_idx in p1.kissing) && kiss_probability() > (1 - 0.9) {
            pv := p.position - p.position_old
            pv1 := p1.position - p1.position_old

            p.position_old = p.position - pv * (2 * p.mass) / (p.mass + p1.mass)
            p1.position_old = p1.position - pv1 * (2 * p1.mass) / (p.mass + p1.mass)
        }
    } else {
        if !incel && kiss_probability() > (1 - 0.1) {
            delete_key(&p.kissing, p1_idx)
            delete_key(&p1.kissing, p_idx)
        }
    }
}

closest_point_on_segment :: proc(point: [2]f64, line_start: [2]f64, line_end: [2]f64) -> [2]f64 {
    line_vec := line_end - line_start
    point_vec := point - line_start
    t := linalg.dot(point_vec, line_vec) / linalg.length2(line_vec)
    t = max(0, min(1, t))
    return line_start + t * line_vec
}

check_particle_wall_collision :: proc(particle: ^Particle, wall: Wall) {
    closest := closest_point_on_segment(particle.position, wall.start, wall.end)

    calc_new_velocity :: proc(particle: ^Particle, dist_vec, distance, point_on_wall: [2]f64) {
        normal := linalg.normalize(dist_vec) if linalg.length(dist_vec) != 0 else {1, 0}
        overlap := particle.radius_visual - distance - line_width/2.0
        particle.position = point_on_wall + distance * normal
        particle.position_old = particle.position - linalg.normalize(linalg.normalize(particle.position - particle.position_old) + normal)* linalg.length(particle.position - particle.position_old)
    }

    dist_vec := particle.position - closest
    distance : f64 = linalg.length(dist_vec)

    if distance < particle.radius_visual {
        calc_new_velocity(particle, dist_vec, distance, closest)

        endpoints := [2][2]f64{wall.start, wall.end};
        for endpoint in endpoints {
            dist_vec := particle.position - endpoint
            distance := linalg.length(dist_vec)
            
            if distance < particle.radius_visual {
                calc_new_velocity(particle, dist_vec, distance, endpoint)
            }
        }

        wall.callback(particle)
    }
}

point_in_geometry :: proc(p: [2]f64, walls: []Wall) -> bool {
    n := len(walls)
    inside := false
    p1 := walls[0].start
    for i in 0..<(n+1) {
        p2 := walls[i % n].start
        if p.y > math.min(p1.y, p2.y) {
            if p.y <= max(p1.y, p2.y) {
                if p.x <= max(p1.x, p2.x) {
                    if p1.y != p2.y {
                        xinters := (p.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x
                        if p.x < xinters {
                            inside = !inside
                        }
                    } else if p1.x == p2.x {
                        inside = !inside
                    }
                }
            }
        }
        p1 = p2
    }
    return inside
}

fill_geometry_with_points :: proc(ps: ^[dynamic]Force_Point, walls: []Wall) {
    for _ in 0..<5000 {
        p := [2]f64{rand.float64()*1000.0,rand.float64()*1000.0}
        if point_in_geometry(p, walls) do append(ps, Force_Point{position = p, strength = {0,-5e9}})
    }
}

draw_points :: proc(ps: []Force_Point) {
    for points in ps {
        p := [2]f32{cast(f32)points.position.x, cast(f32)points.position.y}
        s := [2]f32{cast(f32)points.strength.x, cast(f32)points.strength.y}
        rl.DrawCircleV(p, 2, rl.BROWN);
        rl.DrawLineEx(p, p + 10*linalg.normalize(s), 2, rl.BROWN) 
    }
}

change_force_point_strength :: proc(p: ^Force_Point, mouse_pos: [2]f64) {
    if CheckCollisionCircles(mouse_pos, 50, p.position, 1) {
        p.strength = -10000000000* {0,1}//linalg.vector_normalize(p.position - mouse_pos)
    }
}

null_callback :: proc(p: ^Particle) {}

inlet_callback :: proc(p: ^Particle) {}


left_outlet_callback :: proc(p: ^Particle) {
    p.disabled = true
    left_outlet_count += 1
}

right_outlet_callback :: proc(p: ^Particle) {
    right_outlet_count += 1
    p.disabled = true
}


main :: proc() {
    rl.InitWindow(1000, 1000, "particle sim");

    rl.SetTargetFPS(rl.GetMonitorRefreshRate(rl.GetCurrentMonitor()));
    FPS := cast(f64)rl.GetMonitorRefreshRate(rl.GetCurrentMonitor());

    walls : [dynamic]Wall;
    append(&walls,
           Wall {start = {200, 30},end = {97, 183}, callback = null_callback},
           Wall {start = {97, 183},end = {97, 717}, callback = null_callback},
           Wall {start = {97, 717},end = {191, 969}, callback = null_callback},
           Wall {start = {191, 969},end = {241, 969}, callback = inlet_callback}, // inlet 
           Wall {start = {241, 969},end = {513, 202}, callback = null_callback},
           Wall {start = {513, 202},end = {484, 30}, callback = null_callback},
           Wall {start = {484, 30},end = {434, 30}, callback = right_outlet_callback}, // right outlet
           Wall {start = {434, 30},end = {400, 227}, callback = null_callback},
           Wall {start = {400, 227},end = {344, 227}, callback = null_callback},
           Wall {start = {344, 227},end = {344, 170}, callback = null_callback},
           Wall {start = {344, 170},end = {250, 30}, callback = null_callback},
           Wall {start = {250, 30},end = {200, 30}, callback = left_outlet_callback}, // left outlet
          )

    force_points : [dynamic]Force_Point
    fill_geometry_with_points(&force_points, walls[:])


    max_particles :: 500
    ps : [dynamic]Particle;

    particle_template : Particle = {
        density = density_iron,
        radius = 100, // nm
        radius_visual = 5,
    }
    particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density
    fmt.println(particle_template.mass)

    if (max_particles > 1) {
        for i in 0..<max_particles {
            particle_template.radius = 10 * (10 + rand.float64())
            particle_template.radius_visual = 5 * (1 + rand.float64())
            particle_template.position_old = {rand.float64() * 1000, rand.float64() * 1000}
            particle_template.position = particle_template.position_old
            particle_template.kissing = make(type_of(particle_template.kissing))
            if point_in_geometry(particle_template.position, walls[:]) {
                please_append := true;
                for j in 0..<len(ps) {
                    if CheckCollisionCircles(particle_template.position, particle_template.radius_visual,
                                             ps[j].position, ps[j].radius_visual) {
                        please_append = false
                        break
                    }
                }
                if please_append do append(&ps, particle_template)
            }
        }
    } else if max_particles == 1{
        particle_template.radius = 100
        particle_template.radius_visual = 5
        particle_template.position_old = {300, 285}
        particle_template.position = particle_template.position_old - {-.0625,0}//{-0.0625,0.0625}
        particle_template.kissing = make(type_of(particle_template.kissing))
        append(&ps, particle_template)

        particle_template.position_old = {350, 290}
        particle_template.position = particle_template.position_old - {0,0}
        particle_template.kissing = make(type_of(particle_template.kissing))

        append(&ps, particle_template)
    }

    held_mouse_pos : [2]f64;

    mouse_particle_idx := -1

    magnet: [2]f64 = {700, 100};
    radius_magnet :: 20
    magnet_selected := false

    
    for !rl.WindowShouldClose() {
        dt := cast(f64)rl.GetFrameTime()
        alt_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_ALT);
        ctrl_key_down := rl.IsKeyDown(rl.KeyboardKey.LEFT_CONTROL);
        left_mouse := rl.IsMouseButtonPressed(rl.MouseButton.LEFT);
        mousePos32 := rl.GetMousePosition();
        mousePos : [2]f64 = {cast(f64)mousePos32.x, cast(f64)mousePos32.y}
        

        if left_mouse && CheckCollisionPointCircle(mousePos, magnet, radius_magnet) {
            magnet_selected = !magnet_selected
        }

        if magnet_selected {
            magnet = mousePos
        }

        if rl.IsKeyPressed(rl.KeyboardKey.M) {
            magnet_on = !magnet_on
        }

        for p, i in ps do if p.disabled do unordered_remove(&ps, i)

        for &p, p_idx in ps {
            if CheckCollisionPointCircle(mousePos, p.position, p.radius_visual) {
                if left_mouse do mouse_particle_idx = p_idx
            }
            
            p.force = 0
            avg : [2]f64 = 0
            num := 0

            for &fp in force_points {
                if ctrl_key_down {
                    change_force_point_strength(&fp, mousePos)
                }
                diff := linalg.abs(p.position - fp.position)
                if diff.x < 10 && diff.y < 10 {
                    avg += fp.strength
                    num += 1
                }
            }
            p.force += ((avg / cast(f64)num) if num > 0 else 0) * p.mass
        }
        
        for i in 0..<len(ps) {
            for j in i+1..<len(ps) {
                resolve_collision(ps[:], i, j)
            }

            for j in 0..<len(walls) {
                check_particle_wall_collision(&ps[i], walls[j]);
            }
        }

        for &p in ps {
            add_forces(&p, magnet)
            update(&p, 1.0/1000.0)
        }
        
        rl.BeginDrawing();
        rl.ClearBackground(rl.RAYWHITE);

        if mouse_particle_idx >= len(ps) do mouse_particle_idx = -1
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

        rl.DrawText(rl.TextFormat("left outlet: %d\nRight outlet: %d", left_outlet_count, right_outlet_count), 800, 10, 20, rl.BLACK)

        rl.DrawFPS(14, 14+120+2)

        for p, p_idx in ps {
            pos : [2]f32 = {cast(f32)p.position.x, cast(f32)p.position.y}
            rl.DrawCircleV(pos, cast(f32)p.radius_visual, rl.BLACK);
            // if p_idx == mouse_particle_idx do rl.DrawCircleLinesV(pos, radius_visual+1, rl.RED);
        }

        draw_points(force_points[:])


        for wall in walls {
            draw_wall(wall);
        }

        // draw_magnetic_field(magnet)

        mag := [2]f32{cast(f32)magnet.x, cast(f32)magnet.y}
        rl.DrawCircleV(mag, radius_magnet, rl.BLUE);
        if magnet_selected do rl.DrawCircleLinesV(mag, radius_magnet+1, rl.RED);

        rl.EndDrawing();
    }

    rl.CloseWindow();
}
