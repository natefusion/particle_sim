package particle_sim;
import rl "vendor:raylib"
import "core:fmt"
import "core:math"
import "core:math/linalg"
import "core:math/rand"
import "core:slice"
NM_PER_M : f64 : 1e9
NG_PER_KG : f64 : 1e12
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

magnet_on := false
dt :: 1.0/1000.0


Wall :: struct {
    start: [2]f64,
    end: [2]f64,
    callback : proc(ps: ^Particle),
    no_clip: bool,
}


Big_Rect :: struct {
    width, height, x, y: f64
}

// particle aggregation
// gamma = ((mu_0 · (71Am^2)^2)÷(2pi · (107.4nm)^3 · k_B · 296 K)) = 1.99 * 10^38
// diameter = 107.4 nm
// 107.4nm * 200,000 = 21.48mm = 2.148cm
// 200,000^3 = 8 * 10^15 particles in a (2.148cm)^3 cube

Particle :: struct {
    force        : [2]f64,
    position_old : [2]f64,
    position     : [2]f64,
    angular_position : f64,
    angular_position_old : f64,
    mass         :    f64,
    radius       :    f64,
    radius_visual : f64,
    density      :    f64,
    disabled     : bool,
}

Link :: struct {
    p: int,
    p1: int,
    length: f64,
}

Super_Particle :: [3]Link

Force_Point :: struct {
    position : [2]f64,
    strength : [2]f64,
}

GetSide :: proc(p, a, b: [2]f64) -> int {
    crossProduct := linalg.cross(b - a, p - a)

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

add_forces :: proc(p: ^Particle, magnet: [2]f64) {
    F_d : [2]f64 = F_d(dynamic_viscosity_water, p.radius, (p.position-p.position_old)/dt)
    F_m : [2]f64 = F_m(magnet, p.position) if magnet_on else 0
    p.force += F_d + F_m
}

set_angular_position :: proc(ps: []Particle, sp: ^Super_Particle) {
    // ugleh
    // the location for the three particles in a super particle is specified in make_super_particle, at the bottom of the function
    p1 : ^Particle = &ps[sp[0].p]
    p2 : ^Particle = &ps[sp[1].p]
    p3 : ^Particle = &ps[sp[1].p1]

    a := sp[0].length
    b := sp[1].length
    c := sp[2].length
    center := (a * p3.position + b * p1.position + c * p2.position) / (a + b + c)

    straight := [2]f64 {center.x+1, center.y}

    ref := straight-center
    r1 := p1.position-center
    r2 := p2.position-center
    r3 := p3.position-center
        
    temp := p1.angular_position
    p1.angular_position = math.mod(math.TAU + linalg.angle_between(r1, ref) * math.sign(linalg.cross(r1, ref)), math.TAU)
    p1.angular_position_old = temp

    temp = p2.angular_position
    p2.angular_position = math.mod(math.TAU + linalg.angle_between(r2, ref) * math.sign(linalg.cross(r2, ref)), math.TAU)
    p2.angular_position_old = temp

    temp = p3.angular_position
    p3.angular_position = math.mod(math.TAU + linalg.angle_between(r3, ref) * math.sign(linalg.cross(r3, ref)), math.TAU)
    p3.angular_position_old = temp
}

update :: proc(p: ^Particle, dt: f64) {
    temp := p.position
    a := p.force / p.mass
    x := a * dt * dt / NM_PER_PX
    p.position = 2*p.position - p.position_old + x
    p.position_old = temp
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
        pv_old := (p.position - p.position_old)
        pv1_old := (p1.position - p1.position_old)
        
        diff := p.position - p1.position
        dist := linalg.length(diff)
        n := diff / dist
        delta := 0.5 * (dist - p.radius_visual - p1.radius_visual)
        p.position -= n * delta
        p1.position += n * delta

        pv := pv_old - (2 * p1.mass) / (p.mass + p1.mass) * linalg.dot(pv_old - pv1_old, p.position - p1.position) / linalg.length2(p.position - p1.position) * (p.position - p1.position)
        pv1 := pv1_old - (2 * p.mass) / (p.mass + p1.mass) * linalg.dot(pv1_old - pv_old, p1.position - p.position) / linalg.length2(p1.position - p.position) * (p1.position - p.position)
        
        p.position_old = p.position - pv
        p1.position_old = p1.position - pv1
    }
}

update_link :: proc(ps: []Particle, link : Link) {
    p_pos := &ps[link.p].position
    p1_pos := &ps[link.p1].position

    diff := p_pos^ - p1_pos^
    dist := linalg.length(diff)
    diff_factor := (link.length - dist) / dist
    offset := diff * diff_factor * 0.5

    p_pos^ += offset
    p1_pos^ -= offset
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

circle_in_geometry :: proc(p: [2]f64, radius: f64, walls: []Wall) -> bool {
    for wall in walls {
        if CheckCollisionPointCircle(closest_point_on_segment(p, wall.start, wall.end), p, radius) {
            return false
        }
    }
    
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
        if point_in_geometry(p, walls) do append(ps, Force_Point{position = p, strength = {0,-5e11}})
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
        p.strength = -1e9* {0,1}
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

make_super_particle :: proc(super_particles: ^[dynamic]Super_Particle, ps: ^[dynamic]Particle, center_pos: [2]f64) {
    particle_template : Particle = {density = density_iron,}

    for _ in 0..<3 {
        radius_rand := rand.float64()
        particle_template.radius = 10 + (10 * radius_rand)
        particle_template.radius_visual = 5 * (1 + radius_rand)
        particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density

        append(ps, particle_template)
    }

    last := len(ps) - 1

    distance :: 3
    a := ps[last].radius_visual + ps[last-1].radius_visual + distance
    b := ps[last-1].radius_visual + ps[last-2].radius_visual + distance
    c := ps[last].radius_visual + ps[last-2].radius_visual + distance

    ps[last].position = center_pos - {a - ps[last-1].radius_visual, 0}
    ps[last-1].position = center_pos + {a - ps[last].radius_visual, 0}

    theta := math.acos((c*c - a*a - b*b) / (-2 * a * b))
    height := math.sin(theta) * c
    ps[last-2].position = center_pos - {0, height}

    ps[last].position_old = ps[last].position
    ps[last-1].position_old = ps[last-1].position
    ps[last-2].position_old = ps[last-2].position

    x : Super_Particle
    x[0] = Link {p = last, p1 = last-1, length = linalg.length(ps[last].position - ps[last-1].position)}
    x[1] = Link {p = last-1, p1 = last-2, length = linalg.length(ps[last-1].position - ps[last-2].position)}
    x[2] = Link {p = last, p1 = last-2, length = linalg.length(ps[last].position - ps[last-2].position)}
    append(super_particles, x)
}

main :: proc() {
    rl.InitWindow(1000, 1000, "particle sim");

    rl.SetTargetFPS(rl.GetMonitorRefreshRate(rl.GetCurrentMonitor()));
    FPS := cast(f64)rl.GetMonitorRefreshRate(rl.GetCurrentMonitor());

    walls : [dynamic]Wall;
    append(&walls,
           Wall {start = {200, 30},end = {97, 183}, callback = null_callback, no_clip = false},
           Wall {start = {97, 183},end = {97, 717}, callback = null_callback, no_clip = false},
           Wall {start = {97, 717},end = {191, 969}, callback = null_callback, no_clip = false},
           Wall {start = {191, 969},end = {241, 969}, callback = inlet_callback, no_clip = false}, // inlet 
           Wall {start = {241, 969},end = {513, 202}, callback = null_callback, no_clip = false},
           Wall {start = {513, 202},end = {484, 30}, callback = null_callback, no_clip = false},
           Wall {start = {484, 30},end = {434, 30}, callback = right_outlet_callback, no_clip = true}, // right outlet
           Wall {start = {434, 30},end = {400, 227}, callback = null_callback, no_clip = false},
           Wall {start = {400, 227},end = {344, 227}, callback = null_callback, no_clip = false},
           Wall {start = {344, 227},end = {344, 170}, callback = null_callback, no_clip = false},
           Wall {start = {344, 170},end = {250, 30}, callback = null_callback, no_clip = false},
           Wall {start = {250, 30},end = {200, 30}, callback = left_outlet_callback, no_clip = true}, // left outlet
          )

    force_points : [dynamic]Force_Point
    fill_geometry_with_points(&force_points, walls[:])


    max_particles :: 100
    ps : [dynamic]Particle
    super_particles : [dynamic]Super_Particle

    particle_template : Particle = {
        density = density_iron,
        radius = 100, // nm
        radius_visual = 5,
    }

    for _ in 0..<max_particles {
        make_super_particle(&super_particles, &ps, [2]f64{rand.float64()*416 + 97, rand.float64()*939 + 30})
    }

    // make_super_particle(&super_particles, &ps, [2]f64{0.5*416 + 97, 0.5*939 + 30})

    radius_rand := rand.float64()
    particle_template.radius = 10 + (10 * radius_rand)
    particle_template.radius_visual = 5 * (1 + radius_rand)
    particle_template.mass = sphere_volume(particle_template.radius) * particle_template.density
    particle_template.position = [2]f64{0.5*370 + 97, 0.5*1100 + 30}
    particle_template.position_old = particle_template.position - {0,-1}
    append(&ps, particle_template)


    held_mouse_pos : [2]f64;

    mouse_particle_idx := -1

    magnet: [2]f64 = {700, 100};
    radius_magnet :: 20
    magnet_selected := false

    
    for !rl.WindowShouldClose() {
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

        for &p, p_idx in ps {
            if p.disabled do continue
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
            if ps[i].disabled do continue
            
            for j in i+1..<len(ps) {
                if ps[j].disabled do continue
                resolve_collision(ps[:], i, j)
            }

            for j in 0..<len(walls) {
                check_particle_wall_collision(&ps[i], walls[j]);
            }
        }

        for &p in ps {
            if p.disabled do continue
            add_forces(&p, magnet)
            update(&p, dt)
        }

        for &links in super_particles {
            for link in links do update_link(ps[:], link)
            set_angular_position(ps[:], &links)
        }
        
        rl.BeginDrawing();
        rl.ClearBackground(rl.RAYWHITE);
        
        
        draw_points(force_points[:])


        for wall in walls {
            draw_wall(wall);
        }

        if mouse_particle_idx >= len(ps) do mouse_particle_idx = -1
        if (mouse_particle_idx >= 0) {
            p := ps[mouse_particle_idx]
            v := (p.position - p.position_old) / dt
            a : = p.force / p.mass / NM_PER_PX
            w := (p.angular_position - p.angular_position_old)/dt
            rl.DrawText(rl.TextFormat("Acceleration: (%f, %f) px/s^2", a.x, a.y), 14, 14, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Acceleration: %f px/s^2", linalg.length(a)), 14, 14+20+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Velocity: (%f, %f) px/s", v.x, v.y), 14, 14+40+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Velocity: %f px/s", linalg.length(v)), 14, 14+60+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Angular Velocity: %f rad/s", w), 14, 14+80+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("Position: (%f, %f) px", p.position.x, p.position.y), 14, 14+100+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("mass: %.10f ng", p.mass), 14, 14+120+2, 20, rl.BLACK);
            rl.DrawText(rl.TextFormat("radius: %f nm\tradius_visual: %f px", p.radius, p.radius_visual), 14, 14+140+2, 20, rl.BLACK);
        }

        rl.DrawText(rl.TextFormat("(%f, %f)", mousePos32.x, mousePos32.y), 10, 1000 - 22, 20, rl.BLACK)

        rl.DrawText(rl.TextFormat("left outlet: %d\nRight outlet: %d", left_outlet_count, right_outlet_count), 800, 10, 20, rl.BLACK)

        for p, p_idx in ps {

            pos : [2]f32 = {cast(f32)p.position.x, cast(f32)p.position.y}
            pos_old : [2]f32 = {cast(f32)p.position_old.x, cast(f32)p.position_old.y}
            v := (pos - pos_old) / cast(f32)dt
            rl.DrawCircleV(pos, cast(f32)p.radius_visual, rl.BLACK);
            if p.disabled {
                rl.DrawCircleV(pos, cast(f32)p.radius_visual + 1, rl.RED);
            }
            // rl.DrawLineV(pos, pos + v*10, rl.BLUE)
            // if p_idx == mouse_particle_idx do rl.DrawCircleLinesV(pos, radius_visual+1, rl.RED);
        }

        for &links in super_particles {
            for &link in links {
                pos := ps[link.p].position
                pos1 := ps[link.p1].position
                rl.DrawLineV(cast_vec2_f32(pos), cast_vec2_f32(pos1), rl.BLACK)
            }
        }

        mag := [2]f32{cast(f32)magnet.x, cast(f32)magnet.y}
        rl.DrawCircleV(mag, radius_magnet, rl.BLUE);
        if magnet_selected do rl.DrawCircleLinesV(mag, radius_magnet+1, rl.RED);

        rl.EndDrawing();
    }

    rl.CloseWindow();
}

cast_vec2_f32 :: proc(v: [2]f64) -> [2]f32 {
    return {cast(f32)v.x, cast(f32)v.y}
}
