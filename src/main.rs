mod math;
mod ant_solver;
mod djikstra_solver;
mod fast_solver;
mod endpoint;
mod priority_queue;

use math::*;
use std::{fs::OpenOptions, process};
use stl_io::*;
use ordered_float::*;

use crate::djikstra_solver::*;
use crate::endpoint::*;
use crate::ant_solver::*;

// anyway so maybe it wouldn't be so hard to do a basic pathfinding that explored in straight lines first
// so its ranking search by generation
// could ant shit work? maybe


// maybe theres potential for a heuristic by raycasting on the geometry


fn main() {
    std::env::set_var("RUST_BACKTRACE", "1");
    let mut file = OpenOptions::new().read(true).open("Volume.stl").unwrap();
    let stl = read_stl(&mut file).unwrap();
    
    let tris: Vec<Tri> = stl.faces.iter().map(|face| Tri {
        normal: v3(face.normal[0], face.normal[1], face.normal[2]),
        i1: face.vertices[0],
        i2: face.vertices[1],
        i3: face.vertices[2],
    }).collect();

    let verts: Vec<V3> = stl.vertices.iter().map(|vert| v3(vert[0], vert[1], vert[2])).collect();


    let endpoints = find_endpoints(&verts, &tris);
    // dbg!(&endpoints);

    // ok so we up to here. goddamn find endpoints is working
    // weve got the mesh
    // now we need to try the thing

    // dimensions in meters
    let U_TO_M = 3.6429696;
    let dim_u = v3(1.75, 2.5, 1.5);
    let dim_m = dim_u * U_TO_M;
    let voxel_s = 0.01;
    let voxel_su = voxel_s / U_TO_M;
    let dim_vox = ((dim_m.x / voxel_s) as usize, (dim_m.y / voxel_s) as usize, (dim_m.z / voxel_s) as usize);

    dbg!(&dim_vox);

    let mut voxels = vec![0; dim_vox.0*dim_vox.1*dim_vox.2];

    let z_triangles: Vec<Tri> = tris.iter().filter(|tri| tri.normal.dot(v3(0.0, 0.0, 1.0)) != 0.0).map(|x| *x).collect();

    for i in 0..dim_vox.0 {
        for j in 0..dim_vox.1 {
            // ray intersections down column
            let ray_z_offset = 0.1;
            let ray_origin = v3((i as f32 + 0.5) / dim_vox.0 as f32 * dim_u.x, (j as f32 + 0.5) / dim_vox.1 as f32 * dim_u.y, dim_u.z + ray_z_offset);
            // there was a copy paste error here but i thought it was working
            // now its fixed and saying one of endpoints is at y=1700
            let ray_dir = v3(0.0, 0.0, -1.0);
            let mut intersections = Vec::new();
            for tri in z_triangles.iter() {
                let v0 = verts[tri.i1];
                let v1 = verts[tri.i2];
                let v2 = verts[tri.i3];
                if let Some(d_inter) = ray_triangle_intersection(ray_origin, ray_dir, v0, v1, v2) {
                    intersections.push(OrderedFloat(d_inter - ray_z_offset));
                }
            }
            intersections.sort();

            if i == 0 && j == 0 {
                dbg!(intersections.clone());
            }

            // now make the voxels themselves
            // project from the top and if passing through an odd number of triangles, empty otherwise filled. (Normally it would be the reverse but this volume denotes empty space)
            // hoping the ceiling case is handled
            let mut n_intersections = 0;
            for k in 0..dim_vox.2 {
                let vx_d = k as f32 * voxel_su;
                loop {
                    if n_intersections >= intersections.len() {
                        break;
                    }
                    if vx_d >= intersections[n_intersections].0 {
                        n_intersections += 1;
                    } else {
                        break;
                    }
                }
                if n_intersections % 2 == 1 {
                    voxels[i*dim_vox.1*dim_vox.2 + j*dim_vox.2 + k] = 0;
                } else {
                    voxels[i*dim_vox.1*dim_vox.2 + j*dim_vox.2 + k] = 1;
                }
            }
        }
    }

    // how do I check the correctness of the voxels lmao. export to minecraft.

    println!("done");
    let nfull = voxels.iter().filter(|v| **v == 1).count();
    let nempty = voxels.iter().filter(|v| **v == 0).count();

    println!("nfull {} nempty {} ratio {}", nfull, nempty, nfull as f32 / (nfull + nempty) as f32);

    // well thats what I'd expect so I guess its working

    // what if we add slightly along normal to coord
    // doesnt really matter anyway
    // let hole_voxel_coords: Vec<(usize, usize, usize)> = endpoints.iter().map(|x| x.pos + 0.001*x.normal).map(|v| (((v.x / dim_u.x) * vx_x as f32) as usize, ((v.y / dim_u.y) * vx_y as f32) as usize, ((v.z / dim_u.z) * vx_z as f32) as usize)).collect();


    // lets print out the endpoints
    dbg!(&dim_u);
    dbg!(&endpoints);

    // this could be wrong
    let voxel_endpoints: Vec<VoxelEndpoint> = endpoints.iter().map(|endpoint| {
        VoxelEndpoint {
            x: ((endpoint.pos.x + 0.001*endpoint.normal.x) * dim_vox.0 as f32 / dim_u.x) as usize,
            y: ((endpoint.pos.y + 0.001*endpoint.normal.y) * dim_vox.1 as f32 / dim_u.y) as usize,
            z: ((endpoint.pos.z + 0.001*endpoint.normal.z) * dim_vox.2 as f32 / dim_u.z) as usize,
            nx: endpoint.normal.x.round() as isize,
            ny: endpoint.normal.y.round() as isize,
            nz: endpoint.normal.z.round() as isize,
        }
    }).collect();





    // dbg!(hole_voxel_coords);

    for endpoint in voxel_endpoints.iter() {
        let x = endpoint.x;
        let y = endpoint.y;
        let z = endpoint.z;
        let val = voxels[x * dim_vox.1*dim_vox.2 + y*dim_vox.2 + z];
        dbg!(x, y, z, val);
    }

    // time for god damn pathfinding lol.
    /*
        Using:
            ants?
            cellular?
            vector field?
            generational straight line pathfinder? like greedy bfs but straight lines, i actually like dis 1. Points are prioritized by being after fewer turns
                points in queue need to be unique
     */

    // ok so we also need a 

    // i think lazy pheromone required

    // let mut ant_solver = AntSolver::new(dim_vox, voxel_endpoints, voxels);
    // for i in 0..ant_solver.endpoints.len() {
    //     for n in 0..100000 {
    //         ant_solver.spawn_ant(i);
    //     }
    // }
    // // seems to degenerate like dings go down
    // // also only ding 2
    // dbg!("begin");
    // for i in 0..100000 {
    //     ant_solver.step();
    // }
    // dbg!("done");

    let mut solver = DjikstraSolver::new(dim_vox, voxel_endpoints, voxels);
    dbg!("begin");
    solver.solve_from(0);
    dbg!("done");


    // I think pheromone needs to be per destination and directional
    // and ants need source and target
    // maybe it cant work in 3d because scourge of dimensionality, like the chance of them actually hitting the target is too low
    // or that it goes into suboptimal solution too quickly

}
