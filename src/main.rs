mod math;
mod fast_solver;
mod voxel;
mod mesh;

use math::*;
use voxel::*;

use crate::fast_solver::*;
use crate::mesh::*;
use crate::mesh::IndexedMesh;

fn main() {
    std::env::set_var("RUST_BACKTRACE", "1");

    let mesh = IndexedMesh::from_file("Volume.stl");

    let endpoints = mesh.find_endpoints();
    // dbg!(&endpoints);

    // ok so we up to here. goddamn find endpoints is working
    // weve got the mesh
    // now we need to try the thing

    // dimensions in meters
    let dim_u = v3(1.75, 2.5, 1.5);

    let voxels = Voxels::from_mesh(&mesh, dim_u);
    let dim_vox = voxels.dim;
    let voxels = voxels.voxels;

    // how do I check the correctness of the voxels lmao. export to minecraft.

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

    let voxel_mesh = gen_mesh(&voxels, dim_vox, dim_u, 1);
    voxel_mesh.save("voxmesh.stl");

    let mut solver = FastSolver::new(dim_vox, voxel_endpoints, voxels);
    dbg!("begin");
    solver.solve_from(0);
    dbg!("done");

    let npipe = solver.voxels.iter().filter(|x| **x == 2).count();
    let tot = solver.voxels.len();
    dbg!(npipe, tot, npipe as f32 / tot as f32);

    // I think pheromone needs to be per destination and directional
    // and ants need source and target
    // maybe it cant work in 3d because scourge of dimensionality, like the chance of them actually hitting the target is too low
    // or that it goes into suboptimal solution too quickly

    let pipe_mesh = gen_mesh(&solver.voxels, solver.dim, dim_u, 2);
    pipe_mesh.save("pipes.stl");
    let combined_mesh = mesh.combine(&pipe_mesh);
    combined_mesh.save("combined.stl");
}


// hmm it nearly correct
// the discrepancy might be from the mesh-> voxels step
// I could output a voxel mesh, probably go for a low res one tho otherwise itll be utterly cooked
// and probably should neighbour check to omit quads as well
// maybe just z is flipped or something
// but flip z fails at something 1
// edges n shiz