mod math;
mod fast_solver;
mod voxel;
mod mesh;

use voxel::*;
use std::time::Instant;

use crate::fast_solver::*;
use crate::mesh::IndexedMesh;

fn main() {
    let tstart = Instant::now();
    std::env::set_var("RUST_BACKTRACE", "1");
    let path = "Volume.stl";

    println!("[{:?}] Piperouter starting", Instant::now().duration_since(tstart));
    let mesh = IndexedMesh::from_file(path);
    println!("[{:?}] Loaded mesh {}", Instant::now().duration_since(tstart), path);

    let endpoints = mesh.find_endpoints();
    println!("[{:?}] Computed endpoints {:?}", Instant::now().duration_since(tstart), &endpoints);
    
    let dim_u = mesh.get_dims();
    println!("[{:?}] Computed mesh extent {:?}", Instant::now().duration_since(tstart), dim_u);
    
    let mut voxels = Voxels::from_mesh(&mesh, dim_u);
    println!("[{:?}] Generated voxels from mesh ({} voxels)", Instant::now().duration_since(tstart), voxels.voxels.len());
    
    voxels.dilate();
    println!("[{:?}] Finish dilating walls)", Instant::now().duration_since(tstart));

    let dim_vox = voxels.dim;

    // Convert endpoints in model units to endpoints in voxel coordinates
    let voxel_endpoints: Vec<VoxelEndpoint> = endpoints.iter().map(|endpoint| {
        VoxelEndpoint {
            x: ((endpoint.pos.x + 0.001*endpoint.normal.x) * dim_vox.0 as f32 / dim_u.x) as usize,
            y: ((endpoint.pos.y + 0.001*endpoint.normal.y) * dim_vox.1 as f32 / dim_u.y) as usize,
            z: ((endpoint.pos.z + 0.001*endpoint.normal.z) * dim_vox.2 as f32 / dim_u.z) as usize,
            nx: endpoint.normal.x.round() as i8,
            ny: endpoint.normal.y.round() as i8,
            nz: endpoint.normal.z.round() as i8,
        }
    }).collect();
    println!("[{:?}] Endpoints in voxel coordinates {:?}", Instant::now().duration_since(tstart), voxel_endpoints);
    
    // let voxel_mesh = gen_mesh(&voxels, dim_vox, dim_u, 1);
    // voxel_mesh.save("voxmesh.stl");

    println!("[{:?}] Begin pathfinding", Instant::now().duration_since(tstart));
    let mut solver = FastSolver::new(voxels, voxel_endpoints);
    solver.solve_from(0);
    println!("[{:?}] Finish pathfinding", Instant::now().duration_since(tstart));



    let pipe_mesh = solver.voxels.to_mesh(dim_u, 2);
    pipe_mesh.save("pipes.stl");
    println!("[{:?}] Generated output pipe mesh", Instant::now().duration_since(tstart));

    // voxels to materials
    // first record then delete junctions
    // then calculate the different 4con segments and those are the pipe segments

    // let combined_mesh = mesh.combine(&pipe_mesh);
    // combined_mesh.save("combined.stl");

    // todo voxel mesh dilate
    // todo output 
    // Output a schedule of material, consisting of the count, length, total number
    // of units for each pipe, and connection required
}
