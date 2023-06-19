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
    let mut solver = FastSolver::new(voxels, &voxel_endpoints);
    solver.solve_from(0);
    println!("[{:?}] Finish pathfinding", Instant::now().duration_since(tstart));

    let pipe_mesh = solver.voxels.to_mesh(dim_u, 2);
    pipe_mesh.save("pipes.stl");
    println!("[{:?}] Generated output pipe mesh", Instant::now().duration_since(tstart));

    // schedule of materials
    let mut voxels = solver.voxels;

    let n_pipes = voxels.voxels.iter().filter(|x| **x == 2).count();
    let pipe_m = n_pipes as f32 * VOXEL_S_M;
    println!("Schedule of materials:");
    println!("Total pipe length: {}m", pipe_m);

    // filling in endpoints with pipe so they count as elbows if they immediately turn coming out
    for vox_endpoint in voxel_endpoints.iter() {
        let posi = (vox_endpoint.x as isize + vox_endpoint.nx as isize, vox_endpoint.y as isize + vox_endpoint.ny as isize, vox_endpoint.z as isize + vox_endpoint.nz as isize);
        let endpoint_idx = voxels.get_idx_unchecked_i(posi);
        voxels.voxels[endpoint_idx] = 2;
    }
    
    let elbow_inds = voxels.elbow_inds();
    println!("Elbows: {}", elbow_inds.len());
    for idx in elbow_inds {
        voxels.voxels[idx] = 4;
        let elbow_mesh = voxels.to_mesh(dim_u, 4);
        elbow_mesh.save("elbows.stl")
    }
    
    let t_inds = voxels.t_inds();
    println!("T junctions: {}", t_inds.len());
    // i think it fails to detect the singular T junction, might be because its against an endpoint

    // ran out of time to do individual pipes but I would do it in a similar way to how I assigned the vertices to their corresponding enpoint

    println!("[{:?}] Finished schedule of materials", Instant::now().duration_since(tstart));
}
