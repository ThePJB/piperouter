use crate::mesh::*;
use crate::math::*;
use ordered_float::*;

pub struct VoxelEndpoint {
    pub x: usize,
    pub y: usize,
    pub z: usize,
    pub nx: isize,
    pub ny: isize,
    pub nz: isize,
}

pub const U_TO_M: f32 = 3.6429696;
pub const VOXEL_S_M: f32 = 0.1;
pub const VOXEL_S_U: f32 = VOXEL_S_M / U_TO_M;

pub fn vox_ind(pos: (usize, usize, usize), dim: (usize, usize, usize)) -> usize {
    pos.0*dim.2*dim.1 + pos.1*dim.2 + pos.2
}

pub fn vox_indi(pos: (isize, isize, isize), dim: (usize, usize, usize)) -> usize {
    vox_ind((pos.0 as usize, pos.1 as usize, pos.2 as usize), dim)
}
pub fn vox_ind_ok(pos: (isize, isize, isize), dim: (usize, usize, usize)) -> bool {
    pos.0 >= 0 && pos.1 >= 0 && pos.2 >= 0 &&
    (pos.0 as usize) < dim.0 && (pos.1 as usize) < dim.1 && (pos.2 as usize) < dim.2
}

pub struct Voxels {
    pub voxels: Vec<usize>,
    pub dim: (usize, usize, usize),
}

impl Voxels {
    pub fn from_mesh(mesh: &IndexedMesh, dim_u: V3) -> Self {
        let dim_m = dim_u * U_TO_M;
        let dim_vox: (usize, usize, usize) = ((dim_m.x / VOXEL_S_M) as usize, (dim_m.y / VOXEL_S_M) as usize, (dim_m.z / VOXEL_S_M) as usize);

        let mut voxels = vec![0; dim_vox.0*dim_vox.1*dim_vox.2];

        let z_triangles: Vec<Tri> = mesh.tris.iter().filter(|tri| tri.normal.dot(v3(0.0, 0.0, 1.0)) != 0.0).map(|x| *x).collect();
    
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
                    let v0 = mesh.verts[tri.i1];
                    let v1 = mesh.verts[tri.i2];
                    let v2 = mesh.verts[tri.i3];
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
                    let vx_d = k as f32 * VOXEL_S_U;
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
                    let k_flip = dim_vox.2 - k - 1;
                    if n_intersections % 2 == 1 {
                        voxels[i*dim_vox.1*dim_vox.2 + j*dim_vox.2 + k_flip] = 0;
                    } else {
                        voxels[i*dim_vox.1*dim_vox.2 + j*dim_vox.2 + k_flip] = 1;
                    }
                }
            }
        }

        Voxels { voxels, dim: dim_vox }
    }
}