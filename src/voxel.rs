use crate::mesh::*;
use crate::math::*;
use ordered_float::*;

#[derive(Debug)]
pub struct VoxelEndpoint {
    pub x: usize,
    pub y: usize,
    pub z: usize,
    pub nx: i8,
    pub ny: i8,
    pub nz: i8,
}

pub const U_TO_M: f32 = 3.6429696;
pub const VOXEL_S_M: f32 = 0.01;
pub const VOXEL_S_U: f32 = VOXEL_S_M / U_TO_M;

// pub fn vox_ind(pos: (usize, usize, usize), dim: (usize, usize, usize)) -> usize {
//     pos.0*dim.2*dim.1 + pos.1*dim.2 + pos.2
// }

// pub fn vox_indi(pos: (isize, isize, isize), dim: (usize, usize, usize)) -> usize {
//     vox_ind((pos.0 as usize, pos.1 as usize, pos.2 as usize), dim)
// }
// pub fn vox_ind_ok(pos: (isize, isize, isize), dim: (usize, usize, usize)) -> bool {
//     pos.0 >= 0 && pos.1 >= 0 && pos.2 >= 0 &&
//     (pos.0 as usize) < dim.0 && (pos.1 as usize) < dim.1 && (pos.2 as usize) < dim.2
// }

pub struct Voxels {
    pub voxels: Vec<u8>,
    pub dim: (usize, usize, usize),
}

impl Voxels {
    pub fn get_idx_unchecked_u(&self, pos: (usize, usize, usize)) -> usize {
        pos.0*self.dim.2*self.dim.1 + pos.1*self.dim.2 + pos.2
    }
    pub fn get_idx_unchecked_i(&self, pos: (isize, isize, isize)) -> usize {
        self.get_idx_unchecked_u((pos.0 as usize, pos.1 as usize, pos.2 as usize))
    }
    pub fn get_unchecked_i(&self, pos: (isize, isize, isize)) -> u8 {
        let idx = self.get_idx_unchecked_i(pos);
        self.voxels[idx]
    }
    pub fn get_unchecked_u(&self, pos: (usize, usize, usize)) -> u8 {
        let idx = self.get_idx_unchecked_u(pos);
        self.voxels[idx]
    }
    pub fn pos_in_bounds_i(&self, pos: (isize, isize, isize)) -> bool {
        pos.0 >= 0 && pos.1 >= 0 && pos.2 >= 0 &&
        (pos.0 as usize) < self.dim.0 && (pos.1 as usize) < self.dim.1 && (pos.2 as usize) < self.dim.2
    }

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

        
    // Constructs an indexed mesh of voxels denoted by voxel_value
    // e.g. 0 would be all air, 1 would be all walls, 2 would be all pipes
    // rudimentary mesh optimization only (occluded quads not added)
    pub fn to_mesh(&self, dim_u: V3, voxel_value: u8) -> IndexedMesh {
        let dim = self.dim;

        let mut verts = Vec::new();
        let mut tris = Vec::new();

        let vi = v3(VOXEL_S_U, 0.0, 0.0);
        let vj = v3(0.0, VOXEL_S_U, 0.0);
        let vk = v3(0.0, 0.0, VOXEL_S_U);

        for i in 0..dim.0 {
            for j in 0..dim.1 {
                for k in 0..dim.2 {
                    if self.get_unchecked_u((i,j,k)) == voxel_value {
                        // push verts and tris
                        let corner = v3(i as f32 / dim.0 as f32 * dim_u.x, j as f32 / dim.1 as f32 * dim_u.y, k as f32 / dim.2 as f32 * dim_u.z);
                        verts.push(corner);                 // -8
                        verts.push(corner + vi);            // -7
                        verts.push(corner + vj);            // -6
                        verts.push(corner + vi + vj);       // -5
                        verts.push(vk + corner);            // -4
                        verts.push(vk + corner + vi);       // -3
                        verts.push(vk + corner + vj);       // -2
                        verts.push(vk + corner + vi + vj);  // -1
                        let mi = verts.len();
                        // -z quad
                        if !self.pos_in_bounds_i((i as isize, j as isize, k as isize - 1)) || self.get_unchecked_u((i, j, k-1)) != voxel_value {
                            tris.push(Tri {
                                normal: v3(0.0, 0.0, -1.0),
                                i1: mi - 8,
                                i2: mi - 7,
                                i3: mi - 5,
                            });
                            tris.push(Tri {
                                normal: v3(0.0, 0.0, -1.0),
                                i1: mi - 8,
                                i2: mi - 6,
                                i3: mi - 5,
                            });
                        }
                        // +z quad
                        if !self.pos_in_bounds_i((i as isize, j as isize, k as isize + 1)) || self.get_unchecked_u((i, j, k+1)) != voxel_value {
                            tris.push(Tri {
                                normal: v3(0.0, 0.0, 1.0),
                                i1: mi - 4,
                                i2: mi - 3,
                                i3: mi - 1,
                            });
                            tris.push(Tri {
                                normal: v3(0.0, 0.0, 1.0),
                                i1: mi - 4,
                                i2: mi - 2,
                                i3: mi - 1,
                            });
                        }
                        // -x quad
                        if !self.pos_in_bounds_i((i as isize - 1, j as isize, k as isize)) || self.get_unchecked_u((i-1,j,k)) != voxel_value {
                            tris.push(Tri {
                                normal: v3(-1.0, 0.0, 0.0),
                                i1: mi - 8,
                                i2: mi - 6,
                                i3: mi - 4,
                            });
                            tris.push(Tri {
                                normal: v3(-1.0, 0.0, 0.0),
                                i1: mi - 6,
                                i2: mi - 4,
                                i3: mi - 2,
                            });
                        }
                        // +x quad
                        if !self.pos_in_bounds_i((i as isize + 1, j as isize, k as isize)) || self.get_unchecked_u((i+1,j,k)) != voxel_value {
                            tris.push(Tri {
                                normal: v3(1.0, 0.0, 0.0),
                                i1: mi - 7,
                                i2: mi - 5,
                                i3: mi - 3,
                            });
                            tris.push(Tri {
                                normal: v3(1.0, 0.0, 0.0),
                                i1: mi - 5,
                                i2: mi - 3,
                                i3: mi - 1,
                            });
                        }
                        // -y quad
                        if !self.pos_in_bounds_i((i as isize, j as isize - 1, k as isize)) || self.get_unchecked_u((i,j-1,k)) != voxel_value {
                            tris.push(Tri {
                                normal: v3(0.0, -1.0, 0.0),
                                i1: mi - 8,
                                i2: mi - 7,
                                i3: mi - 4,
                            });                    
                            tris.push(Tri {
                                normal: v3(0.0, -1.0, 0.0),
                                i1: mi - 3,
                                i2: mi - 7,
                                i3: mi - 4,
                            });
                        }
                        // +y quad
                        if !self.pos_in_bounds_i((i as isize, j as isize - 1, k as isize)) || self.get_unchecked_u((i,j+1,k)) != voxel_value {
                            tris.push(Tri {
                                normal: v3(0.0, 1.0, 0.0),
                                i1: mi - 6,
                                i2: mi - 5,
                                i3: mi - 2,
                            });                    
                            tris.push(Tri {
                                normal: v3(0.0, 1.0, 0.0),
                                i1: mi - 1,
                                i2: mi - 5,
                                i3: mi - 2,
                            });
                        }
                    }
                }
            }
        }
        IndexedMesh {tris, verts}
    }
}