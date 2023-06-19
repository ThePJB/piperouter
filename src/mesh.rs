use crate::math::*;
use crate::voxel::*;
use std::fs::OpenOptions;
use std::collections::HashSet;


#[derive(Debug)]
pub struct Endpoint {
    pub pos: V3,
    pub normal: V3,
}

#[derive(Debug, Clone, Copy)]
pub struct Tri {
    pub normal: V3,
    pub i1: usize,
    pub i2: usize,
    pub i3: usize,
}

#[derive(Clone)]
pub struct IndexedMesh {
    pub tris: Vec<Tri>,
    pub verts: Vec<V3>,
}

impl IndexedMesh {
    pub fn from_file(path: &str) -> Self {
        let mut file = OpenOptions::new().read(true).open(path).unwrap();
        let stl = stl_io::read_stl(&mut file).unwrap();
        
        let tris: Vec<Tri> = stl.faces.iter().map(|face| Tri {
            normal: v3(face.normal[0], face.normal[1], face.normal[2]),
            i1: face.vertices[0],
            i2: face.vertices[1],
            i3: face.vertices[2],
        }).collect();
    
        let verts: Vec<V3> = stl.vertices.iter().map(|vert| v3(vert[0], vert[1], vert[2])).collect();
    
        IndexedMesh {tris, verts}
    }

    pub fn save(&self, path: &str) {
        let stl_tris = self.tris.iter().map(|tri| {
            let v1 = self.verts[tri.i1];
            let v2 = self.verts[tri.i2];
            let v3 = self.verts[tri.i3];

            stl_io::Triangle {
                normal: stl_io::Normal::new([tri.normal.x, tri.normal.y, tri.normal.z]),
                vertices: [
                    stl_io::Vertex::new([v1.x, v1.y, v1.z]),
                    stl_io::Vertex::new([v2.x, v2.y, v2.z]),
                    stl_io::Vertex::new([v3.x, v3.y, v3.z]),
                ],
            }
        });

        let mut binary_stl = Vec::<u8>::new();
        stl_io::write_stl(&mut binary_stl, stl_tris).unwrap();
        std::fs::write(path, binary_stl).unwrap();

    }

    pub fn combine(&self, other: &Self) -> Self {
        let mut new_mesh = self.clone();
        let idx_offset = new_mesh.verts.len();
        new_mesh.verts.extend(other.verts.iter());
        new_mesh.tris.extend(other.tris.iter().map(|tri| {
            let mut new_tri = tri.clone();
            new_tri.i1 += idx_offset;
            new_tri.i2 += idx_offset;
            new_tri.i3 += idx_offset;
            new_tri
        }));
        new_mesh
    }

    // computes the endpoint positions and normals from an indexed mesh
    pub fn find_endpoints(&self) -> Vec<Endpoint> {
        let verts = &self.verts;
        let tris = &self.tris;

        // Step 1: Determine endpoints
        // Triangles of the endpoints other than the endcap (slanty triangles) will have non cardinally aligned normals. Test for this is N^2 != |N|
        // These triangles will index either vertices of the endcap or vertices of the hole we wish to plug in the geometry. The centroid of the hole will yield the final location.
        // We collect the indexes on the endpoint triangles.
        // We wish to partition these indexes into endcap and hole. Such that we may delete endcap and plug hole. After deleting slanty triangles, there will exist triangles between endcap indices while there will not exist triangles between hole indices.

        // Find indices of slanty triangles
        let slanty_triangles_indices: Vec<usize> = tris.iter().enumerate()
            .filter(|(idx, tri)| tri.normal * tri.normal != tri.normal.abs())
            .map(|(idx, tri)| idx)
            .collect();

        // dbg!(slanty_triangles_indices);

        let endpoint_vert_index = {
            let mut endpoint_vert_index_set = HashSet::new();
            for i in slanty_triangles_indices.iter() {
                endpoint_vert_index_set.insert(tris[*i].i1);
                endpoint_vert_index_set.insert(tris[*i].i2);
                endpoint_vert_index_set.insert(tris[*i].i3);
            }
            let mut endpoint_vert_index: Vec<usize> = endpoint_vert_index_set.iter().copied().collect();
            endpoint_vert_index.sort();
            endpoint_vert_index
        };
        let endpoint_vert_which = {
            let mut endpoint_vert_which: Vec<usize> = (0..endpoint_vert_index.len()).collect();
            // go through triangles and if a triangle has multiple hole vertices then the holes will be reconciled
            for idx in slanty_triangles_indices.iter().copied() {
                let tri = tris[idx];
                // find positions of referent vertices in endpoints vec
                let p1 = endpoint_vert_index.iter().position(|x| *x == tri.i1).unwrap();
                let p2 = endpoint_vert_index.iter().position(|x| *x == tri.i2).unwrap();
                let p3 = endpoint_vert_index.iter().position(|x| *x == tri.i3).unwrap();
        
                // find keys
                let k1 = endpoint_vert_which[p1];
                let k2 = endpoint_vert_which[p2];
                let k3 = endpoint_vert_which[p3];
        
                // merge their keys to the lowest key
                let lowest_key = k1.min(k2.min(k3));
                // any k1 k2 or k3 replaced with lowest
                for i in 0..endpoint_vert_which.len() {
                    if endpoint_vert_which[i] == k1 || endpoint_vert_which[i] == k2 || endpoint_vert_which[i] == k3 {
                        endpoint_vert_which[i] = lowest_key;
                    }
                }
            }
            endpoint_vert_which
        };

        // dbg!(&endpoint_vert_index.len());
        // dbg!(&endpoint_vert_index);
        // dbg!(&endpoint_vert_which.len());
        // dbg!(&endpoint_vert_which);
        // std::process::exit(0);

        let endpoint_vert_front = {
            let mut endpoint_vert_front = vec![true; endpoint_vert_index.len()];
            // you're the front, but if youre touched by a ENDCAP TRIANGLE you get set to back
            // ENDCAP TRIANGLES - triangles with all 3 vertices from the endpoint group

            // so for each triangle, if all the points are contained within endpoint_vert_index, those points are not front(by index of their occurrence in the table)
            // oh because its for each non slanty triangle

            // so it be that the triangles are all ebing found false

            // y dis wrong

            // for over non slanty triangles
            for tri in tris.iter().filter(|tri| tri.normal*tri.normal == tri.normal.abs()) {
                if let Some(idx1) = endpoint_vert_index.iter().position(|x| *x == tri.i1) {
                    if let Some(idx2) = endpoint_vert_index.iter().position(|x| *x == tri.i2) {
                        if let Some(idx3) = endpoint_vert_index.iter().position(|x| *x == tri.i3) {
                            endpoint_vert_front[idx1] = false;
                            endpoint_vert_front[idx2] = false;
                            endpoint_vert_front[idx3] = false;
                        };
                    };
                };
            }
            endpoint_vert_front
        };

        // dbg!(&endpoint_vert_front.len());
        // dbg!(&endpoint_vert_front);
        // std::process::exit(0);

        // now we would go by key
        // and calculate resulting Vec<Endpoint> for return
        let mut endpoints = Vec::new();

        let keys = {
            let mut keys = endpoint_vert_which.clone();
            keys.sort();
            keys.dedup();
            keys
        };

        let endpoint_vert_records: Vec<(usize, usize, bool)> = (0..endpoint_vert_index.len()).map(|i| (endpoint_vert_index[i], endpoint_vert_which[i], endpoint_vert_front[i])).collect();
        // dbg!(&endpoint_vert_records.len());
        // dbg!(&endpoint_vert_records);
        // std::process::exit(0);
            
        // for each distinct endpoint
        for key in keys {
            // avg of front endpoints will be fold over vertex positions of this specific key
            let avg_front = endpoint_vert_records.iter()
                .filter(|(idx, which, front)| *which == key && *front == true)
                .map(|(idx, _, _)| verts[*idx])
                .fold(v3(0.0, 0.0, 0.0), |acc, v| acc + v) / 5.0;

            let avg_back = endpoint_vert_records.iter()
                .filter(|(idx, which, front)| *which == key && *front == false)
                .map(|(idx, _, _)| verts[*idx])
                .fold(v3(0.0, 0.0, 0.0), |acc, v| acc + v) / 5.0;

            endpoints.push(Endpoint {
                pos: avg_front,
                normal: (avg_front - avg_back).norm(),
            })
        }
        endpoints
    }
}

// Makes a mesh of voxels denoted by voxel_value
// e.g. 0 would be all air, 1 would be all walls, 2 would be all pipes
// rudimentary mesh optimization only (occluded quads not added)
pub fn gen_mesh(voxels: &[usize], dim: (usize, usize, usize), dim_u: V3, voxel_value: usize) -> IndexedMesh {
    let mut verts = Vec::new();
    let mut tris = Vec::new();

    let vi = v3(VOXEL_S_U, 0.0, 0.0);
    let vj = v3(0.0, VOXEL_S_U, 0.0);
    let vk = v3(0.0, 0.0, VOXEL_S_U);

    for i in 0..dim.0 {
        for j in 0..dim.1 {
            for k in 0..dim.2 {
                let idx = vox_ind((i,j,k), dim);
                if voxels[idx] == voxel_value {
                    // push verts and tris
                    let corner = v3(i as f32 / dim.0 as f32 * dim_u.x, j as f32 / dim.1 as f32 * dim_u.y, k as f32 / dim.2 as f32 * dim_u.z);
                    verts.push(corner); // -8
                    verts.push(corner + vi);
                    verts.push(corner + vj);
                    verts.push(corner + vi + vj);
                    verts.push(vk + corner);
                    verts.push(vk + corner + vi);
                    verts.push(vk + corner + vj);   // -2
                    verts.push(vk + corner + vi + vj); // -1
                    let mi = verts.len();
                    // -z quad
                    if !vox_ind_ok((i as isize,j as isize,k as isize-1), dim) || voxels[vox_ind((i,j,k-1), dim)] != voxel_value {
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
                    if !vox_ind_ok((i as isize,j as isize,k as isize+1), dim) || voxels[vox_ind((i,j,k+1), dim)] != voxel_value {
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
                    if !vox_ind_ok((i as isize - 1,j as isize,k as isize), dim) || voxels[vox_ind((i-1,j,k), dim)] != voxel_value {
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
                    if !vox_ind_ok((i as isize + 1,j as isize,k as isize), dim) || voxels[vox_ind((i+1,j,k), dim)] != voxel_value {
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
                    if !vox_ind_ok((i as isize,j as isize - 1,k as isize), dim) || voxels[vox_ind((i,j - 1,k), dim)] != voxel_value {
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
                    if !vox_ind_ok((i as isize,j as isize + 1,k as isize), dim) || voxels[vox_ind((i,j + 1,k), dim)] != voxel_value {
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