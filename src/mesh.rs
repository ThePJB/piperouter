use crate::math::*;
use crate::voxel::*;
use std::fs::OpenOptions;
use stl_io::*;

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
        let mut stl = read_stl(&mut file).unwrap();
        
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
}


// lets make an indexed mesh of voxels
pub fn gen_pipe_mesh(voxels: Vec<usize>, dim: (usize, usize, usize), dim_u: V3) -> IndexedMesh {
    let mut verts = Vec::new();
    let mut tris = Vec::new();

    let vi = v3(VOXEL_S_U, 0.0, 0.0);
    let vj = v3(0.0, VOXEL_S_U, 0.0);
    let vk = v3(0.0, 0.0, VOXEL_S_U);

    for i in 0..dim.0 {
        for j in 0..dim.1 {
            for k in 0..dim.2 {
                let idx = vox_ind((i,j,k), dim);
                if voxels[idx] == 2 {
                    // push verts and tris
                    // verts are duplicated between adjacent voxels - I don't care
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
                    // bottom quad
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
                    // top quad: +4
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
                    // -x quad
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
                    // +x quad
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
                    // -y quad
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
                    // +y quad
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
    IndexedMesh {tris, verts}
}