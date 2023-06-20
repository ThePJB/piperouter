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
        
        // Convert from their types to our types
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

        // Find indices of slanty triangles
        // Slanty Triangles are triangles with normals which do not point in a cardinal direction, this formula will discern them
        let slanty_triangles_indices: Vec<usize> = tris.iter().enumerate()
            .filter(|(idx, tri)| tri.normal * tri.normal != tri.normal.abs())
            .map(|(idx, tri)| idx)
            .collect();

        // Here we compute this table of information about the endpoint vertices (width = number of endpoint vertices):
        // endpoint_vert_index: Their index...
        // endpoint_vert_which: Which endpoint they belong to...
        // endpoint_vert_front: Whether they are a front vertex (around the hole) or not (around the back part)

        // Endpoint vertices are all vertices referred to by the Slanty Triangles
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
            // This algorithm combines the 'which' values of vertices until all vertices belonging to the same endpoint have the same 'which' value
            let mut endpoint_vert_which: Vec<usize> = (0..endpoint_vert_index.len()).collect();
            // go through triangles and if a triangle has multiple endpoint vertices then the endpoints will be reconciled
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

        let endpoint_vert_front = {
            let mut endpoint_vert_front = vec![true; endpoint_vert_index.len()];

            // so for each NON-SLANTY triangle, if all the verts are endpoint verts, those verts are not front
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

        // get a list of the endpoint-keys that the endpoints ended up having (not just 0,1,2.. more like 0,5,9)
        let keys = {
            let mut keys = endpoint_vert_which.clone();
            keys.sort();
            keys.dedup();
            keys
        };

        // the very same table but as a vec of triples for reasons
        let endpoint_vert_records: Vec<(usize, usize, bool)> = (0..endpoint_vert_index.len()).map(|i| (endpoint_vert_index[i], endpoint_vert_which[i], endpoint_vert_front[i])).collect();
            
        // return value
        let mut endpoints = Vec::new();

        // for each distinct endpoint
        for key in keys {
            // avg of front vertices is the center of the front hole
            let avg_front = endpoint_vert_records.iter()
                .filter(|(idx, which, front)| *which == key && *front == true)
                .map(|(idx, _, _)| verts[*idx])
                .fold(v3(0.0, 0.0, 0.0), |acc, v| acc + v) / 5.0;

            // avg of back vertices is center of back face
            let avg_back = endpoint_vert_records.iter()
                .filter(|(idx, which, front)| *which == key && *front == false)
                .map(|(idx, _, _)| verts[*idx])
                .fold(v3(0.0, 0.0, 0.0), |acc, v| acc + v) / 5.0;

            // normal for endpoint given by (front - back) normalized
            endpoints.push(Endpoint {
                pos: avg_front,
                normal: (avg_front - avg_back).norm(),
            })
        }
        endpoints
    }

    // get span of geometry EXCLUDING endpoint markers
    pub fn get_dims(&self) -> V3 {
        let mut dim_max = v3(-f32::INFINITY, -f32::INFINITY, -f32::INFINITY);
        let mut dim_min = v3(f32::INFINITY, f32::INFINITY, f32::INFINITY);
        
        
        // Find indices of slanty triangles
        let slanty_triangles_indices: Vec<usize> = self.tris.iter().enumerate()
                .filter(|(idx, tri)| tri.normal * tri.normal != tri.normal.abs())
                .map(|(idx, tri)| idx)
                .collect();
        
        // Set of vertices to ignore in bounds calculation
        let mut endpoint_vert_index_set = HashSet::new();
        for i in slanty_triangles_indices.iter().copied() {
            endpoint_vert_index_set.insert(self.tris[i].i1);
            endpoint_vert_index_set.insert(self.tris[i].i2);
            endpoint_vert_index_set.insert(self.tris[i].i3);
        }

        for (i, vert) in self.verts.iter().copied().enumerate() {
            if !endpoint_vert_index_set.contains(&i) {
                dim_max = dim_max.max(vert);
                dim_min = dim_min.min(vert);
            }
        }
        
        dim_max - dim_min
    }
}