use crate::math::*;
use std::collections::HashSet;

#[derive(Debug)]
pub struct Endpoint {
    pub pos: V3,
    pub normal: V3,
}

pub struct VoxelEndpoint {
    pub x: usize,
    pub y: usize,
    pub z: usize,
    pub nx: isize,
    pub ny: isize,
    pub nz: isize,
}

#[derive(Debug, Clone, Copy)]
pub struct Tri {
    pub normal: V3,
    pub i1: usize,
    pub i2: usize,
    pub i3: usize,
}

// computes the endpoint positions and normals from an indexed mesh
pub fn find_endpoints(verts: &[V3], tris: &[Tri]) -> Vec<Endpoint> {

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