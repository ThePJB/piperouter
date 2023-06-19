use crate::voxel::*;
use std::collections::VecDeque;

// maybe it should actually search from longest one first? but might not be shortest
// oh rn im poppin so might be a long path
// if inner one was a vecdeq then

pub struct Entry {
    pos: (isize, isize, isize),
    dir: (isize, isize, isize),
}

pub struct FastSolver {
    pub dim: (usize, usize, usize),

    pub endpoints: Vec<VoxelEndpoint>,
    pub endpoint_inds: Vec<usize>,

    pub voxels: Vec<usize>,
    pub backpointers: Vec<Option<(isize, isize, isize)>>,

    pub queues: Vec<VecDeque<Entry>>,
    pub queue_ind: usize,
}

impl FastSolver {
    pub fn new(dim: (usize, usize, usize), endpoints: Vec<VoxelEndpoint>, voxels: Vec<usize>) -> Self {
        let len = voxels.len();
        let endpoint_inds = endpoints.iter().map(|endpoint| (endpoint.x, endpoint.y, endpoint.z)).map(|pos| vox_ind(pos, dim)).collect();
        FastSolver {
            dim,
            endpoints,
            voxels,
            backpointers: vec![None; len],
            endpoint_inds,
            queues: Vec::new(),
            queue_ind: 0,
        }
    }
    pub fn solve_from(&mut self, start_endpoint_idx: usize) {
        let start_entry = Entry {
            pos: (self.endpoints[start_endpoint_idx].x as isize, self.endpoints[start_endpoint_idx].y as isize, self.endpoints[start_endpoint_idx].z as isize),
            dir: (self.endpoints[start_endpoint_idx].nx, self.endpoints[start_endpoint_idx].ny, self.endpoints[start_endpoint_idx].nz),
        };

        let start_ind = vox_indi(start_entry.pos, self.dim);
        self.queues.push(vec![start_entry].into());
        self.backpointers[start_ind] = Some((0, 0, 0));

        loop {
            if self.queue_ind >= self.queues.len() {
                dbg!("finishing -- exhausted all queues", self.queue_ind);
                return;
            }
            if self.endpoint_inds.iter().map(|x| self.backpointers[*x]).all(|x| x.is_some()) {
                // dbg!("finished -- success", self.queue_ind);
                break;
            }
            if let Some(entry) = self.queues[self.queue_ind].pop_front() {
                let ind = vox_indi(entry.pos, self.dim);
                assert!(self.voxels[ind] != 1);
                // check and push front
                // not oob, wall, or prev
                let fwd_pos = (entry.pos.0 + entry.dir.0, entry.pos.1 + entry.dir.1, entry.pos.2 + entry.dir.2);
                if vox_ind_ok(fwd_pos, self.dim) {
                    let fwd_ind = vox_indi(fwd_pos, self.dim);
                    if self.voxels[fwd_ind] == 0 && self.backpointers[fwd_ind].is_none() {
                        self.queues[self.queue_ind].push_back(Entry {
                            pos: fwd_pos,
                            dir: entry.dir
                        });
                        self.backpointers[fwd_ind] = Some((-entry.dir.0, -entry.dir.1, -entry.dir.2));
                    }
                }
                let orthogonal_dirs = orthogonal_dirs(entry.dir);
                for n_dir in orthogonal_dirs {
                    let n_pos = (entry.pos.0 + n_dir.0, entry.pos.1 + n_dir.1, entry.pos.2 + n_dir.2);
                    if vox_ind_ok(n_pos, self.dim) {
                        let n_ind = vox_indi(n_pos, self.dim);
                        if self.voxels[n_ind] == 0 && self.backpointers[n_ind].is_none() {
                            if self.queues.len() == self.queue_ind + 1 {
                                self.queues.push(VecDeque::new());
                            }
                            self.queues[self.queue_ind + 1].push_back(Entry {
                                pos: n_pos,
                                dir: n_dir
                            });
                            self.backpointers[n_ind] = Some((-n_dir.0, -n_dir.1, -n_dir.2));
                        }
                    }
                }
            } else {
                self.queue_ind += 1;
            }
        }

        // Ok now its done -- backpointers have been set
        // Now starting at the ends, follow backpointers and set voxels to 2 until theres one with backptr of 0,0 - thats the start
        // Then generate the mesh
        // then export to STL
        for endpoint in self.endpoints.iter() {
            let mut pos = (endpoint.x as isize, endpoint.y as isize, endpoint.z as isize);
            loop {
                // infinite loop
                let idx = vox_indi(pos, self.dim);
                assert!(self.voxels[idx] != 1);
                self.voxels[idx] = 2;
                if let Some(dir) = self.backpointers[idx] {
                    if dir.0 == 0 && dir.1 == 0 && dir.2 == 0 {
                        // found the start
                        break;
                    }
                    pos = (pos.0 + dir.0, pos.1 + dir.1, pos.2 + dir.2);
                } else {
                    panic!("invalid trail");
                }
            }
        }
    }
}

pub fn orthogonal_dirs(v: (isize, isize, isize)) -> [(isize, isize, isize); 4] {
    [
        (v.1, v.2, v.0),
        (-v.1, -v.2, -v.0),
        (v.2, v.0, v.1),
        (-v.2, -v.0, -v.1),
    ]
}