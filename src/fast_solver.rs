use crate::voxel::*;
use std::collections::VecDeque;

pub struct Entry {
    pos: (isize, isize, isize),
    dir: (i8, i8, i8),
}

pub struct FastSolver {
    pub voxels: Voxels,

    pub endpoints: Vec<VoxelEndpoint>,
    pub endpoint_inds: Vec<usize>,

    // this could be reduced further obviously
    pub backpointers: Vec<Option<(i8, i8, i8)>>,

    pub queues: Vec<VecDeque<Entry>>,
    pub queue_ind: usize,
}

impl FastSolver {
    pub fn new(voxels: Voxels, endpoints: Vec<VoxelEndpoint>) -> Self {
        let len = voxels.voxels.len();
        let endpoint_inds = endpoints.iter().map(|endpoint| (endpoint.x, endpoint.y, endpoint.z)).map(|pos| voxels.get_idx_unchecked_u(pos)).collect();
        FastSolver {
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
            pos: (
                self.endpoints[start_endpoint_idx].x as isize, 
                self.endpoints[start_endpoint_idx].y as isize, 
                self.endpoints[start_endpoint_idx].z as isize
            ),
            dir: (
                self.endpoints[start_endpoint_idx].nx, 
                self.endpoints[start_endpoint_idx].ny, 
                self.endpoints[start_endpoint_idx].nz
            ),
        };

        let start_ind = self.voxels.get_idx_unchecked_i(start_entry.pos);
        self.queues.push(vec![start_entry].into());
        self.backpointers[start_ind] = Some((0, 0, 0)); // Sentinel for target voxel

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
                // let ind = self.voxels.get_idx_unchecked_i(entry.pos);
                // assert!(self.voxels.voxels[ind] != 1);
                // check and push front
                // not oob, wall, or prev
                let fwd_pos = (entry.pos.0 + entry.dir.0 as isize, entry.pos.1 + entry.dir.1 as isize, entry.pos.2 + entry.dir.2 as isize);
                if self.voxels.pos_in_bounds_i(fwd_pos) {
                    let fwd_idx = self.voxels.get_idx_unchecked_i(fwd_pos);
                    if self.voxels.voxels[fwd_idx] == 0 && self.backpointers[fwd_idx].is_none() {
                        self.queues[self.queue_ind].push_back(Entry {
                            pos: fwd_pos,
                            dir: entry.dir
                        });
                        self.backpointers[fwd_idx] = Some((-entry.dir.0, -entry.dir.1, -entry.dir.2));
                    }
                }
                let orthogonal_dirs = orthogonal_dirs(entry.dir);
                for n_dir in orthogonal_dirs {
                    let n_pos = (entry.pos.0 + n_dir.0 as isize, entry.pos.1 + n_dir.1 as isize, entry.pos.2 + n_dir.2 as isize);
                    if self.voxels.pos_in_bounds_i(n_pos) {
                        let n_idx = self.voxels.get_idx_unchecked_i(n_pos);
                        if self.voxels.voxels[n_idx] == 0 && self.backpointers[n_idx].is_none() {
                            if self.queues.len() == self.queue_ind + 1 {
                                self.queues.push(VecDeque::new());
                            }
                            self.queues[self.queue_ind + 1].push_back(Entry {
                                pos: n_pos,
                                dir: n_dir
                            });
                            self.backpointers[n_idx] = Some((-n_dir.0, -n_dir.1, -n_dir.2));
                        }
                    }
                }
            } else {
                self.queue_ind += 1;
            }
        }

        // Ok now its done -- backpointers have been set
        // Now starting at the ends, follow backpointers and set voxels to 2 until theres one with backptr of 0,0,0 - thats the start
        for endpoint in self.endpoints.iter() {
            let mut pos = (endpoint.x as isize, endpoint.y as isize, endpoint.z as isize);
            loop {
                let idx = self.voxels.get_idx_unchecked_i(pos);
                assert!(self.voxels.voxels[idx] != 1);
                self.voxels.voxels[idx] = 2;
                if let Some(dir) = self.backpointers[idx] {
                    if dir.0 == 0 && dir.1 == 0 && dir.2 == 0 {
                        // found the start
                        break;
                    }
                    pos = (pos.0 + dir.0 as isize, pos.1 + dir.1 as isize, pos.2 + dir.2 as isize);
                } else {
                    panic!("invalid trail");
                }
            }
        }
    }
}

// Given a direction, the 4 directions that are at 90 degrees to it
pub fn orthogonal_dirs(v: (i8, i8, i8)) -> [(i8, i8, i8); 4] {
    [
        (v.1, v.2, v.0),
        (-v.1, -v.2, -v.0),
        (v.2, v.0, v.1),
        (-v.2, -v.0, -v.1),
    ]
}