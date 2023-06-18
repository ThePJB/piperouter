use crate::endpoint::*;
use std::collections::VecDeque;

// maybe it should actually search from longest one first? but might not be shortest
// oh rn im poppin so might be a long path
// if inner one was a vecdeq then

pub struct Entry {
    pos: (isize, isize, isize),
    dir: (isize, isize, isize),
}

pub struct FastSolver {
    dim: (usize, usize, usize),

    endpoints: Vec<VoxelEndpoint>,
    endpoint_inds: Vec<usize>,

    voxels: Vec<usize>,
    backpointers: Vec<Option<(isize, isize, isize)>>,

    queues: Vec<VecDeque<Entry>>,
    queue_ind: usize,
}

fn vox_ind(pos: (usize, usize, usize), dim: (usize, usize, usize)) -> usize {
    pos.0*dim.2*dim.1 + pos.1*dim.2 + pos.2
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

        self.queues.push(vec![start_entry].into());

        loop {
            if self.queue_ind >= self.queues.len() {
                dbg!("finishing -- exhausted all queues");
                return;
            }
            if let Some(entry) = self.queues[self.queue_ind].pop_front() {

            } else {
                self.queue_ind += 1;
            }
        }



        let k = KeyType {gen: 0, dist: 0};
        let pos = (self.endpoints[start_endpoint_idx].x, self.endpoints[start_endpoint_idx].y, self.endpoints[start_endpoint_idx].z);
        let dir = (self.endpoints[start_endpoint_idx].nx, self.endpoints[start_endpoint_idx].ny, self.endpoints[start_endpoint_idx].nz);
        let v = ValType {pos, dir};
        self.pq.push(k, v);
        let start_ind = self.voxel_ind(pos.0, pos.1, pos.2);
        self.backpointers[start_ind] = Some((0, 0, 0));
        
        // this should start at the goddamn thing so no backptr


        while let Some((curr_k, curr_v)) = self.pq.pop() {
            dbg!(self.pq.heap.len());
            dbg!(curr_k.gen);
            std::thread::sleep(std::time::Duration::from_millis(300));
            // visit neighbours, put backptr
            // visit fwd neighbour
            let curr_pos = curr_v.pos;
            let fwd_pos = (
                curr_pos.0 as isize + curr_v.dir.0 as isize,
                curr_pos.1 as isize + curr_v.dir.1 as isize,
                curr_pos.2 as isize + curr_v.dir.2 as isize,
            );
            if let Some(fwd_ind) = self.voxel_ind_check(fwd_pos.0, fwd_pos.1, fwd_pos.2) {
                // if not obstructed
                if self.voxels[fwd_ind] == 0 {
                    self.backpointers[fwd_ind] = Some((-curr_v.dir.0, -curr_v.dir.1, -curr_v.dir.2));
                    let k = KeyType {
                        gen: curr_k.gen,
                        dist: curr_k.dist + 1,
                    };
                    let v = ValType {
                        pos: (fwd_pos.0 as usize, fwd_pos.1 as usize, fwd_pos.2 as usize),
                        dir: curr_v.dir,
                    };
                    self.pq.push(k, v)
                }
            };

            // do same but for neighbours, calculate orthogonal directions
            let orthogonal_dirs = orthogonal_dirs(curr_v.dir);
            for ndir in orthogonal_dirs {
                let npos = (
                    curr_pos.0 as isize + ndir.0,
                    curr_pos.1 as isize + ndir.1,
                    curr_pos.2 as isize + ndir.2,
                );
                if let Some(nind) = self.voxel_ind_check(fwd_pos.0, fwd_pos.1, fwd_pos.2) {
                    // if not obstructed
                    if self.voxels[nind] == 0 {
                        self.backpointers[nind] = Some((-ndir.0, -ndir.1, -ndir.2));
                        let k = KeyType {
                            gen: curr_k.gen + 1,
                            dist: curr_k.dist + 1,
                        };
                        let v = ValType {
                            pos: (npos.0 as usize, npos.1 as usize, npos.2 as usize),
                            dir: curr_v.dir,
                        };
                        self.pq.push(k, v)
                    }
                };
            }
            

            // probably need a dst tracker too


            // put new values with same generation if no turn or worse generation with turn

            // check if theres a path to all and that means its done
            if self.endpoint_inds.iter().map(|x| self.backpointers[*x]).all(|x| x.is_some()) {
                dbg!("yes");
            }
        }
    }
    pub fn voxel_ind(&self, i: usize, j: usize, k: usize) -> usize {
        i*self.dim.2*self.dim.1 + j*self.dim.2 + k
    }
    pub fn voxel_ind_check(&self, i: isize, j: isize, k: isize) -> Option<usize> {
        if i < 0 || j < 0 || k < 0 { return None };
        let i = i as usize;
        let j = j as usize;
        let k = k as usize;
        if i < self.dim.0 && j < self.dim.1 && k < self.dim.2 {
            Some(i*self.dim.2*self.dim.1 + j*self.dim.2 + k)
        } else {
            None
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