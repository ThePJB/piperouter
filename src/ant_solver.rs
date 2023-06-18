use ordered_float::Float;

use crate::endpoint::*;
use crate::math::*;

const PHEROMONE_DEPOSIT_AMOUNT: f32 = 1.0;
const PHEROMONE_DECAY_COEFF: f32 = 0.9999;
const PHEROMONE_STRAIGHT_BONUS: f32 = 3.0; // too strong and turns will get rekt
const PHEROMONE_FLAT_BONUS: f32 = 0.1; // for exploration

// so weve got to hope that the ants dont just degenerate from only single type of pheromone because I also believe in that
// to achieve trunking
// but could have multiple pheromones to rule it out too

pub struct AntSolver {
    pub dim: (usize, usize, usize),

    pub endpoints: Vec<VoxelEndpoint>,

    pub voxels: Vec<usize>,
    pub last_pheromone: Vec<f32>,
    pub last_visitation: Vec<usize>,
    pub agent_x: Vec<usize>,
    pub agent_y: Vec<usize>,
    pub agent_z: Vec<usize>,
    
    pub agent_dir_x: Vec<isize>,
    pub agent_dir_y: Vec<isize>,
    pub agent_dir_z: Vec<isize>,

    // agent_pheromone: Vec<f32>,

    pub agent_src_x: Vec<usize>,
    pub agent_src_y: Vec<usize>,
    pub agent_src_z: Vec<usize>,    

    pub seed: usize,
    pub step: usize,
}

impl AntSolver {
    pub fn new(dim: (usize, usize, usize), endpoints: Vec<VoxelEndpoint>, voxels: Vec<usize>) -> Self {
        let len = voxels.len();
        AntSolver {
            dim,
            endpoints,
            voxels,
            last_pheromone: vec![0.0; len],
            last_visitation: vec![0; len],
            agent_x: vec![],
            agent_y: vec![],
            agent_z: vec![],
            agent_dir_x: vec![],
            agent_dir_y: vec![],
            agent_dir_z: vec![],
            agent_src_x: vec![],
            agent_src_y: vec![],
            agent_src_z: vec![],
            seed: 123123,
            step: 0,
        }
    }
    pub fn spawn_ant(&mut self, src_idx: usize) {
        self.agent_x.push(self.endpoints[src_idx].x);
        self.agent_y.push(self.endpoints[src_idx].y);
        self.agent_z.push(self.endpoints[src_idx].z);
        self.agent_src_x.push(self.endpoints[src_idx].x);
        self.agent_src_y.push(self.endpoints[src_idx].y);
        self.agent_src_z.push(self.endpoints[src_idx].z);
        self.agent_dir_x.push(self.endpoints[src_idx].nx);
        self.agent_dir_y.push(self.endpoints[src_idx].ny);
        self.agent_dir_z.push(self.endpoints[src_idx].nz);
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
    pub fn calc_pheromone(&self, idx: usize) -> f32 {
        let nsteps = self.step - self.last_visitation[idx];
        self.last_pheromone[idx] * PHEROMONE_DECAY_COEFF.powi(nsteps as i32)
    }
    pub fn step(&mut self) {


        // decay pheromone
        // too slow
        // for i in 0..self.pheromone.len() {
        //     self.pheromone[i] *= PHEROMONE_DECAY_COEFF;
        // }

        // dont want to reallocate for each agent
        let mut candidate_inds: Vec<usize> = Vec::new();
        let mut candidate_x: Vec<isize> = Vec::new();
        let mut candidate_y: Vec<isize> = Vec::new();
        let mut candidate_z: Vec<isize> = Vec::new();
        let mut candidate_preferred: Vec<bool> = Vec::new();
        let mut candidate_pheromone: Vec<f32> = Vec::new();
        let mut candidate_weight: Vec<f32> = Vec::new();

        for idx in 0..self.agent_x.len() {
            candidate_inds.clear();
            candidate_x.clear();
            candidate_y.clear();
            candidate_z.clear();
            candidate_preferred.clear();
            candidate_pheromone.clear();
            candidate_weight.clear();

            let i = self.agent_x[idx];
            let j = self.agent_y[idx];
            let k = self.agent_z[idx];

            let agent_idx = self.voxel_ind(i,j,k);

            // deposit pheromone
            let pheromone = self.calc_pheromone(agent_idx) + PHEROMONE_DEPOSIT_AMOUNT;
            self.last_visitation[agent_idx] = self.step;
            self.last_pheromone[agent_idx] = pheromone;

            // if at endpoint we turn around
            for endpt_idx in 0..self.endpoints.len() {
                if self.endpoints[endpt_idx].x == i && self.endpoints[endpt_idx].y == j && self.endpoints[endpt_idx].z == k {
                    // turn around vs. assume normal
                    self.agent_dir_x[idx] *= -1;
                    self.agent_dir_y[idx] *= -1;
                    self.agent_dir_z[idx] *= -1;

                    if i != self.agent_src_x[idx] || j != self.agent_src_y[idx] || k != self.agent_src_z[idx] {
                        println!("ding {}", endpt_idx);
                        self.agent_src_x[idx] = i;
                        self.agent_src_y[idx] = j;
                        self.agent_src_z[idx] = k;
                    }

                    // src and dest must swap too
                    // oh it can repeat scoring going around in a circle here though.
                    // so thats why the encouraging dinging

                    break;
                }
            }

            // move forward by selecting from the 5 available directions based on amount of pheromone and obstruction
            // candidate coords
            // pos + dir
            // pos + everything except -dir

            // wait but me never check for walls!! but shuld still work on djikstra first
            // fn can_go take pos and dir

            // so we want basically these vecs: target index usize
            //                                  preferred    bool
            if let Some(front) = self.voxel_ind_check(i as isize + self.agent_dir_x[idx], j as isize + self.agent_dir_y[idx], k as isize + self.agent_dir_z[idx]) {
                candidate_inds.push(front);
                candidate_x.push(i as isize + self.agent_dir_x[idx]);
                candidate_y.push(j as isize + self.agent_dir_y[idx]);
                candidate_z.push(k as isize + self.agent_dir_z[idx]);
                candidate_pheromone.push(self.calc_pheromone(front));
                candidate_preferred.push(true);
            }
            // rotations and -s of rotations
            if let Some(side) = self.voxel_ind_check(i as isize + self.agent_dir_y[idx], j as isize + self.agent_dir_z[idx], k as isize + self.agent_dir_x[idx]) {
                candidate_inds.push(side);
                candidate_x.push(i as isize + self.agent_dir_y[idx]);
                candidate_y.push(j as isize + self.agent_dir_z[idx]);
                candidate_z.push(k as isize + self.agent_dir_x[idx]);
                candidate_pheromone.push(self.calc_pheromone(side));
                candidate_preferred.push(false);
            }
            // rotations and -s of rotations
            if let Some(side) = self.voxel_ind_check(i as isize - self.agent_dir_y[idx], j as isize - self.agent_dir_z[idx], k as isize - self.agent_dir_x[idx]) {
                candidate_inds.push(side);
                candidate_x.push(i as isize - self.agent_dir_y[idx]);
                candidate_y.push(j as isize - self.agent_dir_z[idx]);
                candidate_z.push(k as isize - self.agent_dir_x[idx]);
                candidate_pheromone.push(self.calc_pheromone(side));
                candidate_preferred.push(false);
            }
            // rotations and -s of rotations
            if let Some(side) = self.voxel_ind_check(i as isize + self.agent_dir_z[idx], j as isize + self.agent_dir_x[idx], k as isize + self.agent_dir_y[idx]) {
                candidate_inds.push(side);
                candidate_x.push(i as isize + self.agent_dir_z[idx]);
                candidate_y.push(j as isize + self.agent_dir_x[idx]);
                candidate_z.push(k as isize + self.agent_dir_y[idx]);
                candidate_pheromone.push(self.calc_pheromone(side));
                candidate_preferred.push(false);
            }
            // rotations and -s of rotations
            if let Some(side) = self.voxel_ind_check(i as isize - self.agent_dir_z[idx], j as isize - self.agent_dir_x[idx], k as isize - self.agent_dir_y[idx]) {
                candidate_inds.push(side);
                candidate_x.push(i as isize - self.agent_dir_z[idx]);
                candidate_y.push(j as isize - self.agent_dir_x[idx]);
                candidate_z.push(k as isize - self.agent_dir_y[idx]);
                candidate_pheromone.push(self.calc_pheromone(side));
                candidate_preferred.push(false);
            }

            let pheromone_total = candidate_pheromone.iter().fold(0.0, |acc, v| acc + v) + 
                candidate_preferred.iter().filter(|p| **p).map(|_| PHEROMONE_STRAIGHT_BONUS).fold(0.0, |acc, v| acc + v) +
                candidate_inds.len() as f32 * PHEROMONE_FLAT_BONUS;
            for idx in 0..candidate_inds.len() {
                let w = (candidate_pheromone[idx] + PHEROMONE_STRAIGHT_BONUS * if candidate_preferred[idx] {1.0} else {0.0} + PHEROMONE_FLAT_BONUS) / pheromone_total;
                candidate_weight.push(w);
            }
            // weights should add to 1, now we roll a random number between 0 and 1;
            let roll = krand(self.seed);
            let selection: usize = {
                let mut cum = 0.0;
                let mut selection = 0;
                for idx in 0..candidate_weight.len() {
                    cum += candidate_weight[idx];
                    if cum > roll {
                        selection = idx;
                        break;
                    }
                }
                selection
            };
            self.seed = khash(self.seed);


            
            // 'preferred = what in terms of weights? if the pheromone is drastically better it should prefer the corner highly'
            // all things being 0.0 it should favour preferred
    
            // optimization: skip bounds checking if we can ensure encasement in walls
    
            // maybe need dst-based pheromone table but then would it trunk. i think yes

            // now we need to consider over inds and preferred
            // and new dst will be x - prev x
            let prev_x = i as isize;
            let prev_y = j as isize;
            let prev_z = k as isize;

            let new_x = candidate_x[selection];
            let new_y = candidate_y[selection];
            let new_z = candidate_z[selection];
            
            let new_dir_x = new_x - prev_x;
            let new_dir_y = new_y - prev_y;
            let new_dir_z = new_z - prev_z;

            self.agent_x[idx] = new_x as usize;
            self.agent_y[idx] = new_y as usize;
            self.agent_z[idx] = new_z as usize;

            let old_dir_x = self.agent_dir_x[idx];
            let old_dir_y = self.agent_dir_y[idx];
            let old_dir_z = self.agent_dir_z[idx];

            self.agent_dir_x[idx] = new_dir_x;
            self.agent_dir_y[idx] = new_dir_y;
            self.agent_dir_z[idx] = new_dir_z;

            // println!("step: {} -- ant pos({},{},{}) facing({},{},{}) moved to ({},{},{})", self.step, prev_x, prev_y, prev_z, old_dir_x, old_dir_y, old_dir_z, new_x, new_y, new_z);

        }
        self.step += 1;
    }
}