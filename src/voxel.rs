pub const U_TO_M: f32 = 3.6429696;
pub const VOXEL_S_M: f32 = 0.01;
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