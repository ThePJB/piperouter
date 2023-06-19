pub use std::f32::consts::PI;

pub fn khash(mut state: usize) -> usize {
    state = (state ^ 2747636419).wrapping_mul(2654435769);
    state = (state ^ (state >> 16)).wrapping_mul(2654435769);
    state = (state ^ (state >> 16)).wrapping_mul(2654435769);
    state
}
pub fn khashu32(mut state: u32) -> u32 {
    state = (state ^ 2747636419).wrapping_mul(2654435769);
    state = (state ^ (state >> 16)).wrapping_mul(2654435769);
    state = (state ^ (state >> 16)).wrapping_mul(2654435769);
    state
}
pub fn krand(seed: usize) -> f32 {
    (khash(seed)&0x00000000FFFFFFFF) as f32 / 4294967295.0
}

pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a * (1.0 - t) + b * t
}

#[derive(Clone, Copy, Debug)]
pub struct V2 {
    pub x: f32,
    pub y: f32,
}
pub fn v2(x: f32, y: f32) -> V2 { V2 { x, y } }
pub fn v2_rt(r: f32, theta: f32) -> V2 { V2{x: r* (theta.sin()), y: r*(theta.cos()) } }
#[derive(Clone, Copy)]
pub struct V3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}
pub fn v3(x: f32, y: f32, z: f32) -> V3 { V3 { x, y, z } }
#[derive(Clone, Copy, Debug)]
pub struct V4 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}
pub fn v4(x: f32, y: f32, z: f32, w: f32) -> V4 { V4 { x, y, z, w } }

impl V2 {
    pub fn dot(&self, other: V2) -> f32 {
        self.x*other.x + self.y * other.y
    }
    pub fn rect_centered(&self, w: f32, h: f32) -> V4 {
        v4(self.x-w/2.0, self.y-h/2.0, w, h)
    }
    pub fn lerp(&self, other: V2, t: f32) -> V2 {
        v2(lerp(self.x, other.x, t), lerp(self.y, other.y, t))
    }
    pub fn homogeneous_transform(&self, mat: &[f32; 9]) -> V2 {
        v2(
            mat[0] * self.x + mat[1] * self.y + mat[2],
            mat[3] * self.x + mat[4] * self.y + mat[5],
        )
    }
    pub fn norm(&self) -> f32 {
        let mag = (self.x*self.x + self.y*self.y).sqrt();
        return mag;
    }
    pub fn normalize(&self) -> V2 {
        let mag = (self.x*self.x + self.y*self.y).sqrt();
        if mag == 0.0 {
            return v2(0.0, 0.0);
        }

        *self / mag
    }
}
impl V3 {
    pub fn dot(&self, other: V3) -> f32 {
        self.x*other.x + self.y * other.y + self.z*other.z
    }
    pub fn abs(&self) -> V3 {
        v3(self.x.abs(), self.y.abs(), self.z.abs())
    }
    pub fn norm(&self) -> V3 {
        *self / self.dot(*self).sqrt()
    }
    pub fn cross(self, other: V3) -> V3 {
        V3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
    pub fn min(self, other: V3) -> V3 {
        V3 {
            x: self.x.min(other.x),
            y: self.y.min(other.y),
            z: self.z.min(other.z),
        }
    }
    pub fn max(self, other: V3) -> V3 {
        V3 {
            x: self.x.max(other.x),
            y: self.y.max(other.y),
            z: self.z.max(other.z),
        }
    }
}
impl std::fmt::Debug for V3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "V3 [{} {} {}]", self.x, self.y, self.z)
    }
}
impl V4 {
    pub fn dot(&self, other: V4) -> f32 {
        self.x*other.x + self.y * other.y + self.z*other.z + self.w*other.w
    }
    pub fn tl(&self) -> V2 {v2(self.x, self.y)}
    pub fn br(&self) -> V2 {v2(self.x + self.z, self.y + self.w)}
    pub fn tr(&self) -> V2 {v2(self.x + self.z, self.y)}
    pub fn bl(&self) -> V2 {v2(self.x, self.y + self.w)}
    pub fn grid_child(&self, i: usize, j: usize, w: usize, h: usize) -> V4 {
        let cw = self.z / w as f32;
        let ch = self.w / h as f32;
        v4(self.x + cw * i as f32, self.y + ch * j as f32, cw, ch)
    }
    pub fn hsv_to_rgb(&self) -> V4 {
        let v = self.z;
        let hh = (self.x % 360.0) / 60.0;
        let i = hh.floor() as i32;
        let ff = hh - i as f32;
        let p = self.z * (1.0 - self.y);
        let q = self.z * (1.0 - self.y * ff);
        let t = self.z * (1.0 - self.y * (1.0 - ff));
        match i {
            0 => v4(v, t, p, self.w),
            1 => v4(q, v, p, self.w),
            2 => v4(p, v, t, self.w),
            3 => v4(p, q, v, self.w),
            4 => v4(t, p, v, self.w),
            5 => v4(v, p, q, self.w),
            _ => panic!("unreachable"),
        }
    }
    fn contains(&self, p: V2) -> bool {
        !(p.x < self.x || p.x > self.x + self.z || p.y < self.y || p.y > self.y + self.w)
    }
    fn point_within(&self, p: V2) -> V2 {
        v2(p.x*self.z+self.x, p.y*self.w+self.y)
    }
    fn point_without(&self, p: V2) -> V2 {
        v2((p.x - self.x) / self.z, (p.y - self.y) / self.w)
    }
    fn fit_aspect(&self, a: f32) -> V4 {
        let a_self = self.z/self.w;

        if a_self > a {
            // parent wider
            v4((self.z - self.z*(1.0/a))/2.0, 0.0, self.z*1.0/a, self.w)
        } else {
            // child wider
            v4(0.0, (self.w - self.w*(1.0/a))/2.0, self.z, self.w*a)
        }
    }
}

pub fn ray_triangle_intersection(ray_origin: V3, ray_dir: V3, v0: V3, v1: V3, v2: V3) -> Option<f32> {
    const EPSILON: f32 = 0.0001;

    let edge1 = v1 - v0;
    let edge2 = v2 - v0;

    let h = ray_dir.cross(edge2);
    let a = edge1.dot(h);

    if a.abs() < EPSILON {
        return None; // Ray is parallel to the triangle
    }

    let f = 1.0 / a;
    let s = ray_origin - v0;
    let u = f * s.dot(h);

    if u < 0.0 || u > 1.0 {
        return None; // Intersection is outside the triangle
    }

    let q = s.cross(edge1);
    let v = f * ray_dir.dot(q);

    if v < 0.0 || u + v > 1.0 {
        return None; // Intersection is outside the triangle
    }

    let t = f * edge2.dot(q);
    if t > EPSILON {
        Some(t)
    } else {
        None // Intersection is behind the ray
    }
}

impl std::ops::Mul<V3> for V3 {
    type Output = V3;

    fn mul(self, _rhs: V3) -> V3 {
        V3 { x: self.x * _rhs.x, y: self.y * _rhs.y, z: self.z * _rhs.z }
    }
}
impl std::ops::Sub<V3> for V3 {
    type Output = V3;

    fn sub(self, _rhs: V3) -> V3 {
        V3 { x: self.x - _rhs.x, y: self.y - _rhs.y, z: self.z - _rhs.z }
    }
}

impl std::ops::Mul<f32> for V3 {
    type Output = V3;

    fn mul(self, _rhs: f32) -> V3 {
        V3 { x: self.x * _rhs, y: self.y * _rhs, z: self.z * _rhs }
    }
}

impl std::ops::Mul<V3> for f32 {
    type Output = V3;

    fn mul(self, _rhs: V3) -> V3 {
        V3 { x: self * _rhs.x, y: self * _rhs.y, z: self * _rhs.z }
    }
}

impl PartialEq<V3> for V3 {
    fn eq(&self, _rhs: &V3) -> bool {
        self.x == _rhs.x && self.y == _rhs.y && self.z == _rhs.z
    }
}

impl std::ops::Sub<V2> for V2 {
    type Output = V2;

    fn sub(self, _rhs: V2) -> V2 {
        V2 { x: self.x - _rhs.x, y: self.y - _rhs.y }
    }
}

impl std::ops::Add<V2> for V2 {
    type Output = V2;

    fn add(self, _rhs: V2) -> V2 {
        V2 { x: self.x + _rhs.x, y: self.y + _rhs.y }
    }
}

impl std::ops::Add<V3> for V3 {
    type Output = V3;

    fn add(self, _rhs: V3) -> V3 {
        V3 { x: self.x + _rhs.x, y: self.y + _rhs.y, z: self.z + _rhs.z }
    }
}

impl std::ops::Mul<f32> for V2 {
    type Output = V2;

    fn mul(self, _rhs: f32) -> V2 {
        V2 { x: self.x * _rhs, y: self.y * _rhs }
    }
}

impl std::ops::Mul<V2> for f32 {
    type Output = V2;

    fn mul(self, _rhs: V2) -> V2 {
        V2 { x: self * _rhs.x, y: self * _rhs.y }
    }
}

impl std::ops::Div<f32> for V2 {
    type Output = V2;

    fn div(self, _rhs: f32) -> V2 {
        V2 { x: self.x / _rhs, y: self.y / _rhs }
    }
}

impl std::ops::Div<f32> for V3 {
    type Output = V3;

    fn div(self, _rhs: f32) -> V3 {
        V3 { x: self.x / _rhs, y: self.y / _rhs, z: self.z / _rhs }
    }
}

impl std::ops::Neg for V2 {
    type Output = V2;

    fn neg(self) -> V2 {
        V2 { x: -1.0 * self.x, y: -1.0 * self.y}
    }
}

pub fn mat_rot_z(theta: f32) -> [f32; 9] {
    let (s, c) = theta.sin_cos();
    [
        c, -s, 0.0, 
        s, c, 0.0,
        0.0, 0.0, 1.0
    ]
}
pub fn mat_trans(x: f32, y: f32) -> [f32; 9] {
    [
        1., 0., x,
        0., 1., y,
        0., 0., 1.,
    ]
}
pub fn mat_scale(s: f32) -> [f32; 9] {
    [
        s, 0., 0.,
        0., s, 0.,
        0., 0., s,
    ]
}

pub fn mat_mul33(a: &[f32; 9], b: &[f32; 9]) -> [f32; 9] {
    [
        a[0] * b[0] + a[1] * b[3] + a[2] * b[6],
        a[0] * b[1] + a[1] * b[4] + a[2] * b[7],
        a[0] * b[2] + a[1] * b[5] + a[2] * b[8],
        a[3] * b[0] + a[4] * b[3] + a[5] * b[6],
        a[3] * b[1] + a[4] * b[4] + a[5] * b[7],
        a[3] * b[2] + a[4] * b[5] + a[5] * b[8],
        a[6] * b[0] + a[7] * b[3] + a[8] * b[6],
        a[6] * b[1] + a[7] * b[4] + a[8] * b[7],
        a[6] * b[2] + a[7] * b[5] + a[8] * b[8],
    ]
}

pub fn mat_srt(x: f32, y: f32, scale: f32, theta: f32) -> [f32; 9] {
    let c = theta.cos();
    let s = theta.sin();
    let a = scale * c;
    let b = scale * s;
    [
        a, -b, x,
        b, a, y,
        0.0, 0.0, 1.0,
    ]
}
