pub use std::f32::consts::PI;

#[derive(Clone, Copy)]
pub struct V3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}
pub fn v3(x: f32, y: f32, z: f32) -> V3 { V3 { x, y, z } }

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

impl std::ops::Add<V3> for V3 {
    type Output = V3;

    fn add(self, _rhs: V3) -> V3 {
        V3 { x: self.x + _rhs.x, y: self.y + _rhs.y, z: self.z + _rhs.z }
    }
}

impl std::ops::Div<f32> for V3 {
    type Output = V3;

    fn div(self, _rhs: f32) -> V3 {
        V3 { x: self.x / _rhs, y: self.y / _rhs, z: self.z / _rhs }
    }
}