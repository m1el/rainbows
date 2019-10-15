/// A trait that defines linear interpolation between two values
pub trait Lerp {
    /// Interpolate linearly between `self` and `other` by `factor`.
    /// `lerp(a, b, t) = a * (1 - t) + b * t`
    /// or optimized, `lerp(a, b, t) = a + t * (b - a)`.
    fn lerp(&self, other: &Self, factor: f64) -> Self;
}

impl Lerp for f64 {
    fn lerp(&self, other: &Self, factor: f64) -> Self {
        *self + factor * (*other - *self)
    }
}

use image::Rgb;
impl Lerp for Rgb<f64> {
    fn lerp(&self, other: &Self, factor: f64) -> Self {
        Rgb([
            self.0[0] + factor * (other.0[0] - self.0[0]),
            self.0[1] + factor * (other.0[1] - self.0[1]),
            self.0[2] + factor * (other.0[2] - self.0[2]),
        ])
    }
}

use crate::colorspace::Xyz;
impl Lerp for Xyz<f64> {
    fn lerp(&self, other: &Self, factor: f64) -> Self {
        Xyz([
            self.0[0] + factor * (other.0[0] - self.0[0]),
            self.0[1] + factor * (other.0[1] - self.0[1]),
            self.0[2] + factor * (other.0[2] - self.0[2]),
        ])
    }
}

/// Struct used to represent and calculate an interpolated function.
/// The function is defined by values for regularly separated
/// arguments on a given range.
///
/// For example,
///     RegularInterpFn {
///         start: 0.0,
///         end: 2.0,
///         data: &[2.0, 1.0, 0.0],
///     }
/// defines a function by three datapoints,
///     f(0.0) = 2.0; f(1.0) = 1.0; f(2.0) = 0.0
pub struct RegularInterpFn<T>
where
    T: Lerp,
{
    pub start: f64,
    pub end: f64,
    pub data: Vec<T>,
}

impl<T> RegularInterpFn<T>
where
    T: Lerp,
{
    /// Given a function argument `x`, calculate index into the dataset
    /// and fraction of segment position between values.
    pub fn index_frac(&self, mut x: f64) -> (usize, f64) {
        let range = self.end - self.start;
        let steps = self.data.len() - 1;
        let step = range / (steps as f64);
        x = (x - self.start).max(0.0);
        let index = (x.div_euclid(step) as usize).min(steps);
        let frac = x.rem_euclid(step) / step;
        (index, frac)
    }

    /// Look up a value for argument `x`,
    /// linearly interpolating between known values
    pub fn lookup_lerp(&self, x: f64) -> T {
        let (idx, frac) = self.index_frac(x);
        let next = (idx + 1).min(self.data.len() - 1);
        self.data[idx].lerp(&self.data[next], frac)
    }
}

/// Struct used to represent and calculate an interpolated function.
/// The function is defined by pairs of argument and data.
///
/// For example,
///     RegularInterpFn {
///         args: vec![0.0, 1.0, 2.0]
///         data: vec![2.0, 1.0, 0.0],
///     }
/// defines a function by three datapoints,
///     f(0.0) = 2.0; f(1.0) = 1.0; f(2.0) = 0.0
pub struct SegmentInterpFn<T>
where
    T: Lerp,
{
    pub args: Vec<f64>,
    pub data: Vec<T>,
}

impl<T> SegmentInterpFn<T>
where
    T: Lerp,
{
    /// Given a function argument `x`, calculate index into the dataset
    /// and fraction of segment position between values.
    pub fn index_frac(&self, x: f64) -> (usize, f64) {
        let mut args = self.args.iter().cloned().enumerate();
        if let Some((index, right)) = args.find(|&(_, arg)| arg >= x) {
            if index == 0 {
                (0, 0.0)
            } else {
                let left = self.args[index - 1];
                let step = right - left;
                let offset = x - left;
                (index - 1, offset / step)
            }
        } else {
            (self.args.len() - 1, 0.0)
        }
    }

    /// Look up a value for argument `x`,
    /// linearly interpolating between known values
    pub fn lookup_lerp(&self, x: f64) -> T {
        let (idx, frac) = self.index_frac(x);
        self.data[idx].lerp(&self.data[idx + 1], frac)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Vec3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl Vec3<f64> {
    pub fn dot(&self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: Self) -> Self {
        Vec3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn len(&self) -> f64 {
        self.dot(*self).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let lensq = self.dot(*self);
        let scale = if lensq <= std::f64::EPSILON {
            1.0
        } else {
            1.0 / lensq.sqrt()
        };
        self.scale(scale)
    }

    pub fn scale(&self, factor: f64) -> Self {
        Vec3 {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
        }
    }
}

impl std::ops::Add for Vec3<f64> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl std::ops::Sub for Vec3<f64> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Ray {
    pub origin: Vec3<f64>,
    pub direction: Vec3<f64>,
}

#[derive(Debug, Clone)]
pub struct Sphere {
    pub center: Vec3<f64>,
    pub radius: f64,
}

impl Sphere {
    pub fn normal_at(&self, point: Vec3<f64>) -> Vec3<f64> {
        (point - self.center).normalize()
    }

    #[allow(dead_code)]
    pub fn intersect(&self, ray: &Ray, inside: bool) -> Option<Vec3<f64>> {
        // ((o - c) + t*d)^2 = r^2
        // t^2 * (d.d) + 2*t*(d.(o-c)) ((o-c).(o-c)) - r^2 = 0
        // a = d.d, b = 2*(d.(o-c)), c = ((o-c).(o-c)) - r^2
        // t = (-b+-sqrt(b^2 - 4ac)) / 2a
        let delta_origin = ray.origin - self.center;
        let a = ray.direction.dot(ray.direction);
        let b = 2.0 * ray.direction.dot(delta_origin);
        let c = delta_origin.dot(delta_origin) - self.radius.powi(2);
        let discriminant = b * b - 4.0 * a * c;
        if discriminant <= 0.0 {
            return None;
        }

        let sqrtd = discriminant.sqrt();

        let left = (-b - sqrtd) / (2.0 * a);
        if !inside && left > 0.0 {
            return Some(ray.origin + ray.direction.scale(left));
        }

        let right = (-b + sqrtd) / (2.0 * a);
        if inside && right > 0.0 {
            return Some(ray.origin + ray.direction.scale(right));
        }

        None
    }
}
