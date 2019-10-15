use std::fs::File;
//use std::io::{BufWriter};
use std::error::Error;
//use std::io::{Write as IoWrite};
//use std::fmt::{Write as FmtWrite};

use image::{Pixel, Rgb};
use rand::Rng;
use rayon::prelude::*;
use serde::Deserialize;

use crate::colorspace::{Xyz, linear_to_u8};
use crate::math_helpers::{*};


type LazyResult<T> = Result<T, Box<dyn Error>>;

/// Calculate Fresnel refractance for s- and p-polarized photons
fn fresnel_r(n: f64, cos_i: f64, cos_t: f64) -> (f64, f64) {
    let r_s = (((cos_i / n) - cos_t) / ((cos_i / n) + cos_t)).powi(2);
    let r_p = (((cos_t / n) - cos_i) / ((cos_t / n) + cos_i)).powi(2);
    (r_s, r_p)
}

/// Trace a ray going through the sphrere,
/// where `ray` is incoming ray and output is outgoing ray and probability
fn trace_sphere<R: Rng>(rng: &mut R, sphere: &Sphere, n: f64, ray: &Ray) -> (Ray, f64) {
    let mut ray = ray.clone();
    let mut inside = false;
    let mut p_s = 0.5_f64;
    let mut p_p = 0.5_f64;
    let mut bounces = 0;
    //let mut path = String::new();
    //write!(path, "M{:.3} {:.3}", ray.origin.x * 100.0, ray.origin.y * 100.0);

    while let Some(intersection) = sphere.intersect(&ray, inside) {
        //write!(path, "L{:.3} {:.3}", intersection.x * 100.0, intersection.y * 100.0);
        let n_e = if inside { 1.0 / n } else { n };
        let normal = sphere.normal_at(intersection);
        let cos_i = ray.direction.dot(normal);
        let sin_i = ray.direction.cross(normal).len();
        let sin_t = sin_i / n_e;
        let cos_t;
        let r_s;
        let r_p;
        let will_reflect = if sin_t >= 1.0 {
            cos_t = 0.0;
            r_s = 1.0;
            r_p = 1.0;
            // Total Inner Reflection
            true
        } else {
            cos_t = (1.0 - sin_t.powi(2)).sqrt();
            let r = fresnel_r(n_e, cos_i.abs(), cos_t.abs());
            r_s = r.0;
            r_p = r.1;
            let reflectance = (p_s * r_s + p_p * r_p) / p_p.hypot(p_s);
            rng.gen::<f64>() < reflectance
        };

        //println!("reflectance: {}, will_reflect: {}", reflectance, will_reflect);
        if will_reflect {
            let direction = ray.direction - normal.scale(cos_i * 2.0);
            p_s *= r_s;
            p_p *= r_p;
            ray = Ray {
                origin: intersection,
                direction,
            };
        } else {
            let tangent = (ray.direction - normal.scale(cos_i)).scale(1.0 / n_e);
            let normal = normal.scale(cos_i.signum()).scale(cos_t);
            ray = Ray {
                origin: intersection,
                direction: tangent + normal,
            };
            p_s *= 1.0 - r_s;
            p_p *= 1.0 - r_p;
            inside = !inside;
        }
        bounces += 1;
        if bounces > 100 {
            p_s = 0.0;
            p_p = 0.0;
            break;
        }
    }
    //let end = ray.origin + ray.direction;
    //write!(path, "L{:.3} {:.3}", end.x * 100.0, end.y * 100.0);

    //println!("{}", path);
    (ray, (p_s + p_p).sqrt())
}

const C: f64 = 299_792_458.0;
/// Calculate the intensity density of a black body
/// with given temperature on a given wavelength
fn bb_spectrum(nu: f64, temp: f64) -> f64 {
    const H: f64 = 6.626_070_15e-34;
    const K: f64 = 1.380_648_52e-23;
    (2.0 * H * nu.powi(3) * C.powi(-2))
        / (((H * nu) / (K * temp)).exp() - 1.0)
}

/// Read color matching function data from CSV
fn read_cmf_xyz() -> LazyResult<RegularInterpFn<Xyz<f64>>> {
    /// Record format used in CMF data
    #[derive(Deserialize)]
    struct CMFRecord {
        wavelength: f64,
        x: f64,
        y: f64,
        z: f64,
    }

    // Read CMF data
    let spectrum_file = File::open("data/lin2012xyz2e_1_7sf.csv")?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_reader(spectrum_file);
    let mut first = None;
    let mut last = None;
    let mut records = Vec::new();
    for result in rdr.deserialize() {
        let record: CMFRecord = result?;
        first = first.or(Some(record.wavelength));
        last.replace(record.wavelength);
        records.push(Xyz([record.x, record.y, record.z]));
    }

    let lam_first = first.unwrap();
    let lam_last = last.unwrap();

    // Create an interpolated funtion based on CMF data
    Ok(RegularInterpFn {
        start: lam_first,
        end: lam_last,
        data: records,
    })
}

/// Read water refraction index function from a csv file
fn read_water_rindex() -> LazyResult<SegmentInterpFn<f64>> {
    #[derive(Deserialize)]
    struct RefractRecord {
        wl: f64,
        n: f64,
    }
    // Read refraction data
    //let refraction_file = File::open("data/Daimon-19.0C.csv")?;
    let refraction_file = File::open("data/Martonchik-liquid-90K.csv")?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(refraction_file);
    let mut args = Vec::new();
    let mut data = Vec::new();
    for result in rdr.deserialize() {
        let record: RefractRecord = result?;
        args.push(record.wl * 1000.0);
        data.push(record.n);
    }

    // Create an interpolated funtion based on CMF data
    Ok(SegmentInterpFn { args, data })
}

/// Calculate dispersion functions for each visible (start..=end)
/// Calculation is done by raytracing a sphere,
/// then collecting output ray distribution over different angles.
/// The output is a 2d function represented as a Vec of 1D functions.
/// Vec index is wavelength, and function argument is an angle.
pub fn calc_dispersion(start: usize, end: usize)
    -> LazyResult<Vec<RegularInterpFn<f64>>>
{
    use std::f64::consts::PI;
    let rindex = read_water_rindex()?;
    const N_BINS: usize = 8096;
    const N_RAYS: usize = 10_000_000;
    const MULT: f64 = 3e6 / (N_RAYS as f64);

    Ok((start..=end)
        .into_par_iter()
        .map(|lambda| {
            let rin = rindex.lookup_lerp(lambda as f64);
            let sphere = Sphere {
                center: Vec3 {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                },
                radius: 1.0,
            };
            let mut ray = Ray {
                origin: Vec3 {
                    x: -1.0,
                    y: 0.0,
                    z: 0.0,
                },
                direction: Vec3 {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
            };
            let mut distribution = vec![0.0; N_BINS];
            let mut rng = rand::thread_rng();
            for _ in 0..N_RAYS {
                //ray.origin.y = rng.gen();
                ray.origin.y = rng.gen::<f64>().sqrt();
                let (out, probability) = trace_sphere(&mut rng, &sphere, rin, &ray);
                let angle = out.direction.y.atan2(out.direction.x).abs();
                let slot = (1.0 - (angle / PI)) * (N_BINS as f64);
                let slot = (slot as usize).min(N_BINS - 1);
                distribution[slot] += probability;
            }

            for el in distribution.iter_mut() {
                *el *= MULT;
            }

            //let fname = format!("data/distr/distr_{}.txt", lambda);
            //let out = File::create(fname)
            //    .expect("couldn't open distribution file?");
            //let mut writer = BufWriter::new(out);

            //for (i, d) in distribution.iter().skip(N_BINS / 2).rev().enumerate() {
            //    let angle = 180.0 / (N_BINS as f64) * (i as f64);
            //    writeln!(writer, "{:.4},{:.4}", angle, d)
            //        .expect("couldn't write?");
            //}
            RegularInterpFn {
                start: 0.0,
                end: PI,
                data: distribution,
            }
        })
        .collect())
}

/*
fn spectrum_xyz<F>(spectrum: F, cmf: &[Xyz<f64>]) -> Xyz<f64>
    where F: Fn(f64) -> f64
{
    let mut xyz = Xyz([0.0, 0.0, 0.0]);
    for lambda in 390..=830 {
        let lambda_f = lambda as f64;
        let nu = C / (lambda_f * 1e-9);
        let intensity = spectrum(nu);
        xyz = xyz + cmf[lambda - 390].scale(intensity);
    }
    xyz
}
*/

pub fn attempt_03_refractance() -> LazyResult<()> {
    let cmf_xyz = read_cmf_xyz()?;
    let lam_first = 390;
    let lam_last = 830;

    let dispersion_data = calc_dispersion(lam_first, lam_last)?;

    let width = 1280;
    let height = 640;
    let scale = 1.0 / f64::from(height);
    let sky_color = Rgb([0.5294, 0.80784, 0.92157]).lerp(&Rgb([0., 0., 0.]), 1.0);
    let gamma = 2.2;
    let temp = 6500.0;
    let brightness = 500.0;

    let mut imgbuf = image::ImageBuffer::new(width, height);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        use std::f64::consts::PI;
        let mut fx = (f64::from(x) * scale) - 1.0;
        let mut fy = 1.0 - (f64::from(y) * scale);
        // zoom at an interesting location
        fx = fx / 10.0 + 0.48;
        fy = fy / 10.0;
        //fx /= 1.6;
        //fy /= 1.6;
        let r = fx.hypot(fy);
        let alpha = r * 0.5 * PI;
        // let phi = fy.atan2(fx) - PI;

        let mut color = sky_color;
        if alpha < PI * 0.5 {
            let mut xyz = Xyz([0.0, 0.0, 0.0]);
            for lambda in lam_first..=lam_last {
                let band = &dispersion_data[lambda - lam_first];
                let lambda_f = lambda as f64;
                let nu = C / (lambda_f * 1e-9);
                let intensity = bb_spectrum(nu, temp)
                    * band.lookup_lerp(alpha) * brightness;
                xyz = xyz + cmf_xyz.lookup_lerp(lambda_f).scale(intensity);
            }
            let rgb = xyz.to_srgb();
            color.apply2(&rgb, |x, y| x + y);
        }
        *pixel = linear_to_u8(color, gamma);
    }

    imgbuf.save("attempt_03_dispersion.png")?;
    Ok(())
}
