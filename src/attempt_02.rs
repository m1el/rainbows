use std::error::Error;
use std::fs::File;

use image::Rgb;
use serde::Deserialize;

use crate::colorspace::{linear_to_u8, Xyz};
use crate::math_helpers::{Lerp, RegularInterpFn};

type LazyResult<T> = Result<T, Box<dyn Error>>;

/// Calculate the intensity density of a black body
/// with given temperature on a given wavelength
fn bb_spectrum(lambda: f64, temp: f64) -> f64 {
    const H: f64 = 6.626_070_15e-34;
    const C: f64 = 299_792_458.0;
    const K: f64 = 1.380_648_52e-23;
    let nu = C / lambda;
    let a = 2.0 * H * nu.powi(3) * C.powi(-2);
    let b = ((H * nu) / (K * temp)).exp() - 1.0;
    a / b
}

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
    let spectrum_data = File::open("data/lin2012xyz2e_1_7sf.csv")?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_reader(spectrum_data);
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

/// An attempt to draw a rainbow by drawing spectrum
pub fn attempt_02_spectrum() -> LazyResult<()> {
    let cmf_xyz = read_cmf_xyz()?;
    let lam_first = cmf_xyz.start;
    let lam_last = cmf_xyz.end;

    let width = 1280;
    let height = 640;
    let scale = 1.0 / f64::from(height);
    let sky_color = Rgb([0.5294, 0.80784, 0.92157]).lerp(&Rgb([0., 0., 0.]), 0.3);
    let gamma = 2.2;
    let temp = 6500.0;
    let brightness = 4e6;

    let mut imgbuf = image::ImageBuffer::new(width, height);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let fx = (f64::from(x) * scale) - 1.0;
        let fy = 1.0 - (f64::from(y) * scale);
        let r = fx.hypot(fy);
        // let phi = fy.atan2(fx) - PI;

        let mut color = sky_color;
        if r > 0.8 && r < 1.0 {
            let pos = (r - 0.8) / (1.0 - 0.8);
            // estimate lambda depending on polar coordinates
            let lambda = lam_first.lerp(&lam_last, pos);
            // calculate the intensity based on Sun's spectrum
            let intensity = bb_spectrum(lambda / 1e9, temp) * brightness;
            let xyz = cmf_xyz.lookup_lerp(lambda);
            let rgb = xyz.to_srgb();
            color.0[0] += rgb.0[0] * intensity;
            color.0[1] += rgb.0[1] * intensity;
            color.0[2] += rgb.0[2] * intensity;
        }
        *pixel = linear_to_u8(color, gamma);
    }

    imgbuf.save("attempt_02_spectrum.png")?;
    Ok(())
}
