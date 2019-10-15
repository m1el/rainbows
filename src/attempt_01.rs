use image::Rgb;
use std::error::Error;

use crate::colorspace::linear_to_u8;
use crate::math_helpers::RegularInterpFn;

/// A naive attempt to draw a rainbow by enumerating rainbow colors
pub fn attempt_01_naive() -> Result<(), Box<dyn Error>> {
    let width = 1280;
    let height = 640;
    let scale = 1.0 / f64::from(height);
    let rainbow_colors = RegularInterpFn {
        start: 0.0,
        end: 1.0,
        data: vec![
            Rgb([1.0, 0.0, 0.0]), // red
            Rgb([1.0, 0.5, 0.0]), // orange
            Rgb([0.8, 1.0, 0.0]), // yellow
            Rgb([0.0, 1.0, 0.0]), // green
            Rgb([0.0, 0.8, 0.8]), // sky blue
            Rgb([0.0, 0.0, 1.0]), // dark blue
            Rgb([0.6, 0.0, 0.8]), // violet
        ],
    };
    let sky_color = Rgb([0.5294, 0.80784, 0.92157]);
    let gamma = 2.0;

    let mut imgbuf = image::ImageBuffer::new(width, height);
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let fx = (f64::from(x) * scale) - 1.0;
        let fy = 1.0 - (f64::from(y) * scale);
        let r = fx.hypot(fy);
        // let phi = fy.atan2(fx) - PI;

        let color = if r > 0.8 && r < 1.0 {
            let pos = (1.0 - r) / (1.0 - 0.8);
            rainbow_colors.lookup_lerp(pos)
        } else {
            sky_color
        };
        *pixel = linear_to_u8(color, gamma);
    }

    imgbuf.save("attempt_01_naive.png")?;
    Ok(())
}
