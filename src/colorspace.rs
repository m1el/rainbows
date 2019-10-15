use image::{Pixel, Rgb};

#[derive(Debug)]
pub struct Xyz<T>(pub [T; 3]);

impl Xyz<f64> {
    pub fn to_srgb(&self) -> Rgb<f64> {
        // Adobe RGB
        // const M: [[f64; 3]; 3] = [
        //     [2.0413690, -0.5649464, -0.3446944],
        //     [-0.9692660, 1.8760108, 0.0415560],
        //     [0.0134474, -0.1183897, 1.0154096],
        // ];
        // CIE RGB
        // const M: [[f64; 3]; 3] = [
        //     [2.3706743, -0.9000405, -0.4706338],
        //     [-0.5138850, 1.4253036, 0.0885814],
        //     [0.0052982, -0.0146949, 1.0093968],
        // ];
        // sRGB
        const M: [[f64; 3]; 3] = [
            [3.240_454_2, -1.537_138_5, -0.498_531_4],
            [-0.969_266_0, 1.876_010_8, 0.041_556_0],
            [0.055_643_4, -0.204_025_9, 1.057_225_2],
        ];
        // linear light srgb?..
        //const M: [[f64; 3]; 3] = [
        //    [3.1338561, -1.6168667, -0.4906146],
        //    [-0.9787684, 1.9161415, 0.0334540],
        //    [0.0719453, -0.2289914, 1.4052427],
        //];
        let xyz = &self.0;
        Rgb([
            M[0][0] * xyz[0] + M[0][1] * xyz[1] + M[0][2] * xyz[2],
            M[1][0] * xyz[0] + M[1][1] * xyz[1] + M[1][2] * xyz[2],
            M[2][0] * xyz[0] + M[2][1] * xyz[1] + M[2][2] * xyz[2],
        ])
    }

    pub fn scale(&self, factor: f64) -> Self {
        Xyz([
            self.0[0] * factor,
            self.0[1] * factor,
            self.0[2] * factor,
        ])
    }
}

impl std::ops::Add for Xyz<f64> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Xyz([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
        ])
    }
}

/// Convert a linear floating point RGB color to u8 sRGB color
pub fn linear_to_u8(rgb: Rgb<f64>, gamma: f64) -> Rgb<u8> {
    let mapped = rgb.map(|x| x.max(0.0).min(1.0).powf(1.0 / gamma) * 255.5);
    Rgb([mapped.0[0] as u8, mapped.0[1] as u8, mapped.0[2] as u8])
}
