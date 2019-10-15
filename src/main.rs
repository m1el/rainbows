extern crate csv;
extern crate image;
extern crate rand;
extern crate rayon;
extern crate serde;

mod attempt_01;
mod attempt_02;
mod attempt_03;
mod colorspace;
mod math_helpers;

use crate::attempt_01::*;
use crate::attempt_02::*;
use crate::attempt_03::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    attempt_01_naive()?;
    attempt_02_spectrum()?;
    attempt_03_refractance()?;
    Ok(())
}
