pub struct ColourSystem {
    pub name: &'static str,
    pub x_red: f64,
    pub y_red: f64,
    pub x_green: f64,
    pub y_green: f64,
    pub x_blue: f64,
    pub y_blue: f64,
    pub x_white: f64,
    pub y_white: f64,
    pub gamma: f64,
}

/* For NTSC television */
const ILLUMINANT_C: [f64; 2] = [0.3101, 0.3162];
/* For EBU and SMPTE */
const ILLUMINANT_D65: [f64; 2] = [0.3127, 0.3291];
/* CIE equal-energy illuminant */
const ILLUMINANT_E: [f64; 2] = [0.33333333, 0.33333333];

const GAMMA_REC709: f64 = 0.0;

#[allow(dead_code)]
pub const NSTC_SYSTEM: ColourSystem = ColourSystem {
    name: "NSTC",
    x_red: 0.67,
    y_red: 0.33,
    x_green: 0.21,
    y_green: 0.71,
    x_blue: 0.14,
    y_blue: 0.08,
    x_white: ILLUMINANT_C[0],
    y_white: ILLUMINANT_C[1],
    gamma: GAMMA_REC709,
};

#[allow(dead_code)]
pub const EBU_SYSTEM: ColourSystem = ColourSystem {
    name: "EBU (PAL/SECAM)",
    x_red: 0.64,
    y_red: 0.33,
    x_green: 0.29,
    y_green: 0.60,
    x_blue: 0.14,
    y_blue: 0.08,
    x_white: ILLUMINANT_C[0],
    y_white: ILLUMINANT_C[1],
    gamma: GAMMA_REC709,
};

#[allow(dead_code)]
pub const SMPTE_SYSTEM: ColourSystem = ColourSystem {
    name: "SMPTE",
    x_red: 0.630,
    y_red: 0.340,
    x_green: 0.310,
    y_green: 0.595,
    x_blue: 0.155,
    y_blue: 0.070,
    x_white: ILLUMINANT_D65[0],
    y_white: ILLUMINANT_D65[1],
    gamma: GAMMA_REC709,
};

#[allow(dead_code)]
pub const HDTV_SYSTEM: ColourSystem = ColourSystem {
    name: "HDTV",
    x_red: 0.670,
    y_red: 0.330,
    x_green: 0.210,
    y_green: 0.710,
    x_blue: 0.150,
    y_blue: 0.060,
    x_white: ILLUMINANT_D65[0],
    y_white: ILLUMINANT_D65[1],
    gamma: GAMMA_REC709,
};

#[allow(dead_code)]
pub const CIE_SYSTEM: ColourSystem = ColourSystem {
    name: "CIE",
    x_red: 0.7355,
    y_red: 0.2645,
    x_green: 0.2658,
    y_green: 0.7243,
    x_blue: 0.1669,
    y_blue: 0.0085,
    x_white: ILLUMINANT_E[0],
    y_white: ILLUMINANT_E[1],
    gamma: GAMMA_REC709,
};

#[allow(dead_code)]
pub const REC709_SYSTEM: ColourSystem = ColourSystem {
    name: "CIE REC 709",
    x_red: 0.64,
    y_red: 0.33,
    x_green: 0.30,
    y_green: 0.60,
    x_blue: 0.15,
    y_blue: 0.06,
    x_white: ILLUMINANT_D65[0],
    y_white: ILLUMINANT_D65[1],
    gamma: GAMMA_REC709,
};
