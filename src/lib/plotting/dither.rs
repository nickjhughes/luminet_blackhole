use super::{gilbert, Luma16Image};
use clap::ValueEnum;
use rand::Rng;
use rayon::iter::ParallelIterator;

#[derive(Debug, Copy, Clone, ValueEnum)]
pub enum DitherAlgorithm {
    FloydSteinberg,
    Atkinson,
    BlueNoise,
    Random,
    Riemersma,
}

impl std::fmt::Display for DitherAlgorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DitherAlgorithm::FloydSteinberg => write!(f, "floyd-steinberg"),
            DitherAlgorithm::Atkinson => write!(f, "atkinson"),
            DitherAlgorithm::BlueNoise => write!(f, "blue-noise"),
            DitherAlgorithm::Random => write!(f, "random"),
            DitherAlgorithm::Riemersma => write!(f, "riemersma"),
        }
    }
}

pub fn dither(algorithm: DitherAlgorithm, img: &mut Luma16Image) {
    match algorithm {
        DitherAlgorithm::FloydSteinberg => floyd_steinberg(img),
        DitherAlgorithm::Atkinson => atkinson(img),
        DitherAlgorithm::BlueNoise => blue_noise(img),
        DitherAlgorithm::Random => random(img),
        DitherAlgorithm::Riemersma => riemersma(img),
    }
}

fn floyd_steinberg(img: &mut Luma16Image) {
    let m = [
        (1, 7),
        (img.width() as usize - 1, 3),
        (img.width() as usize, 5),
        (img.width() as usize + 1, 1),
    ];
    for i in 0..img.len() {
        let x = f64::from(*img.get(i).unwrap()) / f64::from(u16::MAX);
        let col = if x > 0.5 { 1.0 } else { 0.0 };
        let err = (x - col) / 16.0;
        for (x, y) in &m {
            if let Some(pixel) = img.get_mut(i + x) {
                *pixel += (err * f64::from(*y) * f64::from(u16::MAX)).round() as u16;
            }
        }
        *img.get_mut(i).unwrap() = (col * f64::from(u16::MAX)).round() as u16;
    }
}

fn atkinson(img: &mut Luma16Image) {
    let m = [
        1,
        2,
        img.width() as usize - 1,
        img.width() as usize,
        img.width() as usize + 1,
        (img.width() * 2) as usize,
    ];
    for i in 0..img.len() {
        let x = f64::from(*img.get(i).unwrap()) / f64::from(u16::MAX);
        let col = if x > 0.5 { 1.0 } else { 0.0 };
        let err = (x - col) / 8.0;
        for x in &m {
            if let Some(pixel) = img.get_mut(i + x) {
                *pixel += (err * f64::from(u16::MAX)).round() as u16;
            }
        }
        *img.get_mut(i).unwrap() = (col * f64::from(u16::MAX)).round() as u16;
    }
}

fn blue_noise(img: &mut Luma16Image) {
    let mask = image::io::Reader::open("images/blue_noise.png")
        .expect("can open blue noise image file")
        .decode()
        .expect("can decode blue noise image file")
        .to_luma16();
    img.par_enumerate_pixels_mut().for_each(|(x, y, pixel)| {
        let mask_value = mask.get_pixel(x % mask.width(), y % mask.height()).0[0];
        let pixel_value = pixel.0[0];
        if pixel_value > mask_value {
            pixel.0[0] = u16::MAX;
        } else {
            pixel.0[0] = 0;
        }
    });
}

fn random(img: &mut Luma16Image) {
    img.par_pixels_mut()
        .for_each_init(rand::thread_rng, |rng, pixel| {
            let offset: i16 = rng.gen();
            let pixel_value = pixel.0[0].saturating_add_signed(offset);
            if pixel_value > u16::MAX / 2 {
                pixel.0[0] = u16::MAX;
            } else {
                pixel.0[0] = 0;
            }
        });
}

fn riemersma(img: &mut Luma16Image) {
    const ERROR_FALLOFF: f64 = 1.0 / 4.0;
    const ERROR_LENGTH: usize = 32;

    let weights = (0..ERROR_LENGTH)
        .rev()
        .map(|i| ERROR_FALLOFF.powf((i as f64) / ((ERROR_LENGTH as f64) - 1.0)))
        .collect::<Vec<f64>>();
    let mut errors = vec![0_i32; ERROR_LENGTH];

    let mut curve_idx = 0;
    while curve_idx < img.len() {
        let (x, y) = gilbert::gilbert_d2xy(curve_idx as u32, img.width(), img.height());
        let adjustment = errors
            .iter()
            .zip(weights.iter())
            .map(|(e, w)| ((*e as f64) * w).round() as i32)
            .sum::<i32>();

        let source_pixel = img.get_pixel(x, y).0[0];
        let adjusted_pixel = (source_pixel as u32).saturating_add_signed(adjustment);
        let quantized_pixel = if adjusted_pixel > (u16::MAX / 2) as u32 {
            u16::MAX
        } else {
            0
        };
        let pixel = image::Luma([quantized_pixel as u16]);

        let error = (source_pixel as i32) - (quantized_pixel as i32);
        errors.push(error);
        errors.remove(0);

        img.put_pixel(x, y, pixel);
        curve_idx += 1;
    }
}
