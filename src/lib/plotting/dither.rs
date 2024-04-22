use rayon::iter::ParallelIterator;

pub fn floyd(img: &mut image::GrayImage) {
    let m = [
        (1, 7),
        (img.width() as usize - 1, 3),
        (img.width() as usize, 5),
        (img.width() as usize + 1, 1),
    ];
    for i in 0..img.len() {
        let x = f64::from(*img.get(i).unwrap()) / 255.0;
        let col = if x > 0.5 { 1.0 } else { 0.0 };
        let err = (x - col) / 16.0;
        for (x, y) in &m {
            if let Some(pixel) = img.get_mut(i + x) {
                *pixel += (err * f64::from(*y) * 255.0).round() as u8;
            }
        }
        *img.get_mut(i).unwrap() = (col * 255.0).round() as u8;
    }
}

pub fn atkinson(img: &mut image::GrayImage) {
    let m = [
        1,
        2,
        img.width() as usize - 1,
        img.width() as usize,
        img.width() as usize + 1,
        (img.width() * 2) as usize,
    ];
    for i in 0..img.len() {
        let x = f64::from(*img.get(i).unwrap()) / 255.0;
        let col = if x > 0.5 { 1.0 } else { 0.0 };
        let err = (x - col) / 8.0;
        for x in &m {
            if let Some(pixel) = img.get_mut(i + x) {
                *pixel += (err * 255.0).round() as u8;
            }
        }
        *img.get_mut(i).unwrap() = (col * 255.0).round() as u8;
    }
}

pub fn blue_noise(img: &mut image::GrayImage) {
    let mask = image::io::Reader::open("images/blue_noise.png")
        .expect("can open blue noise image file")
        .decode()
        .expect("can decode blue noise image file")
        .to_luma8();
    img.par_enumerate_pixels_mut().for_each(|(x, y, pixel)| {
        let mask_value = mask.get_pixel(x % mask.width(), y % mask.height()).0[0];
        let pixel_value = pixel.0[0];
        if pixel_value > mask_value {
            pixel.0[0] = 255;
        } else {
            pixel.0[0] = 0;
        }
    });
}
