use crate::{BlackHole, Sample};
use cgmath::{Deg, Rad, Vector2};
use image::Luma;
use indicatif::{ParallelProgressIterator, ProgressBar};
use rayon::iter::ParallelIterator;
use spade::{Barycentric, DelaunayTriangulation, FloatTriangulation, Triangulation};
use std::{cmp::Ordering, f64::consts::PI, ops::RangeInclusive};

pub type Luma16Image = image::ImageBuffer<Luma<u16>, Vec<u16>>;

/// Image order to show at an image pixel.
enum OrderToShow {
    None,
    Direct,
    Ghost,
}

/// Generate a series of images with the given viewer inclination.
///
/// The flux values will be normalized across the whole series of images.
pub fn generate_flux_images_inclinations(
    blackhole: &BlackHole,
    sample_count: usize,
    inclinations: &[Rad<f64>],
    image_width: u32,
    image_height: u32,
) -> Result<Vec<Luma16Image>, Box<dyn std::error::Error>> {
    // Sample flux at all inclinations
    let mut all_samples = Vec::new();
    let mut max_flux = 0.0;
    for &inclination in inclinations {
        let direct_samples = blackhole.sample_flux_at_points(inclination, sample_count, 0);
        let ghost_samples = blackhole.sample_flux_at_points(inclination, sample_count, 1);

        let inclination_max_flux = direct_samples
            .iter()
            .chain(ghost_samples.iter())
            .map(|s| s.observed_flux)
            .max_by(|a, b| a.partial_cmp(b).expect("no NaNs"))
            .expect("non-empty iter of samples");
        if inclination_max_flux > max_flux {
            max_flux = inclination_max_flux;
        }

        all_samples.push((direct_samples, ghost_samples));
    }
    let flux_range = 0.0..=max_flux;

    let mut images = Vec::new();
    for (&inclination, (direct_samples, ghost_samples)) in
        inclinations.iter().zip(all_samples.iter_mut())
    {
        images.push(generate_flux_image_from_samples(
            blackhole,
            inclination,
            direct_samples,
            ghost_samples,
            image_width,
            image_height,
            Some(flux_range.clone()),
        )?);
    }
    Ok(images)
}

/// Generate an image of the observed flux.
pub fn generate_flux_image<A: Into<Rad<f64>>>(
    blackhole: &BlackHole,
    inclination: A,
    sample_count: usize,
    image_width: u32,
    image_height: u32,
    flux_range: Option<RangeInclusive<f64>>,
) -> Result<Luma16Image, Box<dyn std::error::Error>> {
    let inclination: Rad<f64> = inclination.into();
    let mut direct_samples = blackhole.sample_flux_at_points(inclination, sample_count, 0);
    let mut ghost_samples = blackhole.sample_flux_at_points(inclination, sample_count, 1);
    generate_flux_image_from_samples(
        blackhole,
        inclination,
        &mut direct_samples,
        &mut ghost_samples,
        image_width,
        image_height,
        flux_range,
    )
}

/// Generate an image of the observed flux using the supplied samples.
pub fn generate_flux_image_from_samples(
    blackhole: &BlackHole,
    inclination: Rad<f64>,
    direct_samples: &mut [Sample],
    ghost_samples: &mut [Sample],
    image_width: u32,
    image_height: u32,
    flux_range: Option<RangeInclusive<f64>>,
) -> Result<Luma16Image, Box<dyn std::error::Error>> {
    // Rotate points by -90 deg
    let rotation_angle = Rad::from(Deg(-90.0));
    for sample in direct_samples.iter_mut() {
        sample.alpha += rotation_angle;
    }
    for sample in ghost_samples.iter_mut() {
        sample.alpha += rotation_angle;
    }

    let (min_point, max_point) = samples_range(direct_samples.iter().chain(ghost_samples.iter()));
    let flux_range = flux_range.unwrap_or_else(|| {
        let flux_max = direct_samples
            .iter()
            .chain(ghost_samples.iter())
            .map(|s| s.observed_flux)
            .max_by(|a, b| a.partial_cmp(b).expect("no NaNs"))
            .expect("non-empty iter of samples");
        0.0..=flux_max
    });

    // Fit the sampled region into the full width of the image, and use that to define the number
    // of real-world units per pixel (to make sure the aspect ratio of each pixel is equal)
    let units_per_pixel = (max_point.x - min_point.x) / f64::from(image_width);
    let mut img = Luma16Image::new(image_width, image_height);

    // Create Delaunay triangulation so we can linearly interpolate samples on the image pixel grid
    let direct_triangulation = {
        let mut t: DelaunayTriangulation<&Sample> = DelaunayTriangulation::new();
        for sample in direct_samples.iter() {
            t.insert(sample)?;
        }
        t
    };
    let ghost_triangulation = {
        let mut t: DelaunayTriangulation<&Sample> = DelaunayTriangulation::new();
        for sample in ghost_samples.iter() {
            t.insert(sample)?;
        }
        t
    };

    let progress_bar_style = indicatif::ProgressStyle::with_template(
        "{prefix} {bar:60.cyan/blue} {pos:>7}/{len:7} pixels",
    )
    .unwrap();
    let progress_bar = ProgressBar::new(img.len() as u64)
        .with_prefix("Rendering image...")
        .with_style(progress_bar_style);
    img.par_enumerate_pixels_mut()
        .progress_with(progress_bar)
        .for_each_init(
            || {
                (
                    direct_triangulation.barycentric(),
                    ghost_triangulation.barycentric(),
                )
            },
            |(direct_interpolater, ghost_interpolator), (col, row, pixel)| {
                let x = f64::from((col as i32) - ((image_width / 2) as i32)) * units_per_pixel;
                let y = -f64::from((row as i32) - ((image_height / 2) as i32)) * units_per_pixel;

                // Determine which zone we're in:
                //   - Outside the apparent outer edge of the accretion disk -> show ghost image
                //   - Inside the apparent inner edge of the accretion disk -> show ghost image
                //   - Inside the apparent inner edge of the black hole -> set to black
                //   - Otherwise -> show direct image
                let impact_parameter = (x.powi(2) + y.powi(2)).sqrt();
                let alpha = Rad(y.atan2(x) + PI / 2.0);
                let order_to_show = if impact_parameter
                    <= blackhole.apparent_inner_edge_radius(inclination, alpha)
                    || impact_parameter > blackhole.apparent_outer_edge_radius(inclination, alpha)
                {
                    let apparent_inner_edge_impact_parameter = {
                        blackhole
                            .apparent_inner_edge_radius(inclination, alpha)
                            .min(blackhole.critical_impact_parameter())
                    };
                    if impact_parameter < apparent_inner_edge_impact_parameter {
                        OrderToShow::None
                    } else {
                        OrderToShow::Ghost
                    }
                } else {
                    OrderToShow::Direct
                };

                let point = spade::Point2 { x, y };
                match order_to_show {
                    OrderToShow::None => {
                        *pixel = image::Luma([0]);
                    }
                    OrderToShow::Direct => {
                        let flux = interpolate_and_normalize_flux(
                            &point,
                            direct_interpolater,
                            &flux_range,
                        );
                        #[allow(clippy::cast_possible_truncation)]
                        let luma = (flux * f64::from(u16::MAX)).round() as u16;
                        *pixel = image::Luma([luma]);
                    }
                    OrderToShow::Ghost => {
                        let flux =
                            interpolate_and_normalize_flux(&point, ghost_interpolator, &flux_range);
                        #[allow(clippy::cast_possible_truncation)]
                        let luma = (flux * f64::from(u16::MAX)).round() as u16;
                        *pixel = image::Luma([luma]);
                    }
                }
            },
        );

    Ok(img)
}

fn interpolate_and_normalize_flux(
    point: &spade::Point2<f64>,
    interpolator: &mut Barycentric<'_, DelaunayTriangulation<&Sample>>,
    flux_range: &RangeInclusive<f64>,
) -> f64 {
    if let Some(flux) = interpolator.interpolate(|v| v.data().observed_flux, *point) {
        (flux - flux_range.start()) / (flux_range.end() - flux_range.start())
    } else {
        0.0
    }
}

fn samples_range<'a, I>(samples: I) -> (Vector2<f64>, Vector2<f64>)
where
    I: Iterator<Item = &'a Sample>,
{
    let points = samples.map(Sample::observer_position);
    let (min_x, min_y, max_x, max_y) = points.fold(
        (f64::MAX, f64::MAX, f64::MIN, f64::MIN),
        |(mut min_x, mut min_y, mut max_x, mut max_y), pt| {
            if pt.x.partial_cmp(&min_x) == Some(Ordering::Less) {
                min_x = pt.x;
            }
            if pt.x.partial_cmp(&max_x) == Some(Ordering::Greater) {
                max_x = pt.x;
            }
            if pt.y.partial_cmp(&min_y) == Some(Ordering::Less) {
                min_y = pt.y;
            }
            if pt.y.partial_cmp(&max_y) == Some(Ordering::Greater) {
                max_y = pt.y;
            }
            (min_x, min_y, max_x, max_y)
        },
    );
    (Vector2::new(min_x, min_y), Vector2::new(max_x, max_y))
}

#[cfg(test)]
mod tests {
    use super::samples_range;
    use crate::Sample;
    use cgmath::{assert_abs_diff_eq, Deg, Rad, Vector2};

    #[test]
    fn test_samples_range() {
        {
            let samples = vec![Sample {
                radius: 1.0,
                alpha: Rad(0.0),
                impact_parameter: 1.0,
                order: 0,
                redshift_factor: 0.0,
                observed_flux: 0.0,
            }];
            let (min_pt, max_pt) = samples_range(samples.iter());
            assert_eq!(min_pt, Vector2::new(1.0, 0.0));
            assert_eq!(max_pt, Vector2::new(1.0, 0.0));
        }

        {
            let samples = vec![
                Sample {
                    radius: 1.0,
                    alpha: Rad(0.0),
                    impact_parameter: 1.0,
                    order: 0,
                    redshift_factor: 0.0,
                    observed_flux: 0.0,
                },
                Sample {
                    radius: 1.0,
                    alpha: Rad::from(Deg(-180.0)),
                    impact_parameter: 1.0,
                    order: 0,
                    redshift_factor: 0.0,
                    observed_flux: 0.0,
                },
            ];
            let (min_pt, max_pt) = samples_range(samples.iter());
            assert_abs_diff_eq!(min_pt, Vector2::new(-1.0, 0.0));
            assert_abs_diff_eq!(max_pt, Vector2::new(1.0, 0.0));
        }
    }
}
