use crate::{BlackHole, IsoRadial};
use cgmath::{Basis2, Deg, Rad, Rotation, Rotation2};
use plotters::prelude::*;
use std::f64::consts::PI;

const IMAGE_RESOLUTION: (u32, u32) = (1024, 1024);
const ANGLE_COUNT: usize = 360;

/// Plot a set of isoradial curves for the given black hole.
pub fn plot_isoradials<P: AsRef<std::path::Path>, A: Into<Rad<f64>>>(
    blackhole: &BlackHole,
    inclination: A,
    radii: &[(f64, u32)],
    path: P,
) -> Result<(), Box<dyn std::error::Error>> {
    let inclination: Rad<f64> = inclination.into();

    let root = BitMapBackend::new(&path, IMAGE_RESOLUTION).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart =
        ChartBuilder::on(&root).build_cartesian_2d(-35.0_f32..35.0_f32, -35.0_f32..35.0_f32)?;

    // Plot apparent black hole radius
    let angles = (0..u32::try_from(ANGLE_COUNT)?).map(|i| f64::from(i) / 360_f64 * 2.0 * PI);
    #[allow(clippy::cast_possible_truncation)]
    chart.draw_series(LineSeries::new(
        angles.clone().map(|a| {
            let apparent_inner_edge_impact_parameter = blackhole
                .apparent_inner_edge_radius(inclination, Rad(a + PI / 2.0))
                .min(blackhole.critical_impact_parameter());
            (
                (apparent_inner_edge_impact_parameter * a.cos()) as f32,
                (apparent_inner_edge_impact_parameter * a.sin()) as f32,
            )
        }),
        ShapeStyle {
            color: BLACK.mix(1.0),
            filled: false,
            stroke_width: 2,
        },
    ))?;

    // Plot isoradials
    let rotation = Basis2::from_angle(Deg(-90.0));
    for (radius, order) in radii {
        let isoradial = IsoRadial::new(blackhole, *radius, *order);
        let coords = isoradial.calculate_coordinates(inclination, ANGLE_COUNT);
        #[allow(clippy::cast_possible_truncation)]
        chart.draw_series(LineSeries::new(
            coords
                .iter()
                .map(|&pt| {
                    // Rotate points by -90 deg, and vertically flip ghost image points
                    let pt = rotation.rotate_vector(pt);
                    let y = if *order > 0 { -pt.y } else { pt.y };
                    (pt.x as f32, y as f32)
                })
                .collect::<Vec<(f32, f32)>>(),
            ShapeStyle {
                color: BLACK.mix(if *order > 0 { 0.25 } else { 0.5 }),
                filled: false,
                stroke_width: 2,
            },
        ))?;
    }

    root.present()?;
    Ok(())
}
