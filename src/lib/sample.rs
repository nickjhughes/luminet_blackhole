use crate::BlackHole;
use cgmath::{Angle, Deg, Rad, Vector2};
use std::io::Write;

/// A sample of the observed flux from a black hole's accretion disk.
#[derive(Debug)]
pub struct Sample {
    /// The radius of the emitting photon's position in the black hole's frame.
    pub radius: f64,
    /// The angle of the emitting photon's position in the black hole and observer's frames.
    pub alpha: Rad<f64>,
    /// The radial location of the sample on the observer's photographic plate.
    pub impact_parameter: f64,
    /// The image order of the sample (0 = direct, 1+ = ghost).
    pub order: u32,
    /// The redshift factor `1 + z` of the sample.
    pub redshift_factor: f64,
    /// The observed flux `F_O` of the sample.
    pub observed_flux: f64,
}

impl Sample {
    /// Get the position of this sample in the black hole's reference frame.
    #[must_use]
    pub fn black_hole_position(&self) -> Vector2<f64> {
        Vector2::new(
            self.radius * self.alpha.cos(),
            self.radius * self.alpha.sin(),
        )
    }

    /// Get the position of this sample in the observer's reference frame.
    ///
    /// This will returned a flipped y-coordinate for order > 0 samples.
    #[must_use]
    pub fn observer_position(&self) -> Vector2<f64> {
        let y = self.impact_parameter * self.alpha.sin();
        Vector2::new(
            self.impact_parameter * self.alpha.cos(),
            if self.order > 0 { -y } else { y },
        )
    }
}

impl spade::HasPosition for &Sample {
    type Scalar = f64;

    fn position(&self) -> spade::Point2<f64> {
        let position = self.observer_position();
        spade::Point2 {
            x: position.x,
            y: position.y,
        }
    }
}

/// Save a number of flux samples to a CSV file.
#[allow(dead_code)]
pub fn save_samples<P: AsRef<std::path::Path>>(
    blackhole: &BlackHole,
    inclination: Rad<f64>,
    sample_count: usize,
    path: P,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut direct_samples = blackhole.sample_flux_at_points(inclination, sample_count, 0);
    let mut ghost_samples = blackhole.sample_flux_at_points(inclination, sample_count, 1);

    // Rotate points by -90 deg
    let rotation_angle = Rad::from(Deg(-90.0));
    for sample in &mut direct_samples {
        sample.alpha += rotation_angle;
    }
    for sample in &mut ghost_samples {
        sample.alpha += rotation_angle;
    }

    // Write samples to file
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "x,y,r,b,alpha,order,flux")?;
    for sample in direct_samples.iter().chain(ghost_samples.iter()) {
        let observer_point = sample.observer_position();
        writeln!(
            file,
            "{},{},{},{},{},{},{}",
            observer_point.x,
            observer_point.y,
            sample.radius,
            sample.impact_parameter,
            sample.alpha.0,
            sample.order,
            sample.observed_flux
        )?;
    }

    Ok(())
}
