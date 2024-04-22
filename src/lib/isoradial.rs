use crate::{blackhole::BlackHole, solvers::calc_impact_parameter};
use cgmath::{Angle, Rad, Vector2};
use std::f64::consts::PI;

pub struct IsoRadial {
    /// Mass of the associated black hole.
    mass: f64,
    /// Radius of the isoradial line, in units of the black hole's mass.
    pub radius: f64,
    /// Order of the isoradial line.
    pub order: u32,
}

impl IsoRadial {
    #[must_use]
    pub fn new(blackhole: &BlackHole, radius: f64, order: u32) -> Self {
        IsoRadial {
            mass: blackhole.mass,
            radius,
            order,
        }
    }

    /// Calculate the coordinates of this isoradial line as it would appear to the observer.
    #[must_use]
    pub fn calculate_coordinates(
        &self,
        inclination: Rad<f64>,
        num_angles: usize,
    ) -> Vec<Vector2<f64>> {
        (0..num_angles)
            .map(|i| {
                let alpha = Rad((i as f64) / (num_angles as f64) * 2.0 * PI);
                let impact_parameter =
                    calc_impact_parameter(self.radius, inclination, alpha, self.mass, self.order);
                Vector2::new(
                    impact_parameter * alpha.cos(),
                    impact_parameter * alpha.sin(),
                )
            })
            .collect::<Vec<Vector2<f64>>>()
    }

    /// Calculate the impact parameter corresponding to the given angle on this isoradial curve.
    #[must_use]
    pub fn get_impact_parameter_from_alpha(&self, inclination: Rad<f64>, alpha: Rad<f64>) -> f64 {
        calc_impact_parameter(self.radius, inclination, alpha, self.mass, self.order)
    }
}
