use crate::{
    equations::{calc_observed_flux, calc_redshift_factor},
    solvers::calc_impact_parameter,
    IsoRadial, Sample,
};
use cgmath::Rad;
use rand::{distributions::Uniform, prelude::*};
use rayon::prelude::*;
use std::f64::consts::PI;

pub const DEFAULT_ACCRETION_RATE: f64 = 10e-8;
pub const DEFAULT_DISK_OUTER_EDGE: f64 = 50.0;

/// A black hole with with a thin accretion disk.
pub struct BlackHole {
    /// Black hole mass.
    pub mass: f64,
    /// Accretion rate.
    pub accretion_rate: f64,
    /// The outer edge of the accretion disk, in units of black hole mass.
    disk_outer_edge: f64,
}

impl Default for BlackHole {
    fn default() -> Self {
        Self {
            mass: 1.0,
            accretion_rate: DEFAULT_ACCRETION_RATE,
            disk_outer_edge: DEFAULT_DISK_OUTER_EDGE,
        }
    }
}

impl BlackHole {
    #[must_use]
    pub fn new(mass: f64, accretion_rate: f64, disk_outer_edge: f64) -> Self {
        BlackHole {
            mass,
            accretion_rate,
            disk_outer_edge,
        }
    }

    /// Value of the critical impact parameter for this black hole.
    #[must_use]
    pub fn critical_impact_parameter(&self) -> f64 {
        3.0 * 3.0_f64.sqrt() * self.mass
    }

    /// The radius of the outer edge of the accretion disk.
    #[must_use]
    pub fn disk_outer_edge(&self) -> f64 {
        self.disk_outer_edge * self.mass
    }

    /// The radius of the inner edge of the accretion disk.
    #[must_use]
    pub fn disk_inner_edge(&self) -> f64 {
        6.0 * self.mass
    }

    /// Construct an isoradial forming the apparent inner edge of the accretion disk.
    #[must_use]
    pub fn apparent_inner_disk_edge(&self) -> IsoRadial {
        IsoRadial::new(self, self.disk_inner_edge(), 0)
    }

    /// Construct an isoradial forming the apparent outer edge of the accretion disk.
    #[must_use]
    pub fn apparent_outer_disk_edge(&self) -> IsoRadial {
        IsoRadial::new(self, self.disk_outer_edge(), 0)
    }

    /// Calculate the apparent outer edge radius of the black hole at the given angle.
    #[must_use]
    pub fn apparent_outer_edge_radius(&self, inclination: Rad<f64>, alpha: Rad<f64>) -> f64 {
        self.apparent_outer_disk_edge()
            .get_impact_parameter_from_alpha(inclination, alpha)
    }

    /// Calculate the apparent inner edge radius of the black hole at the given angle.
    #[must_use]
    pub fn apparent_inner_edge_radius(&self, inclination: Rad<f64>, alpha: Rad<f64>) -> f64 {
        self.apparent_inner_disk_edge()
            .get_impact_parameter_from_alpha(inclination, alpha)
    }

    /// Sample the observed flux from the accretion disk at a number of random points.
    #[must_use]
    pub fn sample_flux_at_points<A: Into<Rad<f64>>>(
        &self,
        inclination: A,
        num_points: usize,
        order: u32,
    ) -> Vec<Sample> {
        let inclination: Rad<f64> = inclination.into();

        let radius_dist = Uniform::new(self.disk_inner_edge(), self.disk_outer_edge());
        let alpha_dist = Uniform::new(0.0, 2.0 * PI);

        (0..num_points)
            .into_par_iter()
            .map_init(rand::thread_rng, |rng, _| {
                let radius = rng.sample(radius_dist);
                let alpha = Rad(rng.sample(alpha_dist));

                let impact_parameter =
                    calc_impact_parameter(radius, inclination, alpha, self.mass, order);
                let redshift_factor =
                    calc_redshift_factor(radius, alpha, inclination, self.mass, impact_parameter);
                let observed_flux =
                    calc_observed_flux(radius, self.accretion_rate, self.mass, redshift_factor);

                Sample {
                    radius,
                    alpha,
                    impact_parameter,
                    order,
                    redshift_factor,
                    observed_flux,
                }
            })
            .collect::<Vec<Sample>>()
    }
}
