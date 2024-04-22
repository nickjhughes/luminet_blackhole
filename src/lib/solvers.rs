//! Functions for solving for the periastron and impact parameter of a photon emitted
//! from the black hole's accretion disk.

use crate::equations::{
    calc_impact_parameter_from_periastron, calc_one_over_radius_minus_one_over_radius, ellipse,
};
use cgmath::Rad;

/// Solution tolerance to use when solving for the periastron.
const PERIASTRON_TOLERANCE: f64 = 1e-6;
/// The maximum number of iteration of the bisection method to run.
const MAX_BISECTION_ITERS: usize = 100;
/// The minumum periastron value to solve for, in units of black hole mass.
const MIN_PERIASTRON: f64 = 3.001;
/// The maximum periastron value to solve for, in units of black hole radius.
const MAX_PERIASTRON: f64 = 3.0;

/// For a given black hole reference frame `radius` and angle in the observer's frame `alpha`,
/// calculate the periastron for a photon emitted at `radius`.
///
/// This is done by finding a zero of the function `1.0 / calc_one_over_radius - radius` in terms
/// the periastron via the bisection method. Will fail and return None if no solution can be found
/// in the range [`MIN_PERIASTRON` * mass, `MAX_PERIASTRON` * radius].
pub fn calc_periastron(
    radius: f64,
    inclination: Rad<f64>,
    alpha: Rad<f64>,
    mass: f64,
    order: u32,
) -> Option<f64> {
    let periastron_range = (MIN_PERIASTRON * mass)..=(MAX_PERIASTRON * radius);

    // First determine if a solution exists in the valid range
    let val_at_min_periastron = calc_one_over_radius_minus_one_over_radius(
        radius,
        *periastron_range.start(),
        alpha,
        mass,
        inclination,
        order,
    );
    let val_at_max_periastron = calc_one_over_radius_minus_one_over_radius(
        radius,
        *periastron_range.end(),
        alpha,
        mass,
        inclination,
        order,
    );
    if val_at_min_periastron.signum() == val_at_max_periastron.signum() {
        // No solution in the valid range
        return None;
    }

    // Use the bisection to iteratively improve the solution
    let mut periastron_a = *periastron_range.start();
    let mut val_at_periastron_a = val_at_min_periastron;
    let mut periastron_b = *periastron_range.end();
    let mut val_at_periastron_b = val_at_max_periastron;
    debug_assert!(val_at_periastron_a.signum() != val_at_periastron_b.signum());
    let mut iter_count = 0;
    while (periastron_b - periastron_a).abs() > PERIASTRON_TOLERANCE
        && iter_count < MAX_BISECTION_ITERS
    {
        let periastron_c = (periastron_a + periastron_b) / 2.0;
        let val_at_periastron_c = calc_one_over_radius_minus_one_over_radius(
            radius,
            periastron_c,
            alpha,
            mass,
            inclination,
            order,
        );
        if val_at_periastron_a.signum() != val_at_periastron_c.signum() {
            val_at_periastron_b = val_at_periastron_c;
            periastron_b = periastron_c;
        } else if val_at_periastron_b.signum() != val_at_periastron_c.signum() {
            val_at_periastron_a = val_at_periastron_c;
            periastron_a = periastron_c;
        }
        iter_count += 1;
    }

    let result = (periastron_a + periastron_b) / 2.0;
    if result.is_nan() {
        None
    } else {
        Some(result)
    }
}

/// For a given black hole reference frame `radius` and angle in the observer's frame `alpha`,
/// calculate the impact parameter for a photon emitted at `radius`.
///
/// First solves for the photon's periastron using `calc_periastron`, then converts that to a value
/// for the impact parameter. If no solution for the periastron can be found, it will fallback
/// to using the equation for an ellipse.
pub fn calc_impact_parameter(
    radius: f64,
    inclination: Rad<f64>,
    alpha: Rad<f64>,
    mass: f64,
    order: u32,
) -> f64 {
    if let Some(periastron) = calc_periastron(radius, inclination, alpha, mass, order) {
        calc_impact_parameter_from_periastron(periastron, mass)
    } else {
        ellipse(radius, alpha, inclination)
    }
}
