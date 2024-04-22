//! Equations and definitions from the paper.
//!
//! Note that several equation have errors in the paper. These are noted in the relevant
//! function documentation.

use cgmath::{Angle, Rad};
use spec_math::Ellip;
use std::f64::consts::PI;

const INCLINATION_TOLERANCE: Rad<f64> = Rad(1e-5);

/// Calculate `Q` from the periastron `P` (pg 229).
pub fn calc_q(periastron: f64, mass: f64) -> f64 {
    ((periastron - 2.0 * mass) * (periastron + 6.0 * mass)).sqrt()
}

/// Calculate impact parameter `b` from the periastron `P` (eqn 5).
///
/// Note the equation in the paper has an error, it should be `b^2` on the left-hand side,
/// not `b`.
pub fn calc_impact_parameter_from_periastron(periastron: f64, mass: f64) -> f64 {
    (periastron.powi(3) / (periastron - 2.0 * mass)).sqrt()
}

/// Calculate the modulus `k^2` of the elliptic integral (eqn 12).
///
/// While equation 12 in the paper is correct, the definition of `k` on page 229 has an error,
/// the numerator should be in parentheses.
pub fn calc_modulus(periastron: f64, mass: f64, q: Option<f64>) -> f64 {
    let q = q.unwrap_or_else(|| calc_q(periastron, mass));
    (q - periastron + 6.0 * mass) / (2.0 * q)
}

/// Calculate `zeta_inf` for the elliptic integral (eqn 12).
pub fn calc_zeta_inf(periastron: f64, mass: f64, q: Option<f64>) -> f64 {
    let q = q.unwrap_or_else(|| calc_q(periastron, mass));
    ((q - periastron + 2.0 * mass) / (q - periastron + 6.0 * mass))
        .sqrt()
        .asin()
}

/// Calculate the cosine of angle `gamma` (eqn 10).
pub fn calc_cos_gamma(alpha: Rad<f64>, inclination: Rad<f64>) -> f64 {
    if inclination < INCLINATION_TOLERANCE {
        return 0.0;
    }
    alpha.cos() / (alpha.cos().powi(2) + 1.0 / inclination.tan().powi(2)).sqrt()
}

/// Calculate the cosine of the angle in the observer's reference frame `alpha`, from an angle
/// `phi` in the black hole's reference frame (eqn 9).
#[allow(dead_code)]
pub fn calc_cos_alpha(phi: f64, inclination: Rad<f64>) -> f64 {
    (phi.cos() * inclination.cos()) / (1.0 - inclination.sin().powi(2) * phi.cos().powi(2)).sqrt()
}

/// Calculate the reciprocal of `r` (eqn 13).
///
/// Note that the paper has an error in this equation, the `sqrt(P/Q)` factor in the first term
/// in the argument of the elliptic sine `sn` should be in the denominator, not the numerator.
pub fn calc_one_over_radius(
    periastron: f64,
    alpha: Rad<f64>,
    mass: f64,
    inclination: Rad<f64>,
    order: u32,
) -> f64 {
    let q = calc_q(periastron, mass);
    let zeta_inf = calc_zeta_inf(periastron, mass, Some(q));
    let modulus = calc_modulus(periastron, mass, Some(q));
    let elliptic_inf = zeta_inf.ellip_k_inc(modulus);
    let gamma = calc_cos_gamma(alpha, inclination).acos();

    let jacobian_elliptic_arg = if order == 0 {
        // Direct image
        gamma / (2.0 * (periastron / q).sqrt()) + elliptic_inf
    } else {
        // Higher-order image
        let elliptic_k = modulus.ellip_k();
        (gamma - 2.0 * (order as f64) * PI) / (2.0 * (periastron / q).sqrt()) - elliptic_inf
            + 2.0 * elliptic_k
    };
    let elliptic_sine = jacobian_elliptic_arg.ellip_j(modulus).sn;

    -(q - periastron + 2.0 * mass) / (4.0 * mass * periastron)
        + ((q - periastron + 6.0 * mass) / (4.0 * mass * periastron)) * elliptic_sine.powi(2)
}

/// Calculate the difference between the equation for `1/r` (eqn 13) and the reciprocal of the
/// given radius value. This value should be small in magnitude if the supplied periastron value
/// is the actual periastron value for the emitted photon.
///
/// Used for solving for the value of the periastron for given radius and alpha values.
pub fn calc_one_over_radius_minus_one_over_radius(
    radius: f64,
    periastron: f64,
    alpha: Rad<f64>,
    mass: f64,
    inclination: Rad<f64>,
    order: u32,
) -> f64 {
    1.0 - radius * calc_one_over_radius(periastron, alpha, mass, inclination, order)
}

/// Calculate the intrinsic flux of the disk `F_s` (eqn 15).
pub fn calc_intrinsic_flux(radius: f64, accretion_rate: f64, mass: f64) -> f64 {
    let radius_star = radius / mass;
    let log_arg = ((radius_star.sqrt() + 3.0_f64.sqrt()) * (6.0_f64.sqrt() - 3.0_f64.sqrt()))
        / ((radius_star.sqrt() - 3.0_f64.sqrt()) * (6.0_f64.sqrt() + 3.0_f64.sqrt()));
    ((3.0 * mass * accretion_rate) / (8.0 * PI))
        * (1.0 / ((radius_star - 3.0) * radius_star.powf(5.0 / 2.0)))
        * (radius_star.sqrt() - 6.0_f64.sqrt() + (3.0_f64.sqrt() / 3.0) * log_arg.log10())
}

/// Calculate the observed flux `F_O` (pg 233).
pub fn calc_observed_flux(
    radius: f64,
    accretion_rate: f64,
    mass: f64,
    redshift_factor: f64,
) -> f64 {
    calc_intrinsic_flux(radius, accretion_rate, mass) / redshift_factor.powi(4)
}

/// Calculate the gravitational redshift factor `1 + z`, ignoring cosmological redshift (eqn 19).
///
/// Note that while equation 19 is correct, the unlabelled but presumed equation 18 above is missing
/// several terms, it should read `1 + z = (1 - Ω*b*cos(η)) * (-g_tt -2*Ω*g_tϕ - Ω²*g_ϕϕ)^(-1/2)`.
pub fn calc_redshift_factor(
    radius: f64,
    alpha: Rad<f64>,
    inclination: Rad<f64>,
    mass: f64,
    impact_parameter: f64,
) -> f64 {
    (1.0 + (mass / radius.powi(3)).sqrt() * impact_parameter * inclination.sin() * alpha.sin())
        / (1.0 - 3.0 * mass / radius).sqrt()
}

/// The equation of an ellipse based on `cos(gamma)`.
///
/// Used as a fallback for when a periastron solution cannot be found, as isoradials form
/// ellipses in the Newtonian limit.
pub fn ellipse(radius: f64, alpha: Rad<f64>, inclination: Rad<f64>) -> f64 {
    let gamma = calc_cos_gamma(alpha, inclination).acos();
    radius * gamma.sin()
}
