pub use dither::{dither, DitherAlgorithm};
pub use flux::{generate_flux_image, generate_flux_images_inclinations, Luma16Image};
pub use isoradial::plot_isoradials;

mod dither;
mod flux;
mod gilbert;
mod isoradial;
