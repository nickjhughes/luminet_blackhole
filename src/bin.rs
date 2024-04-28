use cgmath::{Deg, Rad};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Generate plots of isoradial curves.
    Isoradials {
        /// Viewer's inclination in degrees above the equatorial plane.
        #[arg(short, long, default_value_t = 80.0)]
        inclination: f64,

        /// Direct (order = 0) radii to plot.
        #[arg(long, default_values_t = vec![6.0, 10.0, 20.0, 30.0])]
        direct_radii: Vec<f64>,

        /// Ghost (order = 1) radii to plot.
        #[arg(long, default_values_t = vec![6.0, 10.0, 30.0, 10000.0])]
        ghost_radii: Vec<f64>,

        /// Black hole's accretion rate.
        #[arg(long, default_value_t = luminet_blackhole_lib::DEFAULT_ACCRETION_RATE)]
        accretion_rate: f64,

        /// Black hole's accretion disk outer edge.
        #[arg(long, default_value_t = luminet_blackhole_lib::DEFAULT_DISK_OUTER_EDGE)]
        disk_outer_edge: f64,

        /// Output file path.
        path: PathBuf,
    },

    /// Generate image of observered flux.
    Flux {
        /// Viewer's inclination in degrees above the equatorial plane.
        #[arg(short, long, default_value_t = 80.0)]
        inclination: f64,

        /// Number of flux samples (more = slower but better quality output).
        #[arg(short, long, default_value_t = 200_000)]
        samples: usize,

        /// Output image width in pixels.
        #[arg(long, default_value_t = 2048)]
        width: u32,

        /// Output image height in pixels.
        #[arg(long, default_value_t = 1080)]
        height: u32,

        /// Black hole's accretion rate.
        #[arg(long, default_value_t = luminet_blackhole_lib::DEFAULT_ACCRETION_RATE)]
        accretion_rate: f64,

        /// Black hole's accretion disk outer edge.
        #[arg(long, default_value_t = luminet_blackhole_lib::DEFAULT_DISK_OUTER_EDGE)]
        disk_outer_edge: f64,

        /// Output file path.
        path: PathBuf,
    },

    /// Generate series of images of observed flux from different inclinations.
    FluxRange {
        /// Start of inclination range.
        #[arg(long, default_value_t = 10.0)]
        start: f64,

        /// End of inclination range.
        #[arg(long, default_value_t = 80.0)]
        end: f64,

        /// Step size of inclination range.
        #[arg(long, default_value_t = 10.0)]
        step: f64,

        /// Number of flux samples (more = slower but better quality output).
        #[arg(short, long, default_value_t = 200_000)]
        samples: usize,

        /// Output image width in pixels.
        #[arg(long, default_value_t = 2048)]
        width: u32,

        /// Output image height in pixels.
        #[arg(long, default_value_t = 1080)]
        height: u32,

        /// Black hole's accretion rate.
        #[arg(long, default_value_t = luminet_blackhole_lib::DEFAULT_ACCRETION_RATE)]
        accretion_rate: f64,

        /// Black hole's accretion disk outer edge.
        #[arg(long, default_value_t = luminet_blackhole_lib::DEFAULT_DISK_OUTER_EDGE)]
        disk_outer_edge: f64,

        /// Output directory path.
        directory: PathBuf,

        /// Output filename prefix.
        filename_prefix: String,
    },

    /// Dither an image.
    Dither {
        /// The dither algorithm to use.
        #[arg(short, long, default_value_t = luminet_blackhole_lib::plotting::DitherAlgorithm::BlueNoise)]
        algorithm: luminet_blackhole_lib::plotting::DitherAlgorithm,

        /// Input image path.
        input_path: PathBuf,

        /// Output image path.
        output_path: PathBuf,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Command::Isoradials {
            inclination,
            direct_radii,
            ghost_radii,
            accretion_rate,
            disk_outer_edge,
            path,
        } => {
            let blackhole =
                luminet_blackhole_lib::BlackHole::new(1.0, accretion_rate, disk_outer_edge);
            let radii = direct_radii
                .iter()
                .map(|&r| (r, 0))
                .chain(ghost_radii.iter().map(|&r| (r, 1)))
                .collect::<Vec<(f64, u32)>>();
            luminet_blackhole_lib::plotting::plot_isoradials(
                &blackhole,
                Deg(inclination),
                &radii,
                path,
            )?;
        }
        Command::Flux {
            inclination,
            samples,
            width,
            height,
            accretion_rate,
            disk_outer_edge,
            path,
        } => {
            let blackhole =
                luminet_blackhole_lib::BlackHole::new(1.0, accretion_rate, disk_outer_edge);
            let img = luminet_blackhole_lib::plotting::generate_flux_image(
                &blackhole,
                Deg(inclination),
                samples,
                width,
                height,
                None,
            )?;
            img.save(path)?;
        }
        Command::FluxRange {
            start,
            end,
            step,
            samples,
            width,
            height,
            accretion_rate,
            disk_outer_edge,
            directory,
            filename_prefix,
        } => {
            assert!(directory.is_dir(), "`directory` must be a directory");

            let inclinations = {
                let mut inclinations = Vec::new();
                let mut i = start;
                while i <= end {
                    inclinations.push(Rad::from(Deg(i)));
                    i += step;
                }
                inclinations
            };
            let blackhole =
                luminet_blackhole_lib::BlackHole::new(1.0, accretion_rate, disk_outer_edge);
            let images = luminet_blackhole_lib::plotting::generate_flux_images_inclinations(
                &blackhole,
                samples,
                &inclinations,
                width,
                height,
            )?;
            for (inclination, img) in inclinations.iter().zip(images.iter()) {
                let filename = {
                    let mut f = filename_prefix.clone();
                    f.push_str(&format!("{:.0}.png", Deg::from(*inclination).0));
                    f
                };
                let path = directory.join(filename);
                img.save(path)?;
            }
        }
        Command::Dither {
            input_path,
            output_path,
            algorithm,
        } => {
            let dynamic_img = image::io::Reader::open(input_path)?.decode()?;
            let mut img = dynamic_img.to_luma16();
            luminet_blackhole_lib::plotting::dither(algorithm, &mut img);
            img.save(&output_path)?;
        }
    }

    Ok(())
}
