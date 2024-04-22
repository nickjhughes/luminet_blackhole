# Luminet Blackhole

Recreating the simulated image of the observed flux from the black hole's accretion disk from ["Image of a spherical black hole with thin accretion disk.", Luminet, J. -P., AAP, 75, 228-235 (1979)](https://ui.adsabs.harvard.edu/abs/1979A%26A....75..228L/abstract).

![Observed flux from an inclination of 80 degrees](output/flux_80.png?raw=true)

[Observed flux from 1 deg to 89 deg of inclination and back](https://github.com/nickjhughes/luminet_blackhole/assets/689543/96a1afa9-d8cd-475e-888c-0a8befe68db6)

## Running

Run `cargo run --release -- --help` to see options.

To generate a 2K resolution image of the observed flux from an inclination of 80 degrees, run:
```bash
cargo run --release -- flux -i 80 flux_80.png
```

## Alternate Implementations

Huge thanks to these projects for identifying the errors in the paper's equations:

- https://github.com/bgmeulem/Luminet
- https://github.com/peytondmurray/bhsim
