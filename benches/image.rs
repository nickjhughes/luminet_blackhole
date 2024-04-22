use cgmath::Deg;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use luminet_blackhole_lib::{plotting::generate_flux_image, BlackHole};

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("generate_flux_image width=256 samples=5000", |b| {
        b.iter(|| {
            let blackhole = BlackHole::default();
            generate_flux_image(
                &blackhole,
                black_box(Deg(80.0)),
                black_box(5000),
                black_box(256),
                black_box(135),
                None,
            )
            .unwrap();
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
