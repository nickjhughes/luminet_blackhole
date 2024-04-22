use cgmath::Deg;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use luminet_blackhole_lib::BlackHole;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("sample_flux_at_points num_points=1000 order=0", |b| {
        b.iter(|| {
            let blackhole = BlackHole::default();
            let _ = blackhole.sample_flux_at_points(
                black_box(Deg(80.0)),
                black_box(1000),
                black_box(0),
            );
        })
    });

    c.bench_function("sample_flux_at_points num_points=1000 order=1", |b| {
        b.iter(|| {
            let blackhole = BlackHole::default();
            let _ = blackhole.sample_flux_at_points(
                black_box(Deg(80.0)),
                black_box(1000),
                black_box(1),
            );
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
