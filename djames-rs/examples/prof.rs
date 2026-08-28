//! Where signing time goes: K-multiply cost vs. one root-finding call.
use djames::gf::Ext;
use djames::params;
use djames::poly;
use djames::symmetric::Xof;
use std::time::Instant;

fn main() {
    for name in ["d-james-128-q2", "d-james-128-q5", "d-james-128-q13"] {
        let p = params::by_name(name).unwrap();
        let k = Ext::new(p.q, p.n, p.fpoly);
        let mut x = Xof::new(&[b"prof"]);
        let mut rnd = || {
            let cs = djames::symmetric::sample_fq(&mut x, p.n, p.q);
            k.from_coords(&cs)
        };
        let (a, b) = (rnd(), rnd());
        let iters = 200000;
        let t = Instant::now();
        let mut acc = k.zero();
        for _ in 0..iters {
            acc = k.mul(&a, &b);
        }
        let per_mul = t.elapsed().as_nanos() as f64 / iters as f64;
        std::hint::black_box(&acc);

        // one root-finding on a random degree-D polynomial
        let mut f: poly::Poly = (0..=p.dd).map(|_| rnd()).collect();
        f[p.dd] = k.one();
        let mut ed = Xof::new(&[b"edf"]);
        let t = Instant::now();
        let reps = 20;
        for _ in 0..reps {
            std::hint::black_box(poly::roots(&k, &f, &mut ed));
        }
        let per_root = t.elapsed().as_secs_f64() * 1e3 / reps as f64;
        println!(
            "{name:<17} q={:<3} n={:<4} D={:<3} mul {:>7.0} ns   roots {:>8.2} ms  = {:>7.0} muls",
            p.q,
            p.n,
            p.dd,
            per_mul,
            per_root,
            per_root * 1e6 / per_mul
        );
    }
}
