//! Timings for keygen / sign / verify. `cargo run --release --example bench`.

use djames::{keygen, params, sign, verify};
use std::time::Instant;

fn main() {
    let which: Vec<&str> = std::env::args()
        .skip(1)
        .collect::<Vec<_>>()
        .iter()
        .map(|s| Box::leak(s.clone().into_boxed_str()) as &str)
        .collect();
    let sets: Vec<&'static params::Params> = if which.is_empty() {
        [
            "d-james-128-q2",
            "james-128-q2",
            "d-james-256-q2",
            "james-256-q2",
            "d-james-128-q4",
            "d-james-128-q5",
            "d-james-128-q13",
            "d-james-128-q23",
        ]
        .iter()
        .filter_map(|n| params::by_name(n))
        .collect()
    } else {
        which.iter().filter_map(|n| params::by_name(n)).collect()
    };
    println!(
        "{:<18} {:>4} {:>5} {:>10} {:>10} {:>10} {:>8} {:>12}",
        "set", "q", "n", "keygen", "sign", "verify", "sig B", "pk B"
    );
    for p in sets {
        let t = Instant::now();
        let (pk, sk) = keygen(p, &[7u8; 32]).unwrap();
        let tk = t.elapsed();
        let raw = pk.to_bytes();
        let t = Instant::now();
        let n = 5;
        let mut sig = Vec::new();
        for i in 0..n {
            sig = sign(p, &sk, format!("message {i}").as_bytes()).unwrap();
        }
        let ts = t.elapsed() / n;
        let t = Instant::now();
        for _ in 0..n {
            assert!(verify(
                p,
                &pk,
                format!("message {}", n - 1).as_bytes(),
                &sig
            ));
        }
        let tv = t.elapsed() / n;
        println!(
            "{:<18} {:>4} {:>5} {:>9.1?} {:>9.1?} {:>9.1?} {:>8} {:>12}",
            p.name,
            p.q,
            p.n,
            tk,
            ts,
            tv,
            sig.len(),
            raw.len()
        );
    }
}
