//! Parameter-table consistency, checked against `d-james-spec.md`.

use djames::params::{by_name, ALL};

#[test]
fn internally_consistent() {
    for p in ALL {
        assert_eq!(p.m, p.n - p.a, "{}: m = n - a", p.name);
        assert_eq!(
            p.q.pow(p.d as u32) + 1,
            p.dd as u32,
            "{}: D = q^d + 1",
            p.name
        );
        assert_eq!(p.fpoly.len(), p.n, "{}: field polynomial degree", p.name);
        assert!(p.fpoly[0] != 0, "{}: t divides f", p.name);
        assert!(
            p.fpoly.iter().all(|&c| (c as u32) < p.q),
            "{}: coefficient range",
            p.name
        );
        assert_eq!(p.is_dragon(), p.ny > 0, "{}: Dragon iff ny", p.name);
        for (index, &(i, j)) in p.monomials.iter().enumerate() {
            assert!(
                !p.monomials[..index].contains(&(i, j)),
                "{}: duplicate monomial",
                p.name
            );
            assert!(i <= j && j <= p.d, "{}: monomial order", p.name);
            assert!(
                p.q.pow(i as u32) + p.q.pow(j as u32) <= p.dd as u32,
                "{}: monomial degree",
                p.name
            );
        }
        // In characteristic 2 (q=2 and q=4 here), a diagonal q-polynomial
        // monomial is additive and F_2-linear, so these sets use no i = j term.
        if p.q == 2 || p.q == 4 {
            assert!(p.monomials.iter().all(|&(i, j)| i != j), "{}", p.name);
        }
    }
}

#[test]
fn tags_are_well_formed() {
    assert_eq!(
        by_name("d-james-128-q2").unwrap().tag(),
        "D-James/v1/d-james/q2/n189/a27/r2/D5/ny256/mon0-1.0-2"
    );
    assert_eq!(
        by_name("james-128-q2").unwrap().tag(),
        "D-James/v1/james/q2/n283/a27/r2/D5/ny0/mon0-1.0-2"
    );
    assert_eq!(
        by_name("d-james-128-q23").unwrap().tag(),
        "D-James/v1/d-james/q23/n74/a21/r2/D24/ny57/mon0-0.0-1"
    );
    assert_eq!(
        by_name("d-james-128-q4").unwrap().tag(),
        "D-James/v1/d-james/q4/n105/a21/r2/D17/ny128/mon0-1.0-2"
    );
    let q4_256 = by_name("d-james-256-q4").unwrap();
    assert_eq!((q4_256.m, q4_256.a, q4_256.n), (170, 53, 223));
    // Distinct sets must never share a tag: keys are bound to it.
    let mut tags: Vec<String> = ALL.iter().map(|p| p.tag()).collect();
    tags.sort();
    let n = tags.len();
    tags.dedup();
    assert_eq!(tags.len(), n, "duplicate parameter tag");
}

#[test]
fn signature_sizes_match_the_spec() {
    // (name, bits, bytes) from d-james-spec.md section 9.2
    let want: &[(&str, usize)] = &[
        ("d-james-128-q2", 24),
        ("d-james-128-q4", 27),
        ("d-james-128-q5", 28),
        ("d-james-128-q13", 36),
        ("d-james-128-q23", 42),
        ("d-james-256-q2", 49),
        ("d-james-256-q4", 56),
        ("d-james-256-q5", 60),
        ("d-james-256-q13", 80),
        ("d-james-256-q23", 93),
        ("james-128-q2", 36),
        ("james-128-q4", 38),
        ("james-128-q5", 39),
        ("james-128-q13", 43),
        ("james-128-q23", 45),
        ("james-256-q2", 73),
        ("james-256-q4", 78),
        ("james-256-q5", 80),
        ("james-256-q13", 89),
        ("james-256-q23", 95),
    ];
    for &(name, bytes) in want {
        assert_eq!(by_name(name).unwrap().sig_bytes(), bytes, "{name}");
    }
}
