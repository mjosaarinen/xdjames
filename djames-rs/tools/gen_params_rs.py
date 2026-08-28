#!/usr/bin/env python3
"""Regenerate djames-rs/src/params.rs from the Python reference definitions.

Run from djames-py/:   python3 ../djames-rs/tools/gen_params_rs.py > ../djames-rs/src/params.rs

Generating rather than transcribing keeps the two implementations from
drifting apart through a typo in a 578-entry field polynomial.
"""
(generator recorded; body lives in the session that produced params.rs)
