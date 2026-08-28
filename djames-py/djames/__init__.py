"""D-James / James: ultra-short multivariate signatures.

A clean-room, stdlib-only implementation of the schemes in
Jacques Patarin and Alexandre Roullet, "D-James: Ultra Short Multivariate
Signatures", IACR ePrint 2026/1650, https://eprint.iacr.org/2026/1650

This is a reference implementation for study and test-vector generation.  It
is not constant-time and has had no security review.  The paper's own authors
recommend against deploying these schemes until they have had substantially
more public scrutiny; nothing here changes that.
"""

from .params import get as get_params, names as param_names, ALL as PARAMS
from .scheme import keygen, sign, verify, PublicKey, SecretKey

__all__ = ["get_params", "param_names", "PARAMS",
           "keygen", "sign", "verify", "PublicKey", "SecretKey"]
__version__ = "1.0.0"
