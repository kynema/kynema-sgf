"""KLAxellSeparation variants of the Ishihara-Zhou ridge and hill cases.

Usage: python3 separation_cases.py [output directory, default ./cases]

Uses make_cases.py, which writes the KLAxell cases with the TabulatedProfile
inflow, to write the same 2D ridge and 3D hill cases for two KLAxellSeparation
setups:
  limiter  pressure-gradient sensor, realizable Cmu limiter, 10 s relaxed gate
  all      limiter plus production cap, destruction boost, curvature correction
           (richardson) and implicit dissipation
The flat approach-flow cases do not depend on the separation treatments and are
not repeated. See SEPARATION.md.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_cases import HERE, generate  # noqa: E402

LIMITER = {
    "turbulence.model": "KLAxellSeparation",
    "KLAxellSeparation.pressure_gradient_sensor": "true",
    "KLAxellSeparation.realizable_cmu": "true",
    "KLAxellSeparation_coeffs.gate_relaxation_time": "10.0",
    "io.outputs": "pressure_gradient_sensor separation_gate",
}
SEPARATION_MODELS = {
    "limiter": LIMITER,
    "all": {
        **LIMITER,
        "KLAxellSeparation.production_cap": "true",
        "KLAxellSeparation.destruction_boost": "true",
        "KLAxellSeparation.curvature_correction": "true",
        "KLAxellSeparation.implicit_dissipation": "true",
    },
}

if __name__ == "__main__":
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "cases"
    generate(out, SEPARATION_MODELS, approach=False)
