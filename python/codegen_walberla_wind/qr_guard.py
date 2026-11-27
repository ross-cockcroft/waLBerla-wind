import sympy as sp
from pystencils import Assignment

from codegen_walberla_wind.guards import require_found


def guard_qr_division(collision_rule):
    c_pi, r_prime = sp.symbols("qr_c_pi qr_r_prime")
    require_found(any(a.lhs == c_pi for a in collision_rule.subexpressions), "the qr_c_pi assignment", "lbmpy.turbulence_models")
    guarded = [Assignment(a.lhs, sp.Piecewise((a.rhs, r_prime > 0), (sp.Float(0), True))) if a.lhs == c_pi else a
               for a in collision_rule.subexpressions]
    return collision_rule.copy(subexpressions=guarded)
