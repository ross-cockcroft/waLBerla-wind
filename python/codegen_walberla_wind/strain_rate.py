import pystencils as ps
import sympy as sp

from codegen_walberla_wind.guards import require_found


def add_strain_rate_output(collision_rule, strain_rate_field, omega_field, density_field):
    n = collision_rule.method.stencil.D ** 2
    qr_stress_symbols = [sp.Symbol(f"sgs_qr_pi{i}") for i in range(n)]
    defined_symbols = {a.lhs for a in collision_rule.subexpressions}
    require_found(all(symbol in defined_symbols for symbol in qr_stress_symbols),
                  f"the sgs_qr_pi0 to sgs_qr_pi{n - 1} subexpressions", "lbmpy.turbulence_models")
    effective_omega = require_found(next((a.rhs for a in collision_rule.main_assignments if a.lhs == omega_field.center), None),
                                    "the omega output assignment", "lbmpy.creationfunctions")
    density_expression = require_found(next((a.rhs for a in collision_rule.main_assignments if a.lhs == density_field.center), None),
                                       "the density output assignment", "lbmpy.creationfunctions")
    collision_rule.main_assignments.extend(
        ps.Assignment(strain_rate_field.center_vector[i], -3 * effective_omega / (2 * density_expression) * qr_stress_symbols[i])
        for i in range(n))
    return collision_rule
