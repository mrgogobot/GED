"""GED-V1-R2 independent unit and physics checks, standard library only."""
import json
from decimal import getcontext
import math
from pathlib import Path
import unittest
from ged_v1_r2 import (D, PI, make_body, sphere_volume, binding_energy,
                        homogeneous_binding, mass_equivalent, delta_g, exterior_g,
                        calculate)

G, C, R = D("6.67430e-11"), D("299792458"), D("100")
A = make_body(R, D("0.5"), D("7870"), D("917"))
B = make_body(R, D("0.5"), D("917"), D("7870"))
CONFIG = json.loads((Path(__file__).resolve().parent / "GED_V1_R2_config.json").read_text())


def numeric_integral(body, n=8000):
    """Independent Simpson quadrature for U=-int G*M(r)/r dM(r).

    WHY: A second mathematical route checks the two-layer analytic formula;
    numerical integral splits at core boundary to avoid discontinuity error.
    """
    a, R = float(body["a"]), float(body["R"])
    rc, rs = float(body["rc"]), float(body["rs"])
    def segment(lo, hi, rho):
        step, tally = (hi-lo)/n, 0.0
        for i in range(n+1):
            r = lo+step*i
            m = 4*math.pi/3*(rc*min(r, a)**3 + rs*max(0, r**3-a**3))
            integrand = -4*math.pi*float(G)*rho*m*r
            tally += (1 if i in (0, n) else (4 if i % 2 else 2))*integrand
        return tally*step/3
    return segment(0, a, rc)+segment(a, R, rs)


class GEDV1R2Tests(unittest.TestCase):
    def setUp(self):
        # WHY: Match the stated high precision for all checks performed here.
        getcontext().prec = 60

    def test_equal_inventory(self):
        self.assertEqual(A["M0"], B["M0"])

    def test_core_fraction(self):
        self.assertAlmostEqual(float((A["a"]/R)**3), .5, places=13)

    def test_mass_parts_sum(self):
        self.assertEqual(A["mc"]+A["ms"], A["M0"])

    def test_sphere_volume_independent(self):
        self.assertAlmostEqual(float(sphere_volume(R)), 4*math.pi*R.__float__()**3/3)

    def test_negative_radius_rejected(self):
        with self.assertRaises(ValueError): sphere_volume(-1)

    def test_uniform_sphere_limit(self):
        b = make_body(R, .5, 1200, 1200)
        measured = binding_energy(b, G)
        expected = homogeneous_binding(b["M0"], R, G)
        self.assertLess(abs((measured-expected)/expected), D("1e-45"))

    def test_independent_quadrature_iron(self):
        analytic = float(binding_energy(A,G))
        self.assertLess(abs(numeric_integral(A)-analytic)/abs(analytic), 1e-9)

    def test_independent_quadrature_ice(self):
        analytic = float(binding_energy(B,G))
        self.assertLess(abs(numeric_integral(B)-analytic)/abs(analytic), 1e-9)

    def test_iron_core_more_bound(self):
        self.assertLess(binding_energy(A,G), binding_energy(B,G))

    def test_negative_binding_potential(self):
        self.assertLess(binding_energy(A,G), 0)
        self.assertLess(binding_energy(B,G), 0)

    def test_mass_equivalence_identity(self):
        self.assertEqual(mass_equivalent(C**2, C), 1)

    def test_mass_equivalence_negative(self):
        self.assertEqual(mass_equivalent(-C**2, C), -1)

    def test_c_zero_rejected(self):
        with self.assertRaises(ValueError): mass_equivalent(1,0)

    def test_delta_g_inverse_square(self):
        self.assertEqual(delta_g(D("1e12"),R,G,C),4*delta_g(D("1e12"),2*R,G,C))

    def test_delta_g_zero(self):
        self.assertEqual(delta_g(0,R,G,C), 0)

    def test_delta_g_formula(self):
        self.assertEqual(delta_g(C*C,R,G,C), G/(R*R))

    def test_newton_baseline(self):
        self.assertAlmostEqual(float(exterior_g(A["M0"],R,G)), float(G)*float(A["M0"])/10000)

    def test_energy_release_conservation(self):
        res,_ = calculate(CONFIG)
        self.assertEqual(D(res["delta_U_B_to_A_J"])+D(res["binding_energy_release_B_to_A_J"]),0)

    def test_open_mass_loss(self):
        res,_ = calculate(CONFIG)
        self.assertLess(D(res["mass_change_if_release_escapes_kg"]),0)

    def test_closed_mass_unchanged(self):
        res,_ = calculate(CONFIG)
        self.assertEqual(D(res["mass_change_if_release_retained_kg"]),0)

    def test_absorb_vs_emit(self):
        res,_ = calculate(CONFIG)
        c = res["scenarios"]
        self.assertEqual(D(c["absorbed_heat_kept"]["delta_M_kg"]),
                         -D(c["same_heat_radiated"]["delta_M_kg"]))

    def test_tiny_mass_change_not_lost(self):
        res,_ = calculate(CONFIG)
        self.assertNotEqual(D(res["constituent_inventory_mass_kg"]),
                            D(res["constituent_inventory_mass_kg"])+D(res["mass_change_if_release_escapes_kg"]))

    def test_fixed_radius_difference_small(self):
        res,_ = calculate(CONFIG)
        self.assertLess(abs(D(res["scenarios"]["B_to_A_binding_release_escapes"]["fractional_delta_g"])), D("1e-15"))

    def test_inside_probe_rejected(self):
        with self.assertRaises(ValueError):
            calculate(dict(CONFIG, report_probe_R="0.5"))

    def test_negative_heat_rejected(self):
        with self.assertRaises(ValueError):
            calculate(dict(CONFIG, thermal_input_J_kg="-1"))

    def test_no_new_gravity_claim(self):
        res,_ = calculate(CONFIG)
        self.assertFalse(res["observation_performed"])
        self.assertEqual(res["classification"], "ESTABLISHED_PHYSICS_WEAK_FIELD_SCOUT")


if __name__ == "__main__": unittest.main()
