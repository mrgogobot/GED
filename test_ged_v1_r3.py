"""Independent mathematical and validation checks for GED-V1-R3.

WHAT: Guard scientific invariants and catch unintended code changes.
WHY: Passing tests validates IMPLEMENTATION, never an unobserved new field.
"""
import copy
from decimal import Decimal, localcontext
import json
from pathlib import Path
import unittest

import ged_v1_r2 as r2
import ged_v1_r3 as ged


CONFIG = json.loads((Path(__file__).parent / "GED_V1_R3_config.json").read_text())


def d(value):
    """WHAT: Obtain exact decimal value for high-precision comparisons."""
    return Decimal(str(value))


class GEDV1R3Tests(unittest.TestCase):
    """WHY: Separate physical conservation checks from executable formatting."""

    def setUp(self):
        # WHAT: Recalculate using each test's own fresh config dictionary.
        # WHY: One test should never alter another test's inputs.
        self.cfg = copy.deepcopy(CONFIG)
        self.core = ged.binding_transition(self.cfg)

    def case(self, fraction, *, stable=False):
        return ged.scenario(self.core, fraction, stable=stable)

    def test_pinned_r2_source(self):
        self.assertEqual(
            ged.sha256(Path(r2.__file__)), CONFIG["r2_dependency_sha256"]
        )

    def test_r2_inventory_matches(self):
        A = r2.make_body("100", "0.5", "7870", "917")
        B = r2.make_body("100", "0.5", "917", "7870")
        self.assertEqual(A["M0"], B["M0"])

    def test_independent_volume_mass(self):
        # WHAT: Cross-check the R2 constituent mass using a direct volume law.
        with localcontext() as ctx:
            ctx.prec = 60
            direct = (d(4) * r2.PI * d(100)**3 / d(3)) * (d(7870)+d(917))/d(2)
            self.assertLess(abs(direct-self.core["constituent_mass_kg"]), d("1e-40"))

    def test_binding_release_positive(self):
        self.assertGreater(self.core["release_J"], 0)

    def test_binding_release_matches_user_r2(self):
        # WHAT: Cross-check the previous locally reported V1-R2 number.
        self.assertLess(abs(self.core["release_J"]-d("39713221.06352")), d("0.00001"))

    def test_binding_potentials_negative(self):
        self.assertLess(self.core["U_A_J"], 0)
        self.assertLess(self.core["U_B_J"], 0)

    def test_iron_core_more_gravitationally_bound(self):
        self.assertLess(self.core["U_A_J"], self.core["U_B_J"])

    def test_capacity_independently_weighted(self):
        # WHAT: Sum m*c_p for iron and ice, independently of any reordering.
        with localcontext() as ctx:
            ctx.prec = 60
            v_half = d(2)*r2.PI*d(100)**3/d(3)
            expected = v_half*(d(7870)*d(450)+d(917)*d(2100))
            self.assertLess(abs(self.core["heat_capacity_J_K"]-expected), d("1e-35"))

    def test_stable_requires_no_exhaust(self):
        x = self.case("0", stable=True)
        self.assertEqual(d(x["energy_to_bath_J"]), 0)
        self.assertEqual(d(x["potential_energy_change_J"]), 0)
        self.assertEqual(d(x["thermal_energy_retained_J"]), 0)

    def test_stable_has_zero_modeled_entropy_change(self):
        self.assertEqual(d(self.case("0", stable=True)["modeled_thermal_entropy_change_J_K"]), 0)

    def test_stable_mass_does_not_change(self):
        self.assertEqual(d(self.case("0", stable=True)["body_mass_change_kg"]), 0)

    def test_closed_keeps_all_release_as_heat(self):
        x = self.case("0")
        self.assertEqual(d(x["thermal_energy_retained_J"]), self.core["release_J"])

    def test_closed_no_export(self):
        self.assertEqual(d(self.case("0")["energy_to_bath_J"]), 0)

    def test_closed_entropy_increases_without_export(self):
        x = self.case("0")
        self.assertGreater(d(x["body_entropy_change_J_K"]), 0)
        self.assertEqual(d(x["bath_entropy_change_J_K"]), 0)

    def test_closed_mass_does_not_change(self):
        self.assertEqual(d(self.case("0")["body_mass_change_kg"]), 0)

    def test_closed_exterior_gravity_does_not_change(self):
        self.assertEqual(d(self.case("0")["exterior_delta_g_m_s2"]), 0)

    def test_full_export_leaves_no_heat(self):
        self.assertEqual(d(self.case("1")["thermal_energy_retained_J"]), 0)

    def test_full_export_does_not_warm_body(self):
        self.assertEqual(d(self.case("1")["final_body_temperature_K"]), d("250"))

    def test_full_export_environment_gains_entropy(self):
        self.assertGreater(d(self.case("1")["bath_entropy_change_J_K"]), 0)

    def test_full_export_body_mass_decreases(self):
        self.assertLess(d(self.case("1")["body_mass_change_kg"]), 0)

    def test_full_export_delta_g_negative(self):
        self.assertLess(d(self.case("1")["exterior_delta_g_m_s2"]), 0)

    def test_full_export_matches_r2_mass_equivalence(self):
        with localcontext() as ctx:
            ctx.prec = 60
            want = -(self.core["release_J"]/d("299792458")**2)
            got = d(self.case("1")["body_mass_change_kg"])
            self.assertLess(abs(want-got), d("1e-65"))

    def test_half_export_split(self):
        x = self.case("0.5")
        self.assertEqual(d(x["thermal_energy_retained_J"]), d(x["energy_to_bath_J"]))

    def test_half_export_mass_half_of_full(self):
        # WHAT: Keep the test's precision equal to the calculation's.
        # WHY: Default Decimal precision (28) rounds tiny, long mass values.
        with localcontext() as ctx:
            ctx.prec = 60
            self.assertEqual(d(self.case("0.5")["body_mass_change_kg"])*2,
                             d(self.case("1")["body_mass_change_kg"]))

    def test_energy_balance_all_fractions(self):
        for f in ("0", "0.25", "0.5", "0.75", "1"):
            with self.subTest(f=f):
                self.assertEqual(d(self.case(f)["energy_balance_residual_J"]), 0)

    def test_global_mass_energy_balance_all_fractions(self):
        for f in ("0", "0.25", "0.5", "0.75", "1"):
            with self.subTest(f=f):
                x = self.case(f)
                self.assertEqual(d(x["global_mass_energy_residual_kg"]), 0)

    def test_entropy_nonnegative_all_fractions(self):
        for f in ("0", "0.25", "0.5", "0.75", "1"):
            with self.subTest(f=f):
                self.assertGreaterEqual(d(self.case(f)["modeled_thermal_entropy_change_J_K"]), 0)

    def test_body_temperature_monotone_with_retention(self):
        self.assertGreater(d(self.case("0")["final_body_temperature_K"]),
                           d(self.case("0.5")["final_body_temperature_K"]))
        self.assertGreater(d(self.case("0.5")["final_body_temperature_K"]),
                           d(self.case("1")["final_body_temperature_K"]))

    def test_thermal_entropy_uses_logarithm(self):
        x = self.case("0.5")
        with localcontext() as ctx:
            ctx.prec = 60
            independent = self.core["heat_capacity_J_K"] * (
                d(x["final_body_temperature_K"]) / d("250")
            ).ln()
            self.assertLess(abs(independent-d(x["body_entropy_change_J_K"])), d("1e-48"))

    def test_external_entropy_formula(self):
        x = self.case("0.5")
        with localcontext() as ctx:
            ctx.prec = 60
            self.assertEqual(d(x["bath_entropy_change_J_K"]),
                             d(x["energy_to_bath_J"]) / d("200"))

    def test_bath_temperature_changes_entropy_not_mass(self):
        core2 = dict(self.core, environment_temperature_K=d("150"))
        a = self.case("0.5")
        b = ged.scenario(core2, "0.5")
        self.assertEqual(a["body_mass_change_kg"], b["body_mass_change_kg"])
        self.assertGreater(d(b["bath_entropy_change_J_K"]), d(a["bath_entropy_change_J_K"]))

    def test_gravity_delta_inverse_square_distance(self):
        with localcontext() as ctx:
            ctx.prec = 60
            farther = dict(self.core, probe_radius_m=self.core["probe_radius_m"]*2)
            a = d(self.case("1")["exterior_delta_g_m_s2"])
            b = d(ged.scenario(farther, "1")["exterior_delta_g_m_s2"])
            # WHAT: Verify scaling to relative 1e-55 precision.
            # WHY: Decimal division rounds the final digit even at 60 digits.
            self.assertLess(abs((a-4*b)/a), d("1e-55"))

    def test_delta_g_matches_independent_formula(self):
        x = self.case("1")
        with localcontext() as ctx:
            ctx.prec = 60
            independent = (d("6.67430e-11") * d(x["body_mass_change_kg"])
                           / d("200")**2)
            self.assertLess(abs(independent-d(x["exterior_delta_g_m_s2"])), d("1e-76"))

    def test_fractional_gravity_matches_mass_fraction(self):
        x = self.case("0.5")
        with localcontext() as ctx:
            ctx.prec = 60
            expected = d(x["body_mass_change_kg"]) / self.core["constituent_mass_kg"]
            self.assertLess(abs(expected-d(x["exterior_fractional_delta_g"])), d("1e-75"))

    def test_reject_fraction_below_zero(self):
        with self.assertRaises(ValueError):
            self.case("-0.01")

    def test_reject_fraction_over_one(self):
        with self.assertRaises(ValueError):
            self.case("1.01")

    def test_reject_nonfinite_fraction(self):
        with self.assertRaises(ValueError):
            self.case("NaN")

    def test_reject_nonpositive_temperature(self):
        self.cfg["initial_temperature_K"] = "0"
        with self.assertRaises(ValueError):
            ged.binding_transition(self.cfg)

    def test_reject_nonpositive_environment_temperature(self):
        self.cfg["environment_temperature_K"] = "-10"
        with self.assertRaises(ValueError):
            ged.binding_transition(self.cfg)

    def test_reject_negative_specific_heat(self):
        self.cfg["heat_capacity_iron_J_kg_K"] = "-1"
        with self.assertRaises(ValueError):
            ged.binding_transition(self.cfg)

    def test_reject_unequal_material_inventory(self):
        self.cfg["core_volume_fraction"] = "0.7"
        with self.assertRaises(ValueError):
            ged.binding_transition(self.cfg)

    def test_reject_probe_inside_asteroid(self):
        self.cfg["report_probe_R"] = "0.8"
        with self.assertRaises(ValueError):
            ged.binding_transition(self.cfg)

    def test_reject_duplicate_sweep(self):
        self.cfg["escaped_fractions"] = ["0", "0.5", "0.5", "1"]
        with self.assertRaises(ValueError):
            ged.calculate(self.cfg)

    def test_reject_unsorted_sweep(self):
        self.cfg["escaped_fractions"] = ["0", "1", "0.5"]
        with self.assertRaises(ValueError):
            ged.calculate(self.cfg)

    def test_reject_missing_endpoint_sweep(self):
        self.cfg["escaped_fractions"] = ["0.25", "0.5", "1"]
        with self.assertRaises(ValueError):
            ged.calculate(self.cfg)

    def test_results_are_unobserved(self):
        result, rows = ged.calculate(self.cfg)
        self.assertFalse(result["observation_performed"])
        self.assertEqual(result["classification"], "ESTABLISHED_PHYSICS_THERMAL_LEDGER_SCOUT")
        self.assertEqual(len(rows), 6)


if __name__ == "__main__":
    unittest.main()
