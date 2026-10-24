"""GED-V1-R1 independent checks; results are only numerical verifications."""
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest

from ged_v1_r1 import (
    LayeredSphere, alpha_sensitivity, calculate, enclosed_mass, gravity,
    load_config, material, melt_fraction, melt_threshold, point_mass_gravity,
    shell_contrast, toy_gravity, volume,
)


class GEDV1R1Tests(unittest.TestCase):
    # WHAT: Read the real configuration for all tests.
    # WHY: The test suite must exercise the actual documented input values.
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(Path(__file__).with_name("GED_V1_R1_config.json"))
        c = cls.cfg
        cls.G = c["G_m3_kg_s2"]
        cls.R = c["radius_m"]
        cls.iron, cls.ice = material("iron", c), material("ice", c)
        cls.a = LayeredSphere(cls.R, 0.5, cls.iron, cls.ice)
        cls.b = LayeredSphere(cls.R, 0.5, cls.ice, cls.iron)

    def test_volume_units(self):
        # 1-metre radius gives 4*pi/3 cubic metres.
        self.assertAlmostEqual(volume(1), 4 * math.pi / 3)

    def test_mass_independent_sum(self):
        # Independently check mass from the two half-volumes.
        want = volume(self.R) * (self.iron.density + self.ice.density) / 2
        self.assertAlmostEqual(self.a.mass / want, 1.0, places=13)

    def test_matched_masses(self):
        self.assertAlmostEqual(self.a.mass / self.b.mass, 1.0, places=13)

    def test_half_core_volume(self):
        self.assertAlmostEqual(volume(self.a.core_radius) / volume(self.R), 0.5, places=13)

    def test_enclosed_total_mass(self):
        self.assertAlmostEqual(enclosed_mass(self.a, self.R) / self.a.mass, 1, places=13)

    def test_central_gravity_zero(self):
        self.assertEqual(gravity(self.a, 0, self.G), 0)

    def test_interior_gravity_differs(self):
        self.assertGreater(gravity(self.a, self.R / 2, self.G),
                           gravity(self.b, self.R / 2, self.G))

    def test_exterior_null(self):
        # Swapping the materials cannot change the ideal spherical exterior.
        for factor in (1, 2, 10):
            self.assertAlmostEqual(
                gravity(self.a, factor*self.R, self.G) / gravity(self.b, factor*self.R, self.G),
                1.0, places=13,
            )

    def test_independent_point_mass(self):
        r = 7*self.R
        self.assertAlmostEqual(
            gravity(self.a, r, self.G), point_mass_gravity(self.a.mass, r, self.G),
            places=17,
        )

    def test_heat_thresholds_differ(self):
        T = self.cfg["temperature_K"]
        self.assertLess(melt_threshold(self.ice, T), melt_threshold(self.iron, T))

    def test_thermal_response_differs(self):
        T, heat = self.cfg["temperature_K"], self.cfg["heat_input_J_kg"]
        self.assertEqual(melt_fraction(self.iron, T, heat), 0)
        self.assertGreater(melt_fraction(self.ice, T, heat), 0.75)

    def test_normalized_contrast(self):
        T = self.cfg["temperature_K"]
        self.assertAlmostEqual(shell_contrast(self.ice, self.iron, self.ice, T), -1)
        self.assertAlmostEqual(shell_contrast(self.iron, self.iron, self.ice, T), +1)

    def test_zero_toy_parameter_returns_null(self):
        r = 2*self.R
        for body, chi in ((self.a, -1), (self.b, +1)):
            self.assertEqual(toy_gravity(body, r, self.G, 0, chi, 2), gravity(body, r, self.G))

    def test_nonzero_toy_prediction(self):
        # Independent algebra: delta g / g0 = 2*alpha*(R/r)^p.
        r, alpha = 2*self.R, 1e-6
        diff = toy_gravity(self.b, r, self.G, alpha, 1, 2) - toy_gravity(self.a, r, self.G, alpha, -1, 2)
        self.assertAlmostEqual(diff / point_mass_gravity(self.a.mass, r, self.G),
                               2*alpha / 4, delta=1e-14)

    def test_radial_falloff(self):
        alpha = 1e-6
        def relative_difference(r):
            high = toy_gravity(self.b, r, self.G, alpha, 1, 2)
            low = toy_gravity(self.a, r, self.G, alpha, -1, 2)
            return (high - low) / gravity(self.a, r, self.G)
        # The ratio uses differences of nearly equal floats (~1e-7 relative).
        # A 1e-8 relative tolerance accounts for subtraction round-off.
        self.assertTrue(math.isclose(
            relative_difference(2*self.R) / relative_difference(4*self.R),
            4.0, rel_tol=1e-8
        ))

    def test_one_ppm_threshold(self):
        self.assertAlmostEqual(alpha_sensitivity(1.0, 2.0, 2), 2e-6)

    def test_output_is_finite_and_unobserved(self):
        result, rows = calculate(self.cfg)
        self.assertFalse(result["observation_performed"])
        self.assertTrue(result["baseline"]["exterior_fields_equal"])
        self.assertTrue(all(math.isfinite(n) for row in rows for n in row.values()))

    def test_reject_inside_toy(self):
        with self.assertRaises(ValueError):
            toy_gravity(self.a, self.R/2, self.G, 1e-6, 1, 2)

    def test_reject_nonfinite_heat(self):
        with self.assertRaises(ValueError):
            melt_fraction(self.ice, 250, float("nan"))

    def test_reject_wrong_volume_fraction(self):
        altered = copy.deepcopy(self.cfg)
        altered["core_volume_fraction"] = 0.45
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            path.write_text(json.dumps(altered), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_config(path)


if __name__ == "__main__":
    unittest.main()
