"""GED-V0-R1 independent numerical/logic checks (stdlib unittest).

Run from this folder with:
    python -m unittest -v test_ged_v0_r1.py

Tests are not themselves scientific confirmation of an information-field theory.
They check that the implemented baseline matches its analytic assumptions.
"""

import math
from pathlib import Path
import unittest

import ged_v0_r1 as ged


# WHAT: Group checks so one command verifies all core assumptions.
# WHY: Incorrect code must not be mistaken for new physics.
class GEDBaselineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # WHY: Use the same declared inputs as the science script, not a hidden
        # second set of constants that could silently disagree.
        cls.cfg = ged.load_config(Path(__file__).with_name("GED_V0_R1_config.json"))
        cls.G = cls.cfg["constants"]["G_m3_kg-1_s-2"]
        cls.iron = ged.read_material(cls.cfg, "iron")
        cls.ice = ged.read_material(cls.cfg, "ice")
        cls.R = cls.cfg["sphere_radius_m"]
        cls.result = ged.calculate(cls.cfg)

    def test_sphere_mass_uses_density_times_volume(self):
        # WHY: Checks the fundamental dimensional relationship M=rho*V.
        sphere = ged.Sphere(self.R, self.ice)
        expected = self.ice.density_kg_m3 * (4/3) * math.pi * self.R**3
        self.assertAlmostEqual(sphere.mass_kg / expected, 1.0, places=13)

    def test_surface_g_equals_newton_formula(self):
        # WHY: Checks a measured gravitational quantity in SI units.
        M = ged.Sphere(self.R, self.iron).mass_kg
        self.assertAlmostEqual(ged.surface_gravity_m_s2(self.G, M, self.R), self.G*M/self.R**2)

    def test_density_ratio_equals_surface_gravity_ratio(self):
        # WHY: If radii match, g scales directly with average density.
        equal = self.result["equal_radius_comparison"]
        observed = equal["iron"]["surface_gravity_m_s2"] / equal["ice"]["surface_gravity_m_s2"]
        expected = self.iron.density_kg_m3 / self.ice.density_kg_m3
        self.assertAlmostEqual(observed/expected, 1.0, places=13)

    def test_layer_swapping_same_external_gravity(self):
        # WHY: Core/shell ordering cannot be deduced from an idealized spherical
        # body's external Newtonian field if total mass and radius are fixed.
        layers = self.result["layered_same_mass_radius"]
        a = layers["iron_core_ice_shell"]["exterior_g_m_s2"]
        b = layers["ice_core_iron_shell"]["exterior_g_m_s2"]
        self.assertTrue(math.isclose(a, b, rel_tol=1e-12))
        self.assertTrue(layers["external_fields_equal_by_spherical_symmetry"])

    def test_layer_exterior_matches_independent_point_mass_formula(self):
        # WHY: Check the core-shell output against the analytic total-mass
        # prediction, not merely against the second core-shell calculation.
        layers = self.result["layered_same_mass_radius"]
        r = layers["exterior_probe_radius_m"]
        expected = self.G * layers["total_mass_each_kg"] / (r*r)
        actual = layers["iron_core_ice_shell"]["exterior_g_m_s2"]
        self.assertTrue(math.isclose(actual, expected, rel_tol=2e-15))

    def test_layer_swapping_changes_interior_gravity(self):
        # WHY: Internal fields DO carry information about enclosed density.
        layers = self.result["layered_same_mass_radius"]
        a = layers["iron_core_ice_shell"]["interior_g_m_s2"]
        b = layers["ice_core_iron_shell"]["interior_g_m_s2"]
        self.assertGreater(a, b)

    def test_thermal_melting_has_known_limits(self):
        # WHY: The same energy input gives different outcomes because materials
        # have different properties, without adding a novel gravitational rule.
        heat = self.result["thermal_response"]
        self.assertEqual(heat["melt_fraction_iron"], 0.0)
        expected = (300000 - (273.15-250)*2100)/333000
        self.assertAlmostEqual(heat["melt_fraction_ice"], expected, places=13)
        self.assertTrue(0 <= heat["melt_fraction_ice"] <= 1)

    def test_mass_energy_equivalent_uses_c_squared(self):
        # WHY: Added retained energy is not exactly gravitationally invisible.
        heat = self.result["thermal_response"]
        c = self.cfg["constants"]["c_m_s-1"]
        self.assertAlmostEqual(heat["equivalent_added_mass_per_kg_kg"], 300000/c**2)

    def test_fragment_geometry_and_mass_conservation(self):
        # WHY: An imposed fragmentation must not create or remove matter.
        frag = self.result["fragmentation_ice"]
        self.assertEqual(frag["count"], 1000)
        self.assertAlmostEqual(frag["radius_per_fragment_m"], 10.0)
        self.assertTrue(frag["mass_conserved"])
        self.assertAlmostEqual(frag["sum_fragment_mass_kg"]/frag["original_ice_mass_kg"], 1.0)

    def test_fragment_gravity_scales_with_radius(self):
        # WHY: With density unchanged, a fragment 1/10 as large in radius
        # has 1/10 the original sphere's individual surface acceleration.
        frag = self.result["fragmentation_ice"]
        original_g = self.result["equal_radius_comparison"]["ice"]["surface_gravity_m_s2"]
        self.assertAlmostEqual(frag["fragment_surface_g_m_s2"]/original_g, 0.1)

    def test_far_field_approaches_total_mass_at_centre(self):
        # WHY: At a sufficiently distant probe, a compact fragment group
        # approximately behaves like a point mass of the combined total mass.
        probes = self.result["fragmentation_ice"]["probes"]
        near_change = probes["near"]["fractional_difference_magnitude"]
        far_change = probes["far"]["fractional_difference_magnitude"]
        self.assertLess(far_change, 1e-5)
        self.assertLess(far_change, near_change)

    def test_invalid_inputs_rejected(self):
        # WHY: Silent acceptance of negative mass or radius poisons experiments.
        with self.assertRaises(ValueError):
            ged.surface_gravity_m_s2(self.G, 1.0, -1.0)
        with self.assertRaises(ValueError):
            ged.melt_fraction(self.ice, -1.0, 1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
