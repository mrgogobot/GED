"""GED-V1-R8 independent verification tests.

WHAT: Unit-test geometry, analytic Newtonian fields, inverse classification,
      exactly non-unique spherical shells, data integrity and scientific scope.
WHY:  Green tests check the implementation, not the novel hypothesis.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import random
import unittest

import ged_v1_r8 as ged


class GEDR8Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # WHAT: Use the shipped configuration, as the local run will.
        # WHY:  Our independent reference formulas use its actual parameters.
        cls.config = json.loads((Path(__file__).resolve().parent / ged.CONFIG_NAME).read_text())
        cls.result = ged.compute(cls.config)
        cls.library = ged.signatures(cls.config)
        cls.G = cls.config["gravitational_constant_m3_kg_s2"]
        cls.M = cls.config["added_energy_j"] / cls.config["speed_of_light_m_s"]**2
        cls.d = cls.config["parcel_offset_radius_m"]

    def with_change(self, key, value):
        altered = copy.deepcopy(self.config)
        altered[key] = value
        return altered

    def test_experiment_id(self):
        self.assertEqual(self.result["experiment_id"], "GED-V1-R8")

    def test_validation_passes_for_frozen_config(self):
        self.assertIsNone(ged.validate(self.config))

    def test_synthetic_disclosure(self):
        self.assertFalse(self.result["observed"])

    def test_no_new_gravity_claim(self):
        self.assertFalse(self.config["claims"]["new_gravity_detected"])
        self.assertTrue(any("No new" in limit for limit in self.result["limits"]))

    def test_known_template_scope_recorded(self):
        self.assertIn("restricted known-template", self.result["recoverability_scope"])

    def test_probes_outside_body(self):
        self.assertTrue(all(ged.vec_length(p) > self.config["nominal_body_radius_m"]
                            for p in ged.probe_locations(self.config)))

    def test_probe_count(self):
        self.assertEqual(len(ged.probe_locations(self.config)), 24)

    def test_three_components_per_probe(self):
        self.assertTrue(all(len(p) == 3 for p in ged.probe_locations(self.config)))

    def test_probe_first_coordinate(self):
        self.assertEqual(ged.probe_locations(self.config)[0], (200., 0., 0.))

    def test_probe_last_coordinate(self):
        self.assertEqual(ged.probe_locations(self.config)[-1], (0., 0., -800.))

    def test_pair_x_has_two_positions(self):
        self.assertEqual(ged.parcel_positions("pair_x", self.d), ((self.d, 0., 0.), (-self.d, 0., 0.)))

    def test_pair_y_has_two_positions(self):
        self.assertEqual(len(ged.parcel_positions("pair_y", self.d)), 2)

    def test_cube_has_eight_positions(self):
        self.assertEqual(len(ged.parcel_positions("cube", self.d)), 8)

    def test_cube_sources_same_radius(self):
        for location in ged.parcel_positions("cube", self.d):
            self.assertAlmostEqual(ged.vec_length(location), self.d, places=12)

    def test_all_parcel_dipoles_zero(self):
        for kind in ged.CLASSES:
            points = ged.parcel_positions(kind, self.d)
            for i in range(3):
                self.assertAlmostEqual(sum(point[i] for point in points), 0., places=12)

    def test_reject_unknown_source(self):
        with self.assertRaises(ValueError):
            ged.parcel_positions("not_a_pattern", self.d)

    def test_total_equivalent_mass_independent_formula(self):
        self.assertAlmostEqual(self.result["equivalent_mass_kg"] / self.M, 1., places=14)

    def test_point_mass_exact_x_reference(self):
        # WHAT: Formula independently evaluated for radial acceleration.
        # WHY:  A circular comparison to another GED function is insufficient.
        R = 200.
        got = ged.point_mass_field((R, 0., 0.), (0., 0., 0.), self.M, self.G)
        expected = -self.G * self.M / (R * R)
        self.assertAlmostEqual(got[0]/expected, 1., places=14)
        self.assertEqual(got[1:], (0., 0.))

    def test_point_mass_field_at_offaxis_point(self):
        got = ged.point_mass_field((3., 4., 0.), (0., 0., 0.), 2., 10.)
        self.assertAlmostEqual(got[0], -20*3/125, places=14)
        self.assertAlmostEqual(got[1], -20*4/125, places=14)

    def test_point_mass_at_source_rejected(self):
        with self.assertRaises(ValueError):
            ged.point_mass_field((0., 0., 0.), (0., 0., 0.), 2., 10.)

    def test_pair_x_axial_closed_form(self):
        r, a = 200., self.d
        expected = -.5*self.G*self.M*(1/(r-a)**2 + 1/(r+a)**2)
        got = ged.pattern_field((r, 0., 0.), "pair_x", self.config)[0]
        self.assertAlmostEqual(got/expected, 1., places=14)

    def test_pair_y_transverse_closed_form(self):
        r, a = 200., self.d
        expected = -self.G*self.M*r/(r*r+a*a)**1.5
        got = ged.pattern_field((r, 0., 0.), "pair_y", self.config)[0]
        self.assertAlmostEqual(got/expected, 1., places=14)

    def test_cube_axial_independent_closed_form(self):
        r, a = 200., self.d/math.sqrt(3.)
        expected = -.5*self.G*self.M*((r-a)/(((r-a)**2 + 2*a*a)**1.5) +
                                        (r+a)/(((r+a)**2 + 2*a*a)**1.5))
        got = ged.pattern_field((r, 0., 0.), "cube", self.config)[0]
        self.assertAlmostEqual(got/expected, 1., places=13)

    def test_pair_x_axial_field_stronger_than_pair_y(self):
        x = (200., 0., 0.)
        self.assertGreater(abs(ged.pattern_field(x, "pair_x", self.config)[0]),
                           abs(ged.pattern_field(x, "pair_y", self.config)[0]))

    def test_pattern_field_rejects_internal_probe(self):
        with self.assertRaises(ValueError):
            ged.pattern_field((90., 0., 0.), "cube", self.config)

    def test_symmetry_zero_transverse_components(self):
        for kind in ged.CLASSES:
            self.assertEqual(ged.pattern_field((200., 0., 0.), kind, self.config)[1:], (0., 0.))

    def test_cube_differs_from_point_monopole(self):
        g = ged.pattern_field((200., 0., 0.), "cube", self.config)[0]
        monopole = -self.G*self.M / 200**2
        self.assertNotEqual(g, monopole)

    def test_far_field_monopole_limit(self):
        r = 1e7
        expected = -self.G*self.M / r**2
        for kind in ged.CLASSES:
            actual = ged.pattern_field((r, 0., 0.), kind, self.config)[0]
            self.assertLess(abs(actual/expected-1), 1e-9)

    def test_source_signatures_seventy_two_components(self):
        self.assertTrue(all(len(sig) == 72 for sig in self.library.values()))

    def test_templates_are_distinct(self):
        self.assertEqual(len(set(self.library.values())), 3)

    def test_each_noiseless_truth_class_recovered(self):
        for kind in ged.CLASSES:
            predicted, _ = ged.classify(self.library[kind], self.library)
            self.assertEqual(predicted, kind)

    def test_matching_residual_zero(self):
        for kind in ged.CLASSES:
            _, residuals = ged.classify(self.library[kind], self.library)
            self.assertEqual(residuals[kind], 0.)

    def test_nonmatching_residual_positive(self):
        for kind in ged.CLASSES:
            _, residuals = ged.classify(self.library[kind], self.library)
            self.assertTrue(all(residuals[other] > 0 for other in ged.CLASSES if other != kind))

    def test_template_distance_symmetric(self):
        a, b = self.library["pair_x"], self.library["cube"]
        self.assertEqual(ged.squared_distance(a, b), ged.squared_distance(b, a))

    def test_squared_distance_zero_self(self):
        self.assertEqual(ged.squared_distance(self.library["pair_y"], self.library["pair_y"]), 0.)

    def test_squared_distance_rejects_wrong_length(self):
        with self.assertRaises(ValueError):
            ged.squared_distance((1.,), (1., 2.))

    def test_template_separation_positive(self):
        self.assertTrue(all(value > 0 for value in self.result["template_separation_l2_m_s2"].values()))

    def test_zero_noise_is_no_change(self):
        self.assertEqual(ged.add_noise((1., 2.), 0., random.Random(3)), (1., 2.))

    def test_repeated_random_seed_matches(self):
        first = ged.add_noise(self.library["pair_x"], 1e-24, random.Random(42))
        second = ged.add_noise(self.library["pair_x"], 1e-24, random.Random(42))
        self.assertEqual(first, second)

    def test_zero_noise_accuracy_all_classes(self):
        self.assertEqual(self.result["noise_sweep"][0]["accuracy"], 1.)

    def test_low_noise_accuracy_all_classes(self):
        self.assertEqual(self.result["noise_sweep"][1]["accuracy"], 1.)

    def test_high_noise_much_worse_than_low_noise(self):
        self.assertLess(self.result["noise_sweep"][-1]["accuracy"], .7)

    def test_high_noise_not_zero_success(self):
        self.assertGreater(self.result["noise_sweep"][-1]["accuracy"], .1)

    def test_noise_trials_per_class(self):
        trials = self.config["trials_per_class_per_noise_level"]
        for sweep in self.result["noise_sweep"]:
            for truth in ged.CLASSES:
                self.assertEqual(sum(sweep["confusion_matrix"][truth].values()), trials)

    def test_accuracy_matches_confusion_matrix(self):
        for sweep in self.result["noise_sweep"]:
            correct = sum(sweep["confusion_matrix"][k][k] for k in ged.CLASSES)
            self.assertEqual(correct, sweep["correct"])
            self.assertAlmostEqual(correct/sweep["total"], sweep["accuracy"])

    def test_noise_sweep_repeatable(self):
        self.assertEqual(ged.monte_carlo(self.config, self.library), self.result["noise_sweep"])

    def test_shell_30_and_70_exterior_exact(self):
        self.assertTrue(self.result["shell_null_control"]["external_signatures_equal"])
        self.assertEqual(self.result["shell_null_control"]["external_max_absolute_difference_m_s2"], 0.)

    def test_shell_field_outside_matches_monopole(self):
        p = (200., 0., 0.)
        field = ged.spherical_shell_field(p, 30., self.M, self.G)
        self.assertEqual(field, ged.point_mass_field(p, (0., 0., 0.), self.M, self.G))

    def test_shell_interior_zero(self):
        self.assertEqual(ged.spherical_shell_field((50., 0., 0.), 70., self.M, self.G), (0., 0., 0.))

    def test_inner_shell_at_50_nonzero(self):
        self.assertLess(ged.spherical_shell_field((50., 0., 0.), 30., self.M, self.G)[0], 0.)

    def test_shell_surface_rejected(self):
        with self.assertRaises(ValueError):
            ged.spherical_shell_field((30., 0., 0.), 30., self.M, self.G)

    def test_shell_inside_probe_distinguishes_arrangements(self):
        self.assertNotEqual(self.result["shell_null_control"]["inside_field_a_x_m_s2"],
                            self.result["shell_null_control"]["inside_field_b_x_m_s2"])

    def test_shell_equality_all_exterior_directions(self):
        for point in ged.probe_locations(self.config):
            a = ged.spherical_shell_field(point, 30., self.M, self.G)
            b = ged.spherical_shell_field(point, 70., self.M, self.G)
            self.assertEqual(a, b)

    def test_invalid_experiment_id_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("experiment_id", "GED-V1-R7"))

    def test_invalid_classes_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("candidate_classes", ["cube", "pair_x", "pair_y"]))

    def test_invalid_negative_energy_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("added_energy_j", -1.))

    def test_invalid_boolean_G_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("gravitational_constant_m3_kg_s2", True))

    def test_invalid_nan_speed_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("speed_of_light_m_s", float("nan")))

    def test_invalid_parcel_outside_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("parcel_offset_radius_m", 120.))

    def test_invalid_shell_geometry_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("shell_radius_a_m", 60.))

    def test_invalid_unsorted_probe_radii_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("exterior_probe_radii_m", [300., 200.]))

    def test_invalid_probe_inside_body_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("exterior_probe_radii_m", [50., 200.]))

    def test_invalid_duplicate_noise_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("gaussian_noise_std_m_s2", [0., 1e-26, 1e-26]))

    def test_invalid_negative_noise_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("gaussian_noise_std_m_s2", [0., -1e-26]))

    def test_invalid_boolean_trials_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("trials_per_class_per_noise_level", True))

    def test_invalid_zero_trials_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("trials_per_class_per_noise_level", 0))

    def test_invalid_negative_random_seed_rejected(self):
        with self.assertRaises(ValueError):
            ged.validate(self.with_change("random_seed", -1))

    def test_disclosure_required(self):
        altered = copy.deepcopy(self.config)
        altered["claims"]["observational_dataset"] = True
        with self.assertRaises(ValueError):
            ged.validate(altered)

    def test_report_mentions_exact_nonuniqueness(self):
        self.assertIn("EXACT EXTERIOR NONUNIQUENESS", ged.write_report(self.result))

    def test_report_disclaims_new_gravity(self):
        self.assertIn("No new gravity channel", ged.write_report(self.result))

    def test_result_json_serializable_without_nan(self):
        self.assertIn('"experiment_id": "GED-V1-R8"', json.dumps(self.result, allow_nan=False))


if __name__ == "__main__":
    unittest.main(verbosity=2)
