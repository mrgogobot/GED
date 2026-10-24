#!/usr/bin/env python3
"""GED-V1-R9: dimensional, algebraic, physical and provenance regression tests.

WHY: Test identities independently where possible, test failures explicitly,
     and distinguish forward-problem correctness from inverse identifiability.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import random
import unittest

import ged_v1_r9 as ged

BASE=Path(__file__).resolve().parent
CONFIG=json.loads((BASE/ged.CONFIG_NAME).read_text())


class GEDR9Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=CONFIG
        cls.truth=CONFIG["unknown_truth_energy_fractions"]
        cls.sites=ged.source_sites(CONFIG["source_offset_m"])
        cls.all_matrix=ged.measurement_matrix(CONFIG["exterior_probes_m"],CONFIG["source_offset_m"],CONFIG["reference_radius_m"])
        cls.result=ged.sampling_experiment(CONFIG)

    def revised(self, key, value):
        updated=copy.deepcopy(self.c)
        updated[key]=value
        return updated

    # WHAT: Lock configuration and scientific-claim boundaries.
    def test_experiment_identifier(self):
        self.assertEqual(ged.EXPERIMENT_ID,"GED-V1-R9")

    def test_config_is_valid(self):
        ged.validate(self.c)

    def test_six_unknown_energy_weights(self):
        self.assertEqual(len(self.truth),6)

    def test_energy_weights_sum_unity(self):
        self.assertAlmostEqual(math.fsum(self.truth),1.0,places=14)

    def test_truth_fractions_strictly_positive(self):
        self.assertTrue(all(x>0 for x in self.truth))

    def test_no_observation_claim(self):
        self.assertFalse(self.c["scientific_claims"]["observation"])

    def test_no_new_gravity_claim(self):
        self.assertFalse(self.c["scientific_claims"]["new_gravity"])

    def test_source_positions_are_assumed(self):
        self.assertTrue(self.c["scientific_claims"]["source_positions_assumed"])

    def test_background_subtraction_is_idealized(self):
        self.assertTrue(self.c["scientific_claims"]["background_subtraction_assumed_exact"])

    def test_full_gr_is_not_claimed(self):
        self.assertFalse(self.c["scientific_claims"]["full_gr"])

    def test_reject_wrong_id(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("experiment_id","GED-V1-R8"))

    def test_reject_boolean_energy(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("total_added_energy_j",True))

    def test_reject_negative_gravitational_constant(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("G_m3_kg_s2",-1))

    def test_reject_zero_speed_of_light(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("c_m_s",0))

    def test_reject_nan_source_radius(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("source_offset_m",float("nan")))

    def test_reject_source_outside_body(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("source_offset_m",150))

    def test_reject_unsupported_truth_length(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("unknown_truth_energy_fractions",[0.5,0.5]))

    def test_reject_nonunit_truth(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("unknown_truth_energy_fractions",[0.1]*6))

    def test_reject_negative_truth_weight(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("unknown_truth_energy_fractions",[0.35,-0.08,0.19,0.14,0.11,0.29]))

    def test_reject_boolean_truth_weight(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("unknown_truth_energy_fractions",[True,0.08,0.19,0.14,0.11,0.13]))

    def test_reject_duplicate_exterior_probes(self):
        points=copy.deepcopy(self.c["exterior_probes_m"]); points[1]=points[0]
        with self.assertRaises(ValueError):
            ged.validate(self.revised("exterior_probes_m",points))

    def test_reject_interior_probe(self):
        points=copy.deepcopy(self.c["exterior_probes_m"]); points[0]=[20,0,0]
        with self.assertRaises(ValueError):
            ged.validate(self.revised("exterior_probes_m",points))

    def test_reject_duplicate_axis_radii(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("axis_probe_radii_m",[200,250,250,450,600,850]))

    def test_reject_unordered_shells(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("shell_radii_m",[45,20,75]))

    def test_reject_identical_shells(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("shell_energy_fractions",[[0.7,0.2,0.1]]*2))

    def test_reject_negative_noise(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("noise_sigma_in_reference_g",[0,-1,0.1]))

    def test_reject_unsorted_noise(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("noise_sigma_in_reference_g",[0,0.01,0.001]))

    def test_reject_boolean_trials(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("trials_per_level",True))

    def test_reject_negative_rng_seed(self):
        with self.assertRaises(ValueError):
            ged.validate(self.revised("random_seed",-3))

    def test_reject_claims_changed(self):
        change=copy.deepcopy(self.c["scientific_claims"]);change["new_gravity"]=True
        with self.assertRaises(ValueError):
            ged.validate(self.revised("scientific_claims",change))

    # WHAT: Test the forward Newtonian model and physical units.
    def test_norm_3_4_5(self):
        self.assertEqual(ged.norm([3,4,0]),5)

    def test_six_axis_source_sites(self):
        self.assertEqual(len(self.sites),6)

    def test_sources_opposite_pairs(self):
        for a,b in ((0,1),(2,3),(4,5)):
            self.assertEqual(tuple(-x for x in self.sites[a]),self.sites[b])

    def test_source_site_radii(self):
        for site in self.sites:
            self.assertEqual(ged.norm(site),50.)

    def test_all_probes_exterior(self):
        self.assertTrue(all(ged.norm(p)>100 for p in self.c["exterior_probes_m"]))

    def test_full_vector_matrix_dimensions(self):
        self.assertEqual(len(self.all_matrix),24)
        self.assertTrue(all(len(row)==6 for row in self.all_matrix))

    def test_scalar_matrix_dimensions(self):
        axis=[[r,0,0] for r in self.c["axis_probe_radii_m"]]
        mat=ged.measurement_matrix(axis,50,200,component=0)
        self.assertEqual((len(mat),len(mat[0])),(6,6))

    def test_invalid_component(self):
        with self.assertRaises(ValueError):
            ged.measurement_matrix(self.c["exterior_probes_m"],50,200,component=3)

    def test_coincident_probe_fails(self):
        with self.assertRaises(ValueError):
            ged.signature_column([[50,0,0]],[50,0,0],200)

    def test_forward_reject_bad_weights(self):
        with self.assertRaises(ValueError):
            ged.forward(self.all_matrix,[1])

    def test_reference_g_independent_formula(self):
        expected=6.67430e-11*(1e8/(299792458.**2))/(200.**2)
        self.assertAlmostEqual(ged.reference_acceleration(self.c),expected,places=38)

    def test_reference_units_scaling(self):
        value=0.2
        self.assertAlmostEqual(ged.scaled_gravity_m_s2(value,self.c),value*1.856540067279666e-24,delta=1e-38)

    def test_physical_gravity_independent_single_source(self):
        # WHY: Calculate one point's field directly from Newton, independent
        #      of the dimensionless measurement-matrix implementation.
        point=self.c["exterior_probes_m"][0]
        src=self.sites[0]
        dist=math.dist(point,src)
        mass=1e8/(299792458.**2)
        gx=-6.67430e-11*mass*(point[0]-src[0])/(dist**3)
        got=ged.scaled_gravity_m_s2(self.all_matrix[0][0],self.c)
        self.assertAlmostEqual(got,gx,delta=2e-39)

    def test_linear_superposition_identity(self):
        unit0=[1,0,0,0,0,0]; unit1=[0,1,0,0,0,0]
        a=ged.forward(self.all_matrix,unit0)
        b=ged.forward(self.all_matrix,unit1)
        mixed=ged.forward(self.all_matrix,[0.25,0.75,0,0,0,0])
        self.assertLess(max(abs(z-(.25*x+.75*y)) for x,y,z in zip(a,b,mixed)),1e-14)

    def test_monopole_scaling_inverse_square(self):
        p=[200.,0.,0.]
        g200=ged.field_for_shells(p,[20,45,75],[0.7,0.2,0.1],self.c)[0]
        g400=ged.field_for_shells([400.,0,0],[20,45,75],[0.7,0.2,0.1],self.c)[0]
        self.assertAlmostEqual(g200/g400,4,places=12)

    def test_shell_interior_zero(self):
        self.assertEqual(ged.field_for_shells([10,0,0],[20,45,75],[0.7,0.2,0.1],self.c),[0,0,0])

    def test_shell_origin_zero(self):
        self.assertEqual(ged.field_for_shells([0,0,0],[20,45,75],[0.7,0.2,0.1],self.c),[0,0,0])

    def test_shell_surface_rejected(self):
        with self.assertRaises(ValueError):
            ged.field_for_shells([45,0,0],[20,45,75],[0.7,0.2,0.1],self.c)

    def test_shell_external_exact_nonuniqueness(self):
        self.assertEqual(self.result["spherical_shells"]["max_exterior_dimensionless_difference"],0)

    def test_shell_inside_distinguishes(self):
        self.assertGreater(self.result["spherical_shells"]["interior_at_35m_dimensionless_difference"],1)

    # WHAT: Algebraic identifiability and nullspace tests.
    def test_rank_of_full_vector_geometry(self):
        self.assertEqual(len(ged.rref(ged.augment_total_known(self.all_matrix))[1]),6)

    def test_one_probe_augmented_rank_four(self):
        self.assertEqual(self.result["sparse_one_probe"]["augmented_rank"],4)

    def test_one_probe_nullity_two(self):
        self.assertEqual(self.result["sparse_one_probe"]["nullity"],2)

    def test_sparse_null_direction_nonzero(self):
        v=ged.null_vector(ged.augment_total_known(self.all_matrix[:3]))
        self.assertAlmostEqual(max(abs(x) for x in v),1.0,places=12)

    def test_sparse_null_direction_residual(self):
        mat=ged.augment_total_known(self.all_matrix[:3])
        v=ged.null_vector(mat)
        self.assertLess(max(abs(x) for x in ged.forward(mat,v)),1e-12)

    def test_sparse_alternative_nonnegative(self):
        alt=self.result["sparse_one_probe"]["alternative_fractions"]
        self.assertTrue(all(v>=0 for v in alt))

    def test_sparse_alternative_total_conserved(self):
        alt=self.result["sparse_one_probe"]["alternative_fractions"]
        self.assertAlmostEqual(math.fsum(alt),1.0,places=12)

    def test_sparse_alternative_not_identical_to_truth(self):
        self.assertGreater(self.result["sparse_one_probe"]["l1_share_difference"],0.2)

    def test_sparse_alternative_same_measurements(self):
        self.assertLess(self.result["sparse_one_probe"]["maximum_dimensionless_field_residual"],1e-12)

    def test_axis_scalar_geometry_rank_three(self):
        self.assertEqual(self.result["radial_scalar"]["augmented_rank"],3)

    def test_axis_only_never_full_rank_here(self):
        self.assertGreater(self.result["radial_scalar"]["nullity"],0)

    def test_axis_plus_y_equals_minus_y_columns(self):
        axis=[[r,0,0] for r in self.c["axis_probe_radii_m"]]
        mat=ged.measurement_matrix(axis,50,200,component=0)
        self.assertTrue(all(row[2]==row[3] for row in mat))

    def test_axis_y_vs_z_columns_indistinguishable(self):
        axis=[[r,0,0] for r in self.c["axis_probe_radii_m"]]
        mat=ged.measurement_matrix(axis,50,200,component=0)
        self.assertTrue(all(row[2]==row[4] for row in mat))

    def test_rank_must_not_decrease_as_probes_added(self):
        ranks=[row["rank_with_total"] for row in self.result["distributed_geometry"]]
        self.assertEqual(ranks,sorted(ranks))

    def test_two_well_placed_vector_probes_full_rank(self):
        self.assertEqual(self.result["distributed_geometry"][1]["rank_with_total"],6)

    def test_inverse_rref_full_rank_has_no_nullvector(self):
        with self.assertRaises(ValueError):
            ged.null_vector(ged.augment_total_known(self.all_matrix))

    def test_rref_simple_independent_matrix(self):
        reduced,pivots=ged.rref([[1,0,0],[0,1,0]])
        self.assertEqual(pivots,[0,1])

    def test_rref_reject_jagged_matrix(self):
        with self.assertRaises(ValueError):
            ged.rref([[1,2],[3]])

    def test_solve_square_independent_identity(self):
        self.assertEqual(ged.linear_solve([[1,0],[0,1]],[3,4]),[3,4])

    def test_solve_square_pivoting(self):
        sol=ged.linear_solve([[0,1],[1,0]],[4,3])
        self.assertEqual(sol,[3,4])

    def test_solve_square_singular_rejected(self):
        with self.assertRaises(ValueError):
            ged.linear_solve([[1,1],[2,2]],[1,2])

    def test_fitting_simplex(self):
        fit,_=ged.simplex_least_squares(self.all_matrix,ged.forward(self.all_matrix,self.truth))
        self.assertGreaterEqual(min(fit),0)
        self.assertAlmostEqual(math.fsum(fit),1,places=12)

    def test_noiseless_full_inverse_precision(self):
        self.assertLess(self.result["full_vector_noiseless"]["l1_share_error"],1e-9)

    def test_noiseless_solution_is_not_given_to_optimizer(self):
        first=[0.1,0.2,0.1,0.3,0.2,0.1]
        fit,_=ged.simplex_least_squares(self.all_matrix,ged.forward(self.all_matrix,first))
        self.assertLess(math.fsum(abs(a-b) for a,b in zip(fit,first)),1e-9)

    def test_inversion_residual_nonnegative(self):
        self.assertGreaterEqual(self.result["full_vector_noiseless"]["residual"],0)

    def test_optimizer_reject_bad_measurement_length(self):
        with self.assertRaises(ValueError):
            ged.simplex_least_squares(self.all_matrix,[1])

    # WHAT: Seeded Monte Carlo uncertainties; these are simulations, not instruments.
    def test_seeded_sweep_repeatable(self):
        repeat=ged.sampling_experiment(self.c)
        self.assertEqual(repeat["synthetic_noise_sweep"],self.result["synthetic_noise_sweep"])

    def test_noise_levels_preserved(self):
        self.assertEqual([x["noise_sigma_in_reference_g"] for x in self.result["synthetic_noise_sweep"]],self.c["noise_sigma_in_reference_g"])

    def test_noise_trial_count(self):
        self.assertTrue(all(row["trials"]==100 for row in self.result["synthetic_noise_sweep"]))

    def test_noise_zero_recovers_truth(self):
        self.assertLess(self.result["synthetic_noise_sweep"][0]["mean_l1_share_error"],1e-9)

    def test_low_noise_better_than_high_noise(self):
        rows=self.result["synthetic_noise_sweep"]
        self.assertLess(rows[1]["mean_l1_share_error"],rows[-1]["mean_l1_share_error"])

    def test_noise_scale_in_si_units(self):
        row=self.result["synthetic_noise_sweep"][2]
        self.assertAlmostEqual(row["noise_sigma_m_s2"],0.001*ged.reference_acceleration(self.c),delta=1e-40)

    def test_noise_high_has_nonzero_error(self):
        self.assertGreater(self.result["synthetic_noise_sweep"][-1]["mean_l1_share_error"],0.3)

    def test_noise_low_has_small_error(self):
        self.assertLess(self.result["synthetic_noise_sweep"][1]["mean_l1_share_error"],0.01)

    def test_simulated_rather_than_observed(self):
        self.assertTrue(self.result["synthetic"])

    def test_model_limits_present(self):
        self.assertGreaterEqual(len(self.result["limits"]),8)

    def test_sha256_file_matches_hashlib(self):
        p=BASE/ged.CONFIG_NAME
        self.assertEqual(ged.sha256_file(p),hashlib.sha256(p.read_bytes()).hexdigest())

    def test_required_source_files_present(self):
        self.assertTrue(all((BASE/f).is_file() for f in ged.INPUT_FILES))


if __name__ == "__main__":
    unittest.main()
