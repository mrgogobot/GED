"""GED-V1-R10 regression + independent mathematics controls.

WHAT: Check shell-theorem forward physics, source-vs-inventory degeneracy,
      exact-rational inversion, input validation and export receipts.
WHY:  Code that produces plausible numbers must still withstand independent
      known-answer and adversarial tests, with failures retained.
"""
import copy
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

import ged_v1_r10 as m


class GEDR10Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg=m.load_configuration()
        cls.edges=[m.D(x) for x in cls.cfg["shell_boundaries_m"]]
        cls.A=[m.D(x) for x in cls.cfg["inventory_A"]["shell_total_mass_energy_kg"]]
        cls.B=[m.D(x) for x in cls.cfg["inventory_B"]["shell_total_mass_energy_kg"]]
        cls.C=[m.D(x) for x in cls.cfg["inventory_C"]["shell_total_mass_energy_kg"]]
        cls.G=m.D(cls.cfg["G_m3_kg_s2"])
        cls.c=m.D(cls.cfg["c_m_s"])
        cls.result, cls.rows=m.run_science(cls.cfg)

    # WHAT: Pin the experiment scope; WHY: never silently change what we're testing.
    def test_id(self): self.assertEqual(self.cfg['experiment_id'],'GED-V1-R10')
    def test_new_gravity_not_claimed(self): self.assertFalse(self.result['scientific_scope']['new_gravity_detected'])
    def test_remarks_not_observations(self): self.assertTrue(self.result['scientific_scope']['numerical_only'])
    def test_full_gr_not_claimed(self): self.assertFalse(self.result['scientific_scope']['general_relativity_solved'])
    def test_atomic_microphysics_not_claimed(self): self.assertFalse(self.result['scientific_scope']['real_material_microphysics_modeled'])
    def test_pressure_omission_stated(self): self.assertFalse(self.result['scientific_scope']['equilibrium_or_thermal_pressure_modeled'])
    def test_ledger_not_atomic_measurement(self): self.assertFalse(self.result['scientific_scope']['atomic_binding_measured'])
    def test_prior_evidence_document_present(self): self.assertTrue((m.ROOT/'GED_V1_R10_PRIOR_EVIDENCE.json').is_file())
    def test_prior_evidence_local_vs_audited(self):
        evidence=json.loads((m.ROOT/'GED_V1_R10_PRIOR_EVIDENCE.json').read_text())
        self.assertEqual(len(evidence['prior_experiments']),10)
        self.assertEqual(evidence['prior_experiments'][-1]['status'],'local_pass_pending_archive_audit')
    def test_inputs_present(self):
        for name in m.INPUT_NAMES: self.assertTrue((m.ROOT/name).is_file(),name)

    # WHAT: independently computed shell fractions; WHY: check actual geometry.
    def test_fraction_inside_first_shell(self): self.assertEqual(m.shell_fraction(m.D(10),m.D(0),m.D(20)),m.D('0.125'))
    def test_fraction_middle_shell_analytic(self):
        with localcontext() as ctx:
            ctx.prec=60
            self.assertEqual(m.shell_fraction(m.D(35),m.D(20),m.D(50)),m.D(34875)/m.D(117000))
    def test_fraction_before_shell(self): self.assertEqual(m.shell_fraction(m.D(10),m.D(20),m.D(50)),0)
    def test_fraction_after_shell(self): self.assertEqual(m.shell_fraction(m.D(60),m.D(20),m.D(50)),1)
    def test_fraction_boundary(self): self.assertEqual(m.shell_fraction(m.D(20),m.D(20),m.D(50)),0)
    def test_reject_negative_radius(self):
        with self.assertRaises(ValueError): m.shell_fraction(m.D(-1),m.D(0),m.D(20))
    def test_reject_inverted_shell(self):
        with self.assertRaises(ValueError): m.shell_fraction(m.D(3),m.D(20),m.D(10))
    def test_exact_fraction_inner(self): self.assertEqual(m.shell_fraction_exact(10,0,20),Fraction(1,8))
    def test_exact_fraction_mid(self): self.assertEqual(m.shell_fraction_exact(35,20,50),Fraction(34875,117000))
    def test_center_gravity_zero(self): self.assertEqual(m.gravity_magnitude(m.D(0),self.edges,self.A,self.G),0)
    def test_center_enclosed_mass_zero(self): self.assertEqual(m.enclosed_mass(m.D(0),self.edges,self.A),0)
    def test_total_mass_A(self): self.assertEqual(sum(self.A),m.D('1e10'))
    def test_total_mass_B(self): self.assertEqual(sum(self.B),sum(self.A))
    def test_A_B_different_distribution(self): self.assertNotEqual(self.A,self.B)
    def test_A_C_same_distribution(self): self.assertEqual(self.A,self.C)
    def test_enclosed_mass_r10_independent(self): self.assertEqual(m.enclosed_mass(m.D(10),self.edges,self.A),m.D('2.5e8'))
    def test_enclosed_mass_r10_B(self): self.assertEqual(m.enclosed_mass(m.D(10),self.edges,self.B),m.D('3.75e8'))
    def test_analytic_gravity_at_10(self): self.assertEqual(m.gravity_magnitude(m.D(10),self.edges,self.A,self.G),self.G*m.D('2.5e6'))
    def test_analytic_gravity_at_150(self): self.assertEqual(m.gravity_magnitude(m.D(150),self.edges,self.A,self.G),self.G*m.D('1e10')/m.D(22500))
    def test_exterior_equality_all_samples(self):
        for r in (100,150,300,400,1000):
            self.assertEqual(m.gravity_magnitude(m.D(r),self.edges,self.A,self.G),m.gravity_magnitude(m.D(r),self.edges,self.B,self.G))
    def test_interior_difference(self):
        self.assertNotEqual(m.gravity_magnitude(m.D(10),self.edges,self.A,self.G),m.gravity_magnitude(m.D(10),self.edges,self.B,self.G))
    def test_identity_degeneracy_all_radii(self):
        for r in (0,1,10,20,35,50,75,100,150,300):
            self.assertEqual(m.gravity_magnitude(m.D(r),self.edges,self.A,self.G),m.gravity_magnitude(m.D(r),self.edges,self.C,self.G))
    def test_mass_recomposition_from_binding(self):
        with localcontext() as ctx:
            ctx.prec=60
            ledger=m.constituent_energy_ledger(self.A,[m.D(x) for x in self.cfg['inventory_A']['shell_binding_energy_j']],self.c)
            for item, actual in zip(ledger,self.A):
                self.assertAlmostEqual(m.D(item['recovered_total_mass_equivalent_kg']),actual,delta=m.D('1e-40'))
    def test_bound_rest_masses_differ_A_C(self):
        self.assertNotEqual(self.result['case_inventory_ledgers']['inventory_A'][0]['constituent_rest_equivalent_kg'],
                            self.result['case_inventory_ledgers']['inventory_C'][0]['constituent_rest_equivalent_kg'])
    def test_binding_contribution_negative(self):
        self.assertTrue(all(m.D(x['binding_mass_equivalent_kg']) < 0 for x in self.result['case_inventory_ledgers']['inventory_A']))
    def test_ledger_size_mismatch_rejected(self):
        with self.assertRaises(ValueError): m.constituent_energy_ledger([m.D(1)],[],self.c)
    def test_exact_inverse_fraction_math(self):
        F=m.response_matrix_exact([10,35,75],[0,20,50,100])
        self.assertEqual(F[0],[Fraction(1,8),Fraction(0),Fraction(0)])
        self.assertEqual(F[1][2],0)
        self.assertEqual(m.matrix_rank_exact(F),3)
    def test_exterior_rank_one(self): self.assertEqual(m.matrix_rank_exact(m.response_matrix_exact([150,200,300],[0,20,50,100])),1)
    def test_exterior_nullity_two(self): self.assertEqual(self.result['exterior_nullity'],2)
    def test_interior_rank_three(self): self.assertEqual(self.result['interior_rank'],3)
    def test_interior_reconstruction_exact(self): self.assertTrue(self.result['reconstruction_exact'])
    def test_inverse_solution_expected(self): self.assertEqual(self.result['interior_recovered_shell_mass_kg'],['2000000000','3000000000','5000000000'])
    def test_all_exterior_rows_equal(self):
        self.assertEqual(m.response_matrix_exact([150,200,300],[0,20,50,100]),[[Fraction(1)]*3]*3)
    def test_nonunique_alternative_same_exterior_measurements(self):
        exterior=m.response_matrix_exact([150,200,300],[0,20,50,100])
        for row in exterior:
            self.assertEqual(sum(v*Fraction(x) for v,x in zip(row,self.A)),sum(v*Fraction(x) for v,x in zip(row,self.B)))
    def test_solve_known_integer_matrix(self):
        A=[[Fraction(1),Fraction(0)],[Fraction(1),Fraction(1)]]
        self.assertEqual(m.solve_square_exact(A,[Fraction(2),Fraction(7)]),[Fraction(2),Fraction(5)])
    def test_solve_singular_rejected(self):
        with self.assertRaises(ValueError): m.solve_square_exact([[Fraction(1),Fraction(1)],[Fraction(1),Fraction(1)]],[Fraction(2),Fraction(2)])
    def test_rank_zero_matrix(self): self.assertEqual(m.matrix_rank_exact([[Fraction(0)]*3]*2),0)
    def test_rank_independent_identity(self): self.assertEqual(m.matrix_rank_exact([[1,0,0],[0,1,0],[0,0,1]]),3)
    def test_reject_irregular_matrix(self):
        with self.assertRaises(ValueError): m.matrix_rank_exact([[1,2],[3]])
    def test_probe_samples_six(self): self.assertEqual(len(self.rows),6)
    def test_exterior_in_reports_zero(self): self.assertTrue(self.result['external_difference_exactly_zero'])
    def test_interior_report_diff(self): self.assertTrue(self.result['different_internal_field_A_B'])
    def test_A_C_identical_report(self): self.assertTrue(self.result['same_field_A_C_at_all_probes'])
    def test_unearned_rejected_in_verdict(self): self.assertFalse(self.result['hypothesis_assessment']['new_gravity_evidence'])
    def test_zero_gravity_without_source_mass(self): self.assertEqual(m.gravity_magnitude(m.D(150),self.edges,[m.D(0)]*3,self.G),0)
    def test_input_reject_wrong_id(self):
        c=copy.deepcopy(self.cfg);c['experiment_id']='GED-V1-R11'
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_nonfinite(self):
        c=copy.deepcopy(self.cfg);c['G_m3_kg_s2']='NaN'
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_wrong_shell(self):
        c=copy.deepcopy(self.cfg);c['shell_boundaries_m']=['0','40','20','100']
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_mismatched_inventory_mass(self):
        c=copy.deepcopy(self.cfg);c['inventory_B']['shell_total_mass_energy_kg'][0]='4000000000'
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_bad_binding_sign(self):
        c=copy.deepcopy(self.cfg);c['inventory_C']['shell_binding_energy_j'][0]='500000000'
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_identical_binding_ledgers(self):
        c=copy.deepcopy(self.cfg);c['inventory_C']['shell_binding_energy_j']=c['inventory_A']['shell_binding_energy_j'][:]
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_hidden_science_claim(self):
        c=copy.deepcopy(self.cfg);c['scientific_scope']['new_gravity_detected']=True
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_boolean_constant(self):
        c=copy.deepcopy(self.cfg);c['c_m_s']=True
        with self.assertRaises(ValueError): m.validate(c)
    def test_input_reject_duplicate_probes(self):
        c=copy.deepcopy(self.cfg);c['probe_radii_m']=['10','10']
        with self.assertRaises(ValueError): m.validate(c)
    def test_exact_fraction_reject_invalid(self):
        with self.assertRaises(ValueError): m.shell_fraction_exact(10,20,20)
    def test_source_g_ref_units(self): self.assertEqual(self.G*m.D('1e10')/m.D(150)**2, m.D('0.00002966355555555555555555555556'))
    def test_sha256_matches_standard_library(self):
        import hashlib
        self.assertEqual(m.sha256(m.ROOT/'ged_v1_r10.py'),hashlib.sha256((m.ROOT/'ged_v1_r10.py').read_bytes()).hexdigest())
    def test_result_is_json_serializable(self): json.dumps(self.result,allow_nan=False)
    def test_publication_is_not_part_of_run(self): self.assertFalse('publish' in m.run_science(self.cfg)[0])


if __name__=='__main__':
    unittest.main()
