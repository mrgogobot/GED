"""GED-V1-R4 verification gates — physical controls and edge conditions.

WHAT: Check arithmetic, units, conservation, invariants and interpretation.
WHY: A passing implementation test must not be mistaken for physical evidence.
The tests deliberately include independent analytic/float cross-checks.
"""
from __future__ import annotations
import copy
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import unittest

import ged_v1_r4 as ged

CFG = json.loads((Path(__file__).parent / 'GED_V1_R4_config.json').read_text())


def relative_close(test, got, target, tol='1e-32'):
    """WHAT: Compare model Decimals at high precision.
    WHY: exact residual zero is appropriate only for algebraic identities;
    nonterminating cube roots/ln require a finite tolerance.
    """
    with localcontext() as ctx:
        ctx.prec = 60
        lhs, rhs = ged.D(got), ged.D(target)
        scale = max(abs(lhs), abs(rhs), Decimal('1e-50'))
        test.assertLessEqual(abs(lhs-rhs)/scale, Decimal(tol))


class GEDV1R4Tests(unittest.TestCase):
    """WHAT: Regression suite against both intended results and misuse."""
    @classmethod
    def setUpClass(cls):
        cls.p = ged.read_parameters(CFG)
        cls.body = ged.spherical_model(cls.p)
        cls.result, cls.samples = ged.calculate(CFG)
        cls.A = cls.result['control_A_equal_local_stress_energy']
        cls.B = cls.result['control_B_thermal_redistribution']
        cls.C = cls.result['control_C_equal_entropy_different_interior_distribution']

    def test_experiment_is_not_observational(self):
        self.assertFalse(self.result['observation_performed'])

    def test_does_not_claim_new_gravity(self):
        self.assertFalse(self.result['new_gravity_channel_detected'])

    def test_half_volume_geometry(self):
        with localcontext() as ctx:
            ctx.prec = 60
            relative_close(self, (self.body['a']/self.body['R'])**3, Decimal('0.5'))

    def test_independent_spherical_volume(self):
        independent = (4/3) * math.pi * 100**3
        self.assertAlmostEqual(float(self.body['volume']), independent, places=8)

    def test_independent_reference_mass(self):
        self.assertAlmostEqual(float(self.body['reference_mass']),
                               3000 * 4/3 * math.pi * 100**3, places=4)

    def test_equal_zone_rest_mass(self):
        with localcontext() as ctx:
            ctx.prec = 60
            self.assertEqual(self.body['reference_mass'], 2*self.body['mass_zone'])

    def test_independent_heat_capacity(self):
        with localcontext() as ctx:
            ctx.prec = 60
            relative_close(self, self.body['C'], self.body['mass_zone']*Decimal(1000))

    def test_energy_mass_equivalence_independent(self):
        energy = Decimal('3.14e14')
        with localcontext() as ctx:
            ctx.prec = 60
            self.assertEqual(ged.r2.mass_equivalent(energy, self.p['c_m_s']),
                             energy / self.p['c_m_s']**2)

    def test_entropic_macrostate_aligned_zero(self):
        self.assertEqual(ged.microstate_entropy(256, 0), Decimal(0))

    def test_entropic_macrostate_many_configurations(self):
        self.assertGreater(ged.microstate_entropy(256, 128), Decimal(0))

    def test_microstate_multiplicity_independent(self):
        self.assertEqual(self.A['multiplicity_B'], str(math.comb(256, 128)))

    def test_microstate_entropy_float_crosscheck(self):
        expected = 1.380649e-23 * math.log(math.comb(256,128))
        self.assertAlmostEqual(float(self.A['entropy_B_minus_A_J_K']), expected,
                               delta=abs(expected)*1e-14)

    def test_macrostate_combinatorial_symmetry(self):
        self.assertEqual(ged.microstate_entropy(21, 4),
                         ged.microstate_entropy(21, 17))

    def test_control_A_zero_energy(self):
        self.assertEqual(self.A['energy_difference_J'], '0')

    def test_control_A_zero_gravity(self):
        self.assertEqual(self.A['gravitational_field_difference_m_s2'], '0')

    def test_control_A_nonzero_entropy_difference(self):
        self.assertGreater(ged.D(self.A['entropy_B_minus_A_J_K']), 0)

    def test_thermal_initial_hot_core(self):
        self.assertEqual(self.B['initial']['T_core_K'], '300')

    def test_thermal_initial_cold_shell(self):
        self.assertEqual(self.B['initial']['T_shell_K'], '200')

    def test_thermal_final_uniform(self):
        self.assertEqual(self.B['relaxed']['T_core_K'], '250')
        self.assertEqual(self.B['relaxed']['T_shell_K'], '250')

    def test_thermal_energy_initial_zero_relative_reference(self):
        self.assertEqual(ged.D(self.B['initial']['dE_total_J']), Decimal(0))

    def test_thermal_energy_final_zero_relative_reference(self):
        self.assertEqual(ged.D(self.B['relaxed']['dE_total_J']), Decimal(0))

    def test_energy_transfer_conserved(self):
        self.assertEqual(ged.D(self.B['transfer_energy_balance_J']), Decimal(0))

    def test_heat_stays_inside(self):
        self.assertEqual(self.B['energy_exported_J'], '0')

    def test_core_loses_thermal_energy_on_relaxation(self):
        self.assertLess(ged.D(self.B['core_thermal_energy_transfer_J']), 0)

    def test_shell_gains_same_energy(self):
        self.assertEqual(ged.D(self.B['shell_thermal_energy_transfer_J']),
                         ged.D(self.B['core_thermal_energy_transfer_J']).copy_negate())

    def test_entropy_increases_with_thermal_relaxation(self):
        self.assertGreater(ged.D(self.B['entropy_produced_in_relaxation_J_K']), 0)

    def test_entropy_analytic_independent(self):
        # WHAT: Floating-point law C ln(250²/(300*200)) independently.
        # WHY: Avoid a test that simply calls the same production helper.
        C = float(self.body['C'])
        expected = C * math.log((250*250)/(300*200))
        actual = float(self.B['entropy_produced_in_relaxation_J_K'])
        self.assertAlmostEqual(actual, expected, delta=abs(expected)*1e-13)

    def test_mass_same_at_both_endpoints(self):
        self.assertEqual(ged.D(self.B['model_exterior_mass_difference_kg']), Decimal(0))

    def test_exterior_g_same_at_both_endpoints(self):
        self.assertEqual(ged.D(self.B['exterior_g_delta_at_probe_m_s2']), Decimal(0))

    def test_interior_field_changes(self):
        self.assertLess(ged.D(self.B['interior_g_delta_at_core_boundary_m_s2']), 0)

    def test_interior_delta_g_independent(self):
        # At core boundary only the hot core thermal excess contributes.
        with localcontext() as ctx:
            ctx.prec = 60
            transferred = -self.body['C']*Decimal('50')
            expected = self.p['G_m3_kg_s2'] * transferred / (
                self.p['c_m_s']**2*self.body['a']**2)
            relative_close(self, self.B['interior_g_delta_at_core_boundary_m_s2'],
                           expected)

    def test_newtonian_exterior_analytic(self):
        with localcontext() as ctx:
            ctx.prec = 60
            r = self.body['R']*Decimal(2)
            independent_g = (self.p['G_m3_kg_s2']*
                             self.body['reference_mass']/r**2)
            actual = ged.radial_g(ged.thermal_state(self.body, self.p, 300, 200),
                                  self.body, r, self.p['G_m3_kg_s2'])
            relative_close(self, actual, independent_g)

    def test_field_center_zero(self):
        state = ged.thermal_state(self.body, self.p, 300, 200)
        self.assertEqual(ged.radial_g(state, self.body, 0,
                                      self.p['G_m3_kg_s2']), Decimal(0))

    def test_field_continuity_at_core_boundary(self):
        # In the limit from either piece, M_enclosed tends to M_core.
        s = ged.thermal_state(self.body, self.p, 300, 200)
        self.assertEqual(ged.enclosed_mass(s, self.body, self.body['a']),
                         s['M_core_kg'])

    def test_field_at_outer_boundary(self):
        s = ged.thermal_state(self.body, self.p, 300, 200)
        self.assertEqual(ged.enclosed_mass(s, self.body, self.body['R']),
                         s['M_total_kg'])

    def test_opposite_gradients_equal_entropy(self):
        self.assertEqual(ged.D(self.C['entropy_difference_reverse_minus_initial_J_K']), Decimal(0))

    def test_opposite_gradients_equal_total_energy(self):
        self.assertEqual(ged.D(self.C['energy_difference_reverse_minus_initial_J']), Decimal(0))

    def test_opposite_gradients_equal_exterior_mass(self):
        self.assertEqual(ged.D(self.C['model_exterior_mass_difference_kg']), Decimal(0))

    def test_opposite_gradients_different_interior_g(self):
        self.assertLess(ged.D(self.C['core_boundary_g_difference_reverse_minus_initial_m_s2']), 0)

    def test_reverse_gradient_is_twice_individual_interior_change(self):
        with localcontext() as ctx:
            ctx.prec = 60
            relative_close(self, self.C['core_boundary_g_difference_reverse_minus_initial_m_s2'],
                           2*ged.D(self.B['interior_g_delta_at_core_boundary_m_s2']))

    def test_samples_include_interior_and_exterior(self):
        self.assertEqual(len(self.samples), 4)
        self.assertEqual(ged.D(self.samples[-1]['delta_g_uniform_minus_gradient_m_s2']), Decimal(0))

    def test_reject_zero_temperature(self):
        with self.assertRaises(ValueError):
            ged.thermal_state(self.body, self.p, 0, 250)

    def test_reject_negative_radius(self):
        s = ged.thermal_state(self.body, self.p, 250, 250)
        with self.assertRaises(ValueError):
            ged.enclosed_mass(s, self.body, -1)

    def test_reject_invalid_G(self):
        s = ged.thermal_state(self.body, self.p, 250, 250)
        with self.assertRaises(ValueError):
            ged.radial_g(s, self.body, 100, 0)

    def test_reject_zero_temperature_offset(self):
        bad = copy.deepcopy(CFG)
        bad['temperature_offset_K'] = '0'
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_nonpositive_low_temperature(self):
        bad = copy.deepcopy(CFG)
        bad['temperature_offset_K'] = '250'
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_unequal_zone_volume(self):
        bad = copy.deepcopy(CFG)
        bad['core_volume_fraction'] = '0.4'
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_inner_probe_outside(self):
        bad = copy.deepcopy(CFG)
        bad['probe_interior_fraction_R'] = '1.5'
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_outer_probe_inside(self):
        bad = copy.deepcopy(CFG)
        bad['probe_exterior_multiple_R'] = '0.7'
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_nonfinite_density(self):
        bad = copy.deepcopy(CFG)
        bad['density_kg_m3'] = 'NaN'
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_noninteger_site_count(self):
        bad = copy.deepcopy(CFG)
        bad['degenerate_site_count'] = 256.5
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_invalid_excitation(self):
        bad = copy.deepcopy(CFG)
        bad['macrostate_B_label_one_count'] = 300
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_reject_identical_macrostate_counts(self):
        bad = copy.deepcopy(CFG)
        bad['macrostate_B_label_one_count'] = 0
        with self.assertRaises(ValueError):
            ged.calculate(bad)

    def test_limits_mark_full_gr_not_solved(self):
        self.assertTrue(any('exact GR' in x for x in self.result['limitations']))

    def test_pinned_r2_sha256_matches(self):
        self.assertEqual(ged.sha256(Path(ged.r2.__file__)),
                         CFG['r2_dependency_sha256'])


if __name__ == '__main__':
    unittest.main()
