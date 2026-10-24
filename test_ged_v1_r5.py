"""GED-V1-R5 verification suite — independent physics identities and edge cases.

WHAT: Check numerical consistency, fixed-energy material control, units,
loading/unloading energy ledgers, and provenance dependency pin.
WHY: A passing model should not rely on its own output alone for validation.
"""
import copy
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import unittest

import ged_v1_r5 as m


class GEDV1R5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = m.load_config()
        cls.out = m.calculate(cls.cfg)
        cls.iron = cls.cfg['materials']['polycrystalline_iron']
        cls.ice = cls.cfg['materials']['polycrystalline_ice']

    def assertRel(self, a, b, rel='1e-52'):
        """WHAT: Decimal relative tolerance. WHY: sqrt is rounded at 60 digits."""
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            x, y = m.D(a), m.D(b)
            self.assertLessEqual(abs(x - y), m.D(rel) * max(abs(y), Decimal('1e-60')))

    # ---- Model assumptions and canonical constants. ----
    def test_experiment_id(self):
        self.assertEqual(self.out['experiment_id'], 'GED-V1-R5')

    def test_is_numerical_not_observed(self):
        self.assertFalse(self.out['observational_data'])

    def test_no_new_coupling_asserted(self):
        self.assertEqual(self.out['novel_information_channel'], 'NOT_DERIVED_OR_TESTED')

    def test_pinned_r2_code_sha256(self):
        self.assertEqual(m.check_dependency(), self.cfg['r2_dependency_sha256'])

    def test_moduli_are_distinct_and_positive(self):
        self.assertGreater(m.D(self.iron['young_modulus_Pa']), m.D(self.ice['young_modulus_Pa']))

    def test_elastic_units(self):
        # [Pa*m^3] = [N*m] = [J], used by each elasticity identity.
        self.assertEqual(m.D('1e6') * m.D('1e-6'), Decimal('1'))

    # ---- Direct Hookean independent equal-stress and equal-strain checks. ----
    def test_equal_stress_matches_input(self):
        for row in self.out['cases']['equal_applied_stress'].values():
            self.assertEqual(row['stress_Pa'], m.D(self.cfg['equal_stress_Pa']))

    def test_equal_strain_matches_input(self):
        for row in self.out['cases']['equal_imposed_strain'].values():
            self.assertEqual(row['strain_dimensionless'], m.D(self.cfg['equal_strain_dimensionless']))

    def test_equal_stress_iron_energy_independent(self):
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            Y = m.D(self.iron['young_modulus_Pa'])
            sigma = m.D(self.cfg['equal_stress_Pa'])
            V = m.D(self.cfg['sample_volume_m3'])
            independent = sigma**2 * V / (2 * Y)
            actual = self.out['cases']['equal_applied_stress']['polycrystalline_iron']['elastic_energy_J']
            self.assertRel(actual, independent)

    def test_equal_stress_ice_energy_independent(self):
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            Y = m.D(self.ice['young_modulus_Pa'])
            sigma = m.D(self.cfg['equal_stress_Pa'])
            V = m.D(self.cfg['sample_volume_m3'])
            independent = sigma**2 * V / (2 * Y)
            actual = self.out['cases']['equal_applied_stress']['polycrystalline_ice']['elastic_energy_J']
            self.assertRel(actual, independent)

    def test_equal_strain_iron_energy_independent(self):
        Y = m.D(self.iron['young_modulus_Pa'])
        strain = m.D(self.cfg['equal_strain_dimensionless'])
        V = m.D(self.cfg['sample_volume_m3'])
        expected = Y * strain * strain * V / 2
        actual = self.out['cases']['equal_imposed_strain']['polycrystalline_iron']['elastic_energy_J']
        self.assertRel(actual, expected)

    def test_equal_strain_ice_energy_independent(self):
        Y = m.D(self.ice['young_modulus_Pa'])
        strain = m.D(self.cfg['equal_strain_dimensionless'])
        V = m.D(self.cfg['sample_volume_m3'])
        expected = Y * strain * strain * V / 2
        actual = self.out['cases']['equal_imposed_strain']['polycrystalline_ice']['elastic_energy_J']
        self.assertRel(actual, expected)

    def test_fixed_stress_soft_ice_strains_more(self):
        data = self.out['cases']['equal_applied_stress']
        self.assertGreater(data['polycrystalline_ice']['strain_dimensionless'],
                           data['polycrystalline_iron']['strain_dimensionless'])

    def test_fixed_stress_soft_ice_stores_more(self):
        data = self.out['cases']['equal_applied_stress']
        self.assertGreater(data['polycrystalline_ice']['elastic_energy_J'],
                           data['polycrystalline_iron']['elastic_energy_J'])

    def test_fixed_strain_iron_stress_higher(self):
        data = self.out['cases']['equal_imposed_strain']
        self.assertGreater(data['polycrystalline_iron']['stress_Pa'],
                           data['polycrystalline_ice']['stress_Pa'])

    def test_fixed_strain_iron_stores_more(self):
        data = self.out['cases']['equal_imposed_strain']
        self.assertGreater(data['polycrystalline_iron']['elastic_energy_J'],
                           data['polycrystalline_ice']['elastic_energy_J'])

    def test_energy_scales_with_stress_squared(self):
        low = m.elastic_state('10000000000', '1', stress_Pa='1000')
        high = m.elastic_state('10000000000', '1', stress_Pa='2000')
        self.assertEqual(high['elastic_energy_J'], 4 * low['elastic_energy_J'])

    def test_energy_scales_with_strain_squared(self):
        low = m.elastic_state('10000000000', '1', strain='0.000001')
        high = m.elastic_state('10000000000', '1', strain='0.000002')
        self.assertEqual(high['elastic_energy_J'], 4 * low['elastic_energy_J'])

    def test_energy_scales_with_volume(self):
        a = m.elastic_state('10000000000', '1', strain='0.000001')
        b = m.elastic_state('10000000000', '3', strain='0.000001')
        self.assertEqual(b['elastic_energy_J'], 3 * a['elastic_energy_J'])

    # ---- Matched-energy control isolates energy from stiffness. ----
    def test_target_energy_iron(self):
        row = self.out['cases']['equal_stored_elastic_energy']['polycrystalline_iron']
        self.assertRel(row['elastic_energy_J'], self.cfg['equal_elastic_energy_J'])

    def test_target_energy_ice(self):
        row = self.out['cases']['equal_stored_elastic_energy']['polycrystalline_ice']
        self.assertRel(row['elastic_energy_J'], self.cfg['equal_elastic_energy_J'])

    def test_matched_energy_yet_different_stresses(self):
        fixed = self.out['cases']['equal_stored_elastic_energy']
        self.assertNotEqual(fixed['polycrystalline_iron']['stress_Pa'],
                            fixed['polycrystalline_ice']['stress_Pa'])

    def test_matched_energy_yet_different_strains(self):
        fixed = self.out['cases']['equal_stored_elastic_energy']
        self.assertNotEqual(fixed['polycrystalline_iron']['strain_dimensionless'],
                            fixed['polycrystalline_ice']['strain_dimensionless'])

    def test_matched_energy_monopole_correction_same(self):
        fixed = self.out['cases']['equal_stored_elastic_energy']
        a = fixed['polycrystalline_iron']['delta_g_point_monopole_proxy_m_s2']
        b = fixed['polycrystalline_ice']['delta_g_point_monopole_proxy_m_s2']
        self.assertRel(a, b)

    def test_baseline_gravity_differs_due_to_density(self):
        fixed = self.out['cases']['equal_stored_elastic_energy']
        self.assertNotEqual(fixed['polycrystalline_iron']['point_monopole_baseline_g_m_s2'],
                            fixed['polycrystalline_ice']['point_monopole_baseline_g_m_s2'])

    def test_mass_equivalent_independent_E_over_c_squared(self):
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            row = self.out['cases']['equal_imposed_strain']['polycrystalline_iron']
            independent = row['elastic_energy_J'] / m.D(self.cfg['c_m_s'])**2
            self.assertRel(row['delta_mass_if_energy_retained_kg'], independent)

    def test_gravity_increment_independent_formula(self):
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            row = self.out['cases']['equal_imposed_strain']['polycrystalline_ice']
            independent = (m.D(self.cfg['G_m3_kg_s2']) * row['elastic_energy_J'] /
                           (m.D(self.cfg['c_m_s'])**2 * m.D(self.cfg['far_field_probe_m'])**2))
            self.assertRel(row['delta_g_point_monopole_proxy_m_s2'], independent)

    def test_gravity_increment_inverse_square(self):
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            a = m.gravitational_proxy('10', '1000', '100', '6.67430e-11', '299792458')
            b = m.gravitational_proxy('10', '1000', '200', '6.67430e-11', '299792458')
            self.assertRel(a['delta_g_point_monopole_proxy_m_s2'],
                           4 * b['delta_g_point_monopole_proxy_m_s2'])

    def test_mass_increment_baseline_fraction(self):
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            row = self.out['cases']['equal_imposed_strain']['polycrystalline_iron']
            self.assertRel(row['fractional_mass_increment'],
                           row['delta_mass_if_energy_retained_kg'] / row['base_mass_kg'])

    # ---- Sample + reservoir + environment conservation controls. ----
    def test_energy_ledger_totals_zero(self):
        for case in self.out['cases'].values():
            for row in case.values():
                for ledger in row['energy_ledger'].values():
                    self.assertEqual(ledger['whole_system_balance_J'], 0)

    def test_loaded_sample_increases_energy(self):
        logs = m.energy_ledgers('125')
        self.assertEqual(logs['loaded_energy_stored']['sample_total_change_J'], Decimal(125))
        self.assertEqual(logs['loaded_energy_stored']['loading_reservoir_J'], Decimal(-125))

    def test_reversible_unload_restores_sample_and_reservoir(self):
        logs = m.energy_ledgers('125')
        self.assertEqual(logs['unloaded_reversibly']['sample_total_change_J'], 0)
        self.assertEqual(logs['unloaded_reversibly']['loading_reservoir_J'], 0)

    def test_dissipated_inside_keeps_sample_energy(self):
        logs = m.energy_ledgers('125')
        self.assertEqual(logs['dissipated_as_retained_heat']['sample_total_change_J'], Decimal(125))
        self.assertEqual(logs['dissipated_as_retained_heat']['sample_elastic_J'], 0)

    def test_exported_heat_goes_to_environment(self):
        logs = m.energy_ledgers('125')
        self.assertEqual(logs['heat_exported_to_environment']['sample_total_change_J'], 0)
        self.assertEqual(logs['heat_exported_to_environment']['environment_J'], Decimal(125))

    # ---- Rejection of scientifically invalid inputs. ----
    def test_reject_both_stress_and_strain(self):
        with self.assertRaises(ValueError):
            m.elastic_state('1e9', '1', stress_Pa='10', strain='0.000001')

    def test_reject_neither_stress_nor_strain(self):
        with self.assertRaises(ValueError):
            m.elastic_state('1e9', '1')

    def test_reject_overlarge_strain(self):
        with self.assertRaises(ValueError):
            m.elastic_state('1e9', '1', strain='0.01')

    def test_reject_negative_stress(self):
        with self.assertRaises(ValueError):
            m.elastic_state('1e9', '1', stress_Pa='-1')

    def test_reject_nonfinite(self):
        with self.assertRaises(ValueError):
            m.elastic_state('NaN', '1', strain='0.000001')

    def test_reject_boolean(self):
        with self.assertRaises(ValueError):
            m.elastic_state(True, '1', strain='0.000001')

    def test_reject_negative_target(self):
        with self.assertRaises(ValueError):
            m.state_for_target_energy('1e9', '1', '-1')

    def test_reject_negative_energy_ledger(self):
        with self.assertRaises(ValueError):
            m.energy_ledgers('-1')

    def test_reject_probe_too_close(self):
        wrong = copy.deepcopy(self.cfg)
        wrong['far_field_probe_m'] = '1'
        with self.assertRaises(ValueError):
            m.calculate(wrong)

    def test_reject_bad_config_version(self):
        wrong = copy.deepcopy(self.cfg)
        wrong['experiment_id'] = 'GED-V1-R4'
        with self.assertRaises(ValueError):
            m.calculate(wrong)

    def test_reject_bad_density(self):
        wrong = copy.deepcopy(self.cfg)
        wrong['materials']['polycrystalline_ice']['density_kg_m3'] = '0'
        with self.assertRaises(ValueError):
            m.calculate(wrong)

    def test_reject_invalid_gravitational_constant(self):
        with self.assertRaises(ValueError):
            m.gravitational_proxy('1', '1000', '100', '-1', '299792458')

    def test_reject_bad_dependency(self):
        # WHAT: Hash validation independently computed. WHY: block stale R2.
        expected = hashlib.sha256(Path(m.__file__).with_name('ged_v1_r2.py').read_bytes()).hexdigest()
        self.assertEqual(m.R2_SHA256, expected)


if __name__ == '__main__':
    unittest.main()
