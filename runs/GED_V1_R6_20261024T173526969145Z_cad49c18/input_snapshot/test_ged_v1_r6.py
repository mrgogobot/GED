"""GED-V1-R6 independent math, edge-case and provenance regression tests.

WHAT: Check analytic spherical-shell formulas, conservation, pinned dependencies,
input validation, elastic energy dimensions, and output provenance.
WHY: Passing a program's internal formulas alone cannot validate the physics;
this suite includes standalone algebra and explicitly marks model limitations.
"""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import ged_v1_r6 as ged


class GEDR6Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = ged.load_config()
        cls.params, cls.materials, cls.fractions, cls.radii = ged.validate(cls.cfg)
        cls.result = ged.calculate(cls.cfg)
        cls.geom = cls.result['geometry']

    def test_experiment_id(self):
        self.assertEqual(self.result['experiment_id'], 'GED-V1-R6')

    def test_no_claim_of_detection(self):
        self.assertFalse(self.result['observational_data'])
        self.assertEqual(self.result['novel_gravitational_information_channel'], 'NOT_DERIVED_OR_TESTED')

    def test_pinned_r2_exact_sha256(self):
        self.assertEqual(ged.sha256(ged.ROOT / 'ged_v1_r2.py'), ged.R2_SHA256)

    def test_pinned_r5_exact_sha256(self):
        self.assertEqual(ged.sha256(ged.ROOT / 'ged_v1_r5.py'), ged.R5_SHA256)

    def test_dependency_checker(self):
        self.assertEqual(ged.check_dependencies()['ged_v1_r5.py'], ged.R5_SHA256)

    def test_analytic_core_volume(self):
        with localcontext() as ctx:
            ctx.prec = 70
            expected = Decimal(4) * ged.r2.PI * Decimal(50)**3 / Decimal(3)
            self.assertLess(abs((self.geom['core_volume_m3'] - expected) / expected), Decimal('1e-65'))

    def test_shell_volume_seven_times_core(self):
        # Independently: 100^3 / 50^3 = 8, leaving 7 core volumes outside.
        with localcontext() as ctx:
            ctx.prec = 70
            difference = abs(self.geom['shell_volume_m3'] - self.geom['core_volume_m3'] * 7)
            self.assertLess(difference/self.geom['shell_volume_m3'], Decimal('1e-65'))

    def test_independent_total_rest_mass(self):
        with localcontext() as ctx:
            ctx.prec = 70
            expected = (Decimal(4) * ged.r2.PI / 3 *
                        (Decimal(7870) * Decimal(50)**3 +
                         Decimal(917) * (Decimal(100)**3 - Decimal(50)**3)))
            self.assertLess(abs((self.geom['total_rest_mass_kg'] - expected) / expected), Decimal('1e-65'))

    def test_equal_total_allocated_energy(self):
        for record in self.result['configs'].values():
            self.assertEqual(sum(v['allocated_energy_J'] for v in record.values()), Decimal('1e8'))

    def test_core_weighted_energy(self):
        self.assertEqual(self.result['configs']['core_weighted']['core']['allocated_energy_J'], Decimal('75000000'))

    def test_shell_weighted_energy(self):
        self.assertEqual(self.result['configs']['shell_weighted']['shell']['allocated_energy_J'], Decimal('75000000'))

    def test_each_case_has_same_material_inventory(self):
        self.assertEqual(len(self.result['configs']), 2)
        with localcontext() as ctx:
            ctx.prec = 70
            self.assertEqual(self.geom['total_rest_mass_kg'],
                             self.geom['core_rest_mass_kg'] + self.geom['shell_rest_mass_kg'])

    def test_energy_repartition_conserves(self):
        ledger = self.result['repartition_ledger_J']
        self.assertEqual(ledger['core'], Decimal('-50000000'))
        self.assertEqual(ledger['shell'], Decimal('50000000'))
        self.assertEqual(ledger['whole_body'], Decimal(0))

    def test_iron_strain_small(self):
        for case in self.result['configs'].values():
            self.assertLess(case['core']['strain_dimensionless'], Decimal('0.01'))

    def test_ice_strain_small(self):
        for case in self.result['configs'].values():
            self.assertLess(case['shell']['strain_dimensionless'], Decimal('0.01'))

    def test_elastic_density_in_pascal(self):
        for case in self.result['configs'].values():
            for region in ('core', 'shell'):
                # Pa = J/m^3; E*strain^2/2 should match energy density.
                strain = case[region]['strain_dimensionless']
                E = self.materials[region]['E_Pa']
                with localcontext() as ctx:
                    ctx.prec = 70
                    reconstructed = E * strain**2 / 2
                    expected = case[region]['energy_density_J_m3']
                    self.assertLess(abs((reconstructed - expected)/expected), Decimal('1e-55'))

    def test_center_gravity_zero(self):
        row = self.result['radial_samples'][0]
        for name in ('core_weighted','shell_weighted'):
            self.assertEqual(row[name]['total_g_proxy_m_s2'], Decimal(0))

    def test_center_energy_enclosed_zero(self):
        self.assertEqual(self.result['radial_samples'][0]['core_weighted']['enclosed_elastic_energy_J'], 0)

    def test_enclosed_fraction_inner_bounds(self):
        self.assertEqual(ged.enclosed_fraction(20,50), Decimal('0.064'))

    def test_enclosed_fraction_shell_before_core(self):
        self.assertEqual(ged.enclosed_fraction(20,100,50), Decimal(0))

    def test_enclosed_fraction_shell_outer(self):
        self.assertEqual(ged.enclosed_fraction(100,100,50), Decimal(1))

    def test_shell_fraction_75m_analytic(self):
        # 75^3 - 50^3 divided by 100^3 - 50^3 = 19/56.
        with localcontext() as ctx:
            ctx.prec = 70
            self.assertEqual(ged.enclosed_fraction(75,100,50), +(Decimal(19)/Decimal(56)))

    def test_core_boundary_difference_is_nonzero(self):
        row = next(item for item in self.result['radial_samples'] if item['radius_m']==50)
        self.assertGreater(row['difference_delta_g_proxy_m_s2'], 0)

    def test_core_boundary_independent_formula(self):
        # Difference = G*(.75-.25)*U/(c^2*r^2) at core edge.
        with localcontext() as ctx:
            ctx.prec = 70
            expected = Decimal('6.67430e-11') * Decimal('5e7') / (Decimal(299792458)**2 * Decimal(50)**2)
            row = next(item for item in self.result['radial_samples'] if item['radius_m']==50)
            self.assertLess(abs((row['difference_delta_g_proxy_m_s2']-expected)/expected), Decimal('1e-65'))

    def test_internal_field_difference_at_75m(self):
        row = next(item for item in self.result['radial_samples'] if item['radius_m']==75)
        self.assertGreater(row['difference_delta_g_proxy_m_s2'], 0)

    def test_exterior_delta_fields_equal(self):
        for row in self.result['radial_samples']:
            if row['radius_m'] >= self.params['outer_radius_m']:
                self.assertEqual(row['difference_delta_g_proxy_m_s2'], Decimal(0))

    def test_exterior_full_total_field_equal(self):
        for row in self.result['radial_samples']:
            if row['radius_m'] >= self.params['outer_radius_m']:
                self.assertEqual(row['core_weighted']['total_g_proxy_m_s2'],
                                 row['shell_weighted']['total_g_proxy_m_s2'])

    def test_exterior_monopole_independent_formula(self):
        row = self.result['radial_samples'][-1]
        with localcontext() as ctx:
            ctx.prec = 70
            expected = Decimal('6.67430e-11') * Decimal('1e8') / (Decimal(299792458)**2 * Decimal('1000')**2)
            actual = row['core_weighted']['elastic_delta_g_proxy_m_s2']
            self.assertLess(abs((actual-expected)/expected),Decimal('1e-65'))

    def test_exterior_inverse_square(self):
        a = next(row for row in self.result['radial_samples'] if row['radius_m']==200)
        b = next(row for row in self.result['radial_samples'] if row['radius_m']==1000)
        self.assertEqual(a['core_weighted']['elastic_delta_g_proxy_m_s2'] /
                         b['core_weighted']['elastic_delta_g_proxy_m_s2'], Decimal('25'))

    def test_baseline_same_at_all_radii(self):
        for row in self.result['radial_samples']:
            self.assertEqual(row['core_weighted']['baseline_g_m_s2'],
                             row['shell_weighted']['baseline_g_m_s2'])

    def test_core_loading_has_larger_inner_correction(self):
        row = next(row for row in self.result['radial_samples'] if row['radius_m']==25)
        self.assertGreater(row['core_weighted']['elastic_delta_g_proxy_m_s2'],
                           row['shell_weighted']['elastic_delta_g_proxy_m_s2'])

    def test_mass_equivalence_fixed_at_exterior(self):
        row = self.result['radial_samples'][-1]
        with localcontext() as ctx:
            ctx.prec = 70
            expected = Decimal('1e8') / Decimal(299792458)**2
            self.assertLess(abs((row['core_weighted']['mass_equivalent_enclosed_kg']-expected)/expected), Decimal('1e-65'))

    def test_zero_energy_enclosed_inside_shell_center(self):
        self.assertEqual(ged.enclosed_fraction(50,100,50), 0)

    def test_sample_radii_sorted(self):
        self.assertEqual(self.radii, sorted(self.radii))

    def test_serialization_preserves_small_values(self):
        v = json.loads(json.dumps(ged.json_safe(self.result)))
        x = v['radial_samples'][2]['difference_delta_g_proxy_m_s2']
        self.assertIsInstance(x, str)
        self.assertNotEqual(Decimal(x), 0)

    def test_invalid_negative_radius(self):
        with self.assertRaises(ValueError): ged.enclosed_fraction(-1,100,50)

    def test_invalid_equal_shell_bounds(self):
        with self.assertRaises(ValueError): ged.enclosed_fraction(10,50,50)

    def test_invalid_nonfinite(self):
        with self.assertRaises(ValueError): ged.D('Infinity')

    def test_invalid_config_id(self):
        cfg=copy.deepcopy(self.cfg); cfg['experiment_id']='WRONG'
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_core_radius(self):
        cfg=copy.deepcopy(self.cfg); cfg['core_radius_m']='100'
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_wrong_R5_pin(self):
        cfg=copy.deepcopy(self.cfg); cfg['r5_dependency_sha256']='0'*64
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_duplicate_sample(self):
        cfg=copy.deepcopy(self.cfg); cfg['sample_radii_m'] += ['100']
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_missing_outer_boundary(self):
        cfg=copy.deepcopy(self.cfg); cfg['sample_radii_m'].remove('100')
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_bad_fraction(self):
        cfg=copy.deepcopy(self.cfg); cfg['configurations']['core_weighted']['core_energy_fraction']='1'
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_reversed_fractions(self):
        cfg=copy.deepcopy(self.cfg); cfg['configurations']['core_weighted']['core_energy_fraction']='0.1'
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_negative_modulus(self):
        cfg=copy.deepcopy(self.cfg); cfg['materials']['core']['young_modulus_Pa']='-1'
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_invalid_negative_energy(self):
        cfg=copy.deepcopy(self.cfg); cfg['total_elastic_energy_J']='-1'
        with self.assertRaises(ValueError): ged.validate(cfg)

    def test_pinned_dependency_failure_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'ged_v1_r2.py').write_text('not valid source', encoding='utf-8')
            (root/'ged_v1_r5.py').write_text('not valid source', encoding='utf-8')
            with self.assertRaises(ValueError): ged.check_dependencies(root)

    def test_source_manifest_files_exist(self):
        for name in ged.INPUT_FILES:
            self.assertTrue((ged.ROOT/name).exists(), name)

    def test_scientific_limits_are_present(self):
        self.assertGreaterEqual(len(self.result['limitations']), 5)
        self.assertTrue(any('tensor' in item.lower() for item in self.result['limitations']))

    def test_standalone_clean_directory_run_and_receipts(self):
        # WHAT: Reconstruct a clean working tree from just the official files.
        # WHY: Verify import independence, frozen tests and manifest coverage.
        if os.environ.get('GED_R6_SNAPSHOT_TESTS') == '1':
            self.skipTest('Avoid nested clean-run invocation inside a frozen snapshot')
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for name in ged.INPUT_FILES:
                shutil.copy2(ged.ROOT/name, root/name)
            out=ged.run_experiment(root)
            self.assertTrue((out/'provenance.json').exists())
            self.assertTrue((out/'tests.txt').exists())
            self.assertIn('OK', (out/'tests.txt').read_text())
            manifest=(out/'SHA256SUMS').read_text().splitlines()
            self.assertEqual(len(manifest),len(ged.INPUT_FILES)+5)
            for entry in manifest:
                digest, relative = entry.split('  ',1)
                self.assertEqual(hashlib.sha256((out/relative).read_bytes()).hexdigest(),digest)
            self.assertEqual(json.loads((out/'provenance.json').read_text())['test_exit_status'],0)


if __name__ == '__main__':
    unittest.main()
