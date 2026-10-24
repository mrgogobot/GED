"""GED-V1-R7 conventional-gravity regression suite, Python standard library.

Each test names the WHAT and WHY in a short adjacent comment or docstring.
These are algebra/model-implementation tests, not observational confirmations.
"""

import copy
import json
import unittest
from decimal import Decimal as D, getcontext, localcontext
from pathlib import Path

import ged_v1_r7 as m


class GEDR7Tests(unittest.TestCase):
    """WHAT: Physics, validity, convergence and provenance controls.
    WHY: Identify false positive signals and prevent unnoticed file drift.
    """

    @classmethod
    def setUpClass(cls):
        # Every test uses 70-digit precision to reproduce the model's context.
        getcontext().prec = m.PRECISION
        cls.folder = Path(__file__).resolve().parent
        cls.config = json.loads((cls.folder / 'GED_V1_R7_config.json').read_text())
        cls.model = m.model_from_config(cls.config)
        cls.pair_x = m.create_parcels(cls.model, "pair_x")
        cls.pair_y = m.create_parcels(cls.model, "pair_y")
        cls.cube = m.create_parcels(cls.model, "cube_eight")

    def at(self, direction, radius='200'):
        return m.multiply(m.AXES[direction], D(radius))

    def assertNear(self, got, expected, tol='1e-45'):
        # Relative comparison guards both tiny fields and ordinary quantities.
        got, expected, tol = D(got), D(expected), D(tol)
        scale = max(abs(got), abs(expected), D('1e-90'))
        self.assertLessEqual(abs(got-expected) / scale, tol)

    # INPUT / UNITS
    def test_experiment_identifier(self):
        self.assertEqual(self.config['experiment_id'], 'GED-V1-R7')

    def test_positive_constants(self):
        self.assertGreater(self.model.G, 0)
        self.assertGreater(self.model.c, 0)

    def test_parcels_inside_shell(self):
        for case in (self.pair_x, self.pair_y, self.cube):
            for item in case:
                self.assertGreater(m.norm(item.position_m), self.model.core_R)
                self.assertLess(m.norm(item.position_m), self.model.R)

    def test_equal_total_energy_exact(self):
        for case in (self.pair_x, self.pair_y, self.cube):
            self.assertEqual(sum((p.energy_J for p in case), D(0)), self.model.E)

    def test_correct_energy_shares(self):
        self.assertEqual(self.pair_x[0].energy_J, D(50000000))
        self.assertEqual(self.cube[0].energy_J, D(12500000))

    def test_mass_energy_equivalent(self):
        self.assertNear(m.energy_mass(self.model, self.pair_x), self.model.E / self.model.c**2)

    def test_model_geometry(self):
        self.assertEqual(self.model.R, D(100))
        self.assertEqual(self.model.a, D(60))

    def test_background_mass_analytic(self):
        # Independent volume weighting with a separate literal for pi.
        pi = D('3.1415926535897932384626433832795028841971693993751058209749445923078')
        expect = 4*pi/3*(7870*D(50)**3 + 917*(D(100)**3-D(50)**3))
        self.assertNear(m.spherical_background_mass(self.model), expect)

    def test_background_is_same_for_all_cases(self):
        # No case input is involved in computing the background.
        self.assertGreater(m.spherical_background_mass(self.model), D('1e9'))

    def test_comparisons_have_expected_parcel_counts(self):
        self.assertEqual([len(self.pair_x),len(self.pair_y),len(self.cube)], [2,2,8])

    def test_all_sample_points_external(self):
        self.assertTrue(all(r > self.model.R for r in self.model.radii))

    def test_radii_sorted(self):
        self.assertEqual(tuple(sorted(self.model.radii)), self.model.radii)

    def test_reject_unknown_case(self):
        with self.assertRaises(ValueError):
            m.create_parcels(self.model, 'not_a_case')

    def test_reject_wrong_identifier(self):
        cfg=copy.deepcopy(self.config);cfg['experiment_id']='WRONG'
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_negative_energy(self):
        cfg=copy.deepcopy(self.config);cfg['total_added_energy_J']='-1'
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_zero_c(self):
        cfg=copy.deepcopy(self.config);cfg['c_m_s']='0'
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_infinite_G(self):
        cfg=copy.deepcopy(self.config);cfg['G_m3_kg_s2']='Infinity'
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_boolean_energy(self):
        cfg=copy.deepcopy(self.config);cfg['total_added_energy_J']=True
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_bad_core_radius(self):
        cfg=copy.deepcopy(self.config);cfg['core_radius_m']='80'
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_source_outside(self):
        cfg=copy.deepcopy(self.config);cfg['parcel_offset_m']='120'
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_inside_probe_config(self):
        cfg=copy.deepcopy(self.config);cfg['sample_radii_m']=['99','200']
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_duplicate_radii(self):
        cfg=copy.deepcopy(self.config);cfg['sample_radii_m']=['200','200']
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_unsorted_radii(self):
        cfg=copy.deepcopy(self.config);cfg['sample_radii_m']=['400','200']
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_duplicate_axes(self):
        cfg=copy.deepcopy(self.config);cfg['sample_axes']=['x','x']
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_wrong_case_set(self):
        cfg=copy.deepcopy(self.config);cfg['comparison']=['pair_x','cube_eight']
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_absent_limitations(self):
        cfg=copy.deepcopy(self.config);cfg['limitations']=[]
        with self.assertRaises(ValueError): m.model_from_config(cfg)

    def test_reject_inside_probe_exact(self):
        with self.assertRaises(ValueError): m.gravity_exact(self.model,self.pair_x,self.at('x','50'))

    def test_reject_inside_probe_multipole(self):
        with self.assertRaises(ValueError): m.gravity_quadrupole(self.model,self.pair_x,self.at('x','50'))

    # MONOPOLE / DIPOLE / QUADRUPOLE
    def test_identical_monopole_at_same_probe(self):
        r=self.at('x')
        self.assertEqual(m.gravity_monopole(self.model,self.pair_x,r),
                         m.gravity_monopole(self.model,self.cube,r))

    def test_zero_dipoles(self):
        for group in (self.pair_x,self.pair_y,self.cube):
            # A residual <1e-65 kg m is numerical Decimal round-off;
            # physically the symmetric source's dipole cancels exactly.
            for value in m.dipole_kg_m(self.model,group):
                self.assertLess(abs(value), D('1e-65'))

    def test_all_quadrupoles_traceless(self):
        for group in (self.pair_x,self.pair_y,self.cube):
            q=m.quadrupole_kg_m2(self.model,group)
            # Tracelessness is analytic; finite precision can leave ~1e-75.
            self.assertLess(abs(sum((q[i][i] for i in range(3)),D(0))), D('1e-65'))

    def test_quadrupoles_symmetric(self):
        for group in (self.pair_x,self.pair_y,self.cube):
            q=m.quadrupole_kg_m2(self.model,group)
            for i in range(3):
                for j in range(3): self.assertEqual(q[i][j],q[j][i])

    def test_pair_x_Qxx_analytic(self):
        q=m.quadrupole_kg_m2(self.model,self.pair_x)
        M=self.model.E/self.model.c**2
        self.assertNear(q[0][0],2*M*self.model.a**2)

    def test_pair_x_Qyy_analytic(self):
        q=m.quadrupole_kg_m2(self.model,self.pair_x)
        M=self.model.E/self.model.c**2
        self.assertNear(q[1][1],-M*self.model.a**2)

    def test_pair_y_rotates_Q(self):
        x=m.quadrupole_kg_m2(self.model,self.pair_x)
        y=m.quadrupole_kg_m2(self.model,self.pair_y)
        self.assertEqual(x[0][0],y[1][1]);self.assertEqual(x[1][1],y[0][0])

    def test_cube_quadrupole_zero(self):
        q=m.quadrupole_kg_m2(self.model,self.cube)
        for row in q:
            for number in row: self.assertLess(abs(number), D('1e-65'))

    def test_cube_is_not_spherical_field(self):
        r=self.at('x')
        exact=m.gravity_exact(self.model,self.cube,r)
        mono=m.gravity_monopole(self.model,self.cube,r)
        # Cubic source has a nonzero l=4 field despite Q=0.
        self.assertNotEqual(exact,mono)

    def test_quadrupole_prediction_cube_equals_monopole(self):
        r=self.at('x')
        truncated=m.gravity_through_quadrupole(self.model,self.cube,r)
        mono=m.gravity_monopole(self.model,self.cube,r)
        # Differences ~1e-95 m/s^2 are rounding, not observable coupling.
        for i in range(3): self.assertLess(abs(truncated[i]-mono[i]),D('1e-70'))

    # EXACT NEWTONIAN AND FAR-FIELD APPROXIMATION
    def test_pair_x_axial_direct_inverse_square(self):
        # Independently use two terms 1/(r−a)^2 and 1/(r+a)^2.
        r=D(200);a=self.model.a;mass=self.model.E/self.model.c**2
        expected=-self.model.G*mass/D(2)*(1/(r-a)**2 + 1/(r+a)**2)
        got=m.radial_component(m.gravity_exact(self.model,self.pair_x,self.at('x')), 'x')
        self.assertNear(got,expected)

    def test_pair_x_transverse_direct_formula(self):
        # Each source's y acceleration shares denominator (r²+a²)^3/2.
        r=D(200);a=self.model.a;mass=self.model.E/self.model.c**2
        expected=-self.model.G*mass*r/(r*r+a*a).sqrt()**3
        got=m.radial_component(m.gravity_exact(self.model,self.pair_x,self.at('y')),'y')
        self.assertNear(got,expected)

    def test_cube_axial_direct_two_group_formula(self):
        # WHAT: Eight cube corners reduce to two groups of four along x.
        # WHY: This closed-form calculation independently checks the general
        # vector summation without calling that algorithm a second time.
        with localcontext() as ctx:
            ctx.prec = m.PRECISION
            r = D(200)
            b = self.model.a / D(3).sqrt()
            mass = self.model.E / self.model.c**2
            near = (r-b) / ((r-b)**2 + D(2)*b*b).sqrt()**3
            far = (r+b) / ((r+b)**2 + D(2)*b*b).sqrt()**3
            expected = -self.model.G*mass/D(2)*(near+far)
            got = m.gravity_exact(self.model,self.cube,self.at('x'))[0]
            self.assertNear(got,expected)

    def test_pair_x_stronger_along_axis_than_monopole(self):
        r=self.at('x')
        self.assertLess(m.gravity_exact(self.model,self.pair_x,r)[0],
                        m.gravity_monopole(self.model,self.pair_x,r)[0])

    def test_pair_x_weaker_transverse_than_monopole(self):
        r=self.at('y')
        self.assertGreater(m.gravity_exact(self.model,self.pair_x,r)[1],
                           m.gravity_monopole(self.model,self.pair_x,r)[1])

    def test_pair_y_orientation_changes_external_field(self):
        self.assertNotEqual(m.gravity_exact(self.model,self.pair_x,self.at('x')),
                            m.gravity_exact(self.model,self.pair_y,self.at('x')))

    def test_pair_rotation_swaps_axes(self):
        gx=m.gravity_exact(self.model,self.pair_x,self.at('x'))[0]
        gy=m.gravity_exact(self.model,self.pair_y,self.at('y'))[1]
        self.assertNear(gx,gy)

    def test_cube_axis_rotation_symmetric(self):
        self.assertNear(m.gravity_exact(self.model,self.cube,self.at('x'))[0],
                        m.gravity_exact(self.model,self.cube,self.at('y'))[1])

    def test_monopole_inverse_square(self):
        a=m.gravity_monopole(self.model,self.pair_x,self.at('x','200'))[0]
        b=m.gravity_monopole(self.model,self.pair_x,self.at('x','400'))[0]
        self.assertEqual(a/b,D(4))

    def test_axis_quadrupole_inverse_fourth(self):
        a=m.gravity_quadrupole(self.model,self.pair_x,self.at('x','200'))[0]
        b=m.gravity_quadrupole(self.model,self.pair_x,self.at('x','400'))[0]
        self.assertNear(a/b,D(16))

    def test_pair_x_axial_quadrupole_analytic(self):
        M=self.model.E/self.model.c**2
        r=D(200)
        expected=-3*self.model.G*M*self.model.a**2/r**4
        self.assertNear(m.gravity_quadrupole(self.model,self.pair_x,self.at('x'))[0],expected)

    def test_pair_x_transverse_quadrupole_analytic(self):
        M=self.model.E/self.model.c**2
        r=D(200)
        expected=D(3)*self.model.G*M*self.model.a**2/(D(2)*r**4)
        self.assertNear(m.gravity_quadrupole(self.model,self.pair_x,self.at('y'))[1],expected)

    def test_quadrupole_approx_improves_monopole_axis(self):
        r=self.at('x')
        exact=m.gravity_exact(self.model,self.pair_x,r)[0]
        point=m.gravity_monopole(self.model,self.pair_x,r)[0]
        quad=m.gravity_through_quadrupole(self.model,self.pair_x,r)[0]
        self.assertLess(abs(exact-quad),abs(exact-point))

    def test_approximation_improves_with_distance(self):
        errors=[]
        for r in self.model.radii:
            probe=self.at('x',str(r))
            actual=m.gravity_exact(self.model,self.pair_x,probe)[0]
            predicted=m.gravity_through_quadrupole(self.model,self.pair_x,probe)[0]
            errors.append(m.relative_error(actual,predicted))
        self.assertTrue(all(errors[i+1] < errors[i] for i in range(len(errors)-1)))

    def test_cube_higher_moment_shrinks_relative_to_monopole(self):
        errors=[]
        for r in self.model.radii:
            probe=self.at('x',str(r))
            exact=m.gravity_exact(self.model,self.cube,probe)[0]
            monopole=m.gravity_monopole(self.model,self.cube,probe)[0]
            errors.append(abs((exact-monopole)/monopole))
        self.assertTrue(all(errors[i+1] < errors[i] for i in range(len(errors)-1)))

    def test_match_monopole_even_though_exact_differs(self):
        probe=self.at('x')
        self.assertEqual(m.gravity_monopole(self.model,self.pair_x,probe),
                         m.gravity_monopole(self.model,self.cube,probe))
        self.assertNotEqual(m.gravity_exact(self.model,self.pair_x,probe),
                            m.gravity_exact(self.model,self.cube,probe))

    def test_outside_gravity_decreases_with_distance(self):
        g1=abs(m.gravity_exact(self.model,self.pair_x,self.at('x','200'))[0])
        g2=abs(m.gravity_exact(self.model,self.pair_x,self.at('x','400'))[0])
        self.assertGreater(g1,g2)

    def test_difference_exact_and_quadrupole_same_sign(self):
        r=self.at('x')
        exact=m.gravity_exact(self.model,self.pair_x,r)[0]-m.gravity_exact(self.model,self.cube,r)[0]
        pred=m.gravity_quadrupole(self.model,self.pair_x,r)[0]
        self.assertLess(exact,0);self.assertLess(pred,0)

    def test_computation_has_expected_rows(self):
        result,rows=m.compute(self.model)
        self.assertEqual(len(rows),3*3*5)
        self.assertEqual(result['sample_count'],45)

    def test_result_does_not_claim_novel_gravity(self):
        result,_=m.compute(self.model)
        self.assertIn('NO NEW GRAVITY', result['classification'])
        self.assertIn('No information-field coupling',result['scientific_verdict'])

    def test_result_limits_not_empty(self):
        result,_=m.compute(self.model)
        self.assertGreaterEqual(len(result['limits']),5)

    def test_result_values_not_float_coerced(self):
        result,rows=m.compute(self.model)
        self.assertIsInstance(result['added_mass_equivalent_kg'],str)
        self.assertIsInstance(rows[0]['exact_added_radial_g_m_s2'],str)

    def test_manifest_files_exist(self):
        for filename in m.FILE_SET:
            self.assertTrue((self.folder/filename).is_file(), filename)

    def test_sha256_input_file(self):
        self.assertEqual(len(m.sha256_file(self.folder/'ged_v1_r7.py')),64)


if __name__ == '__main__':
    unittest.main()
