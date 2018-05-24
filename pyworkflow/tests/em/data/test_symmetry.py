# **************************************************************************
# *
# * Authors:     roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
#
from pyworkflow.tests import *
from pyworkflow.em.constants import SYM_CYCLIC, SYM_DIHEDRAL,\
    SYM_OCTAHEDRAL, SYM_TETRAHEDRAL, SYM_TETRAHEDRAL_Z3, SYM_I222,\
    SYM_I222r, SYM_In25, SYM_In25r
from pyworkflow.em.symmetry import getSymmetryMatrices
from pyworkflow.em.transformations import identity_matrix
import numpy


class TestSymmetry(unittest.TestCase):

    def assertArrayAlmostEqual(self, a1, a2):
        try:
            numpy.testing.assert_array_almost_equal(a1, a2, decimal=3)
            res = True
        except AssertionError as err:
            res = False
            print (err)
        self.assertTrue(res)

    def testSymmetryCyclic(self):
        matrices = getSymmetryMatrices(SYM_CYCLIC, n=3)
        m1 = matrices[0]
        m2 = matrices[1]
        m3 = matrices[2]

        #m1 is identity
        r = identity_matrix()
        self.assertArrayAlmostEqual(r, m1)

        #m2 rotated 120 degrees
        c = -0.5
        s = 0.86602540378
        r = [[c, -s, 0, 0],
             [s, c, 0, 0],
             [0, 0, 1.0, 0],
             [0, 0, 0.0, 1.0]]

        self.assertArrayAlmostEqual(r, m2)
        #m3 rotated -120 degrees
        c = -0.5
        s = -0.86602540378
        r = [[c, -s, 0, 0],
             [s, c, 0, 0],
             [0, 0, 1.0, 0],
             [0, 0, 0.0, 1.0]]

        self.assertArrayAlmostEqual(r, m3)

    def testSymmetryDihedral(self):
        matrices = getSymmetryMatrices(SYM_DIHEDRAL, n=3)
        # skip m1, m2 and m3 since there are identical to CYCLIC
        m4 = matrices[3]
        m5 = matrices[4]
        m6 = matrices[5]

        #m4 mirrored
        r = identity_matrix()
        r[1][1] *= -1.
        r[2][2] *= -1.
        self.assertArrayAlmostEqual(r, m4)

        #m5 rotated 120 degrees and mirrored
        c = -0.5
        s = 0.86602540378
        r = [[c, -s, 0, 0],
             [-s, -c, 0, 0],
             [0, 0, -1.0, 0],
             [0, 0, 0.0, 1.0]]
        self.assertArrayAlmostEqual(r, m5)

        # m6 rotated -120 degrees and mirrored
        c = -0.5
        s = -0.86602540378
        r = [[c, -s, 0, 0],
             [-s, -c, 0, 0],
             [0, 0, -1.0, 0],
             [0, 0, 0.0, 1.0]]
        self.assertArrayAlmostEqual(r, m6)

    def testSymmetryOctahedral(self):
        matrices = getSymmetryMatrices(SYM_OCTAHEDRAL)
        refMatrices=[
            [[ 1.,  0.,  0., 0.],  # 0
             [ 0.,  1.,  0., 0.],
             [ 0.,  0.,  1., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0., -1.,  0., 0.],  # 1
             [ 1.,  0.,  0., 0.],
             [ 0.,  0.,  1., 0.],
             [ 0.,  0.,  0., 1.]],
            [[-1.,  0.,  0., 0.],  # 2
             [ 0., -1.,  0., 0.],
             [ 0.,  0.,  1., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  1.,  0., 0.],  # 3
             [-1.,  0.,  0., 0.],
             [ 0.,  0.,  1., 0.],
             [ 0.,  0.,  0., 1.]],

            [[ 1.,  0.,  0., 0.],  # 4
             [ 0.,  0., -1., 0.],
             [ 0.,  1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0., -1.,  0., 0.],  # 5
             [ 0.,  0., -1., 0.],
             [ 1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[-1.,  0.,  0., 0.],  # 6
             [ 0.,  0., -1., 0.],
             [ 0., -1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  1.,  0., 0.],  # 7
             [ 0.,  0., -1., 0.],
             [-1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],

            [[  1.,  0.,  0., 0.],  # 8
             [  0., -1.,  0., 0.],
             [  0.,  0., -1., 0.],
             [  0.,  0.,  0., 1.]],
            [[  0., -1.,  0., 0.],  # 9
             [ -1.,  0.,  0., 0.],
             [  0.,  0., -1., 0.],
             [  0.,  0.,  0., 1.]],
            [[-1.,  0.,  0., 0.],  # 10
             [ 0.,  1.,  0., 0.],
             [ 0.,  0., -1., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  1.,  0., 0.],  # 11
             [ 1.,  0.,  0., 0.],
             [ 0.,  0., -1., 0.],
             [ 0.,  0.,  0., 1.]],


            [[ 1.,  0.,  0., 0.],  # 12
             [ 0.,  0.,  1., 0.],
             [ 0., -1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0., -1.,  0., 0.],  # 13
             [ 0.,  0.,  1., 0.],
             [-1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[-1.,  0.,  0., 0.],  # 14
             [ 0.,  0.,  1., 0.],
             [ 0.,  1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  1.,  0., 0.],  # 15
             [ 0.,  0.,  1., 0.],
             [ 1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],

            [[ 0.,  0.,  1., 0.],  # 16
             [ 0.,  1.,  0., 0.],
             [-1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  0.,  1., 0.],  # 17
             [ 1.,  0.,  0., 0.],
             [ 0.,  1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  0.,  1., 0.],  # 18
             [ 0., -1.,  0., 0.],
             [ 1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  0.,  1., 0.],  # 19
             [-1.,  0.,  0., 0.],
             [ 0., -1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],

            [[  0.,  0., -1., 0.],  # 20
             [  0.,  1.,  0., 0.],
             [  1.,  0.,  0., 0.],
             [  0.,  0.,  0., 1.]],
            [[  0.,  0., -1., 0.],  # 21
             [  1.,  0.,  0., 0.],
             [  0., -1.,  0., 0.],
             [  0.,  0.,  0., 1.]],
            [[ 0.,  0., -1., 0.],  # 22
             [ 0., -1.,  0., 0.],
             [-1.,  0.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
            [[ 0.,  0., -1., 0.],  # 23
             [-1.,  0.,  0., 0.],
             [ 0.,  1.,  0., 0.],
             [ 0.,  0.,  0., 1.]],
       ]

        for m1, m2 in zip(matrices[:len(refMatrices)],refMatrices):
            self.assertArrayAlmostEqual(m1, m2)

    def testSymmetryTetrahedral(self):
        matrices = getSymmetryMatrices(SYM_TETRAHEDRAL)
        refMatrices= [
            [[ 1.00000000e+00,  6.40987562e-17, -4.53246652e-17, 0.00000000e+00],
             [ 5.55111512e-17,  1.00000000e+00, -5.55111512e-17, 0.00000000e+00],
             [ 0.00000000e+00, -5.55111512e-17,  1.00000000e+00, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[ 2.22044605e-16,  5.77350269e-01,  8.16496581e-01, 0.00000000e+00],
             [ 5.77350269e-01, -6.66666667e-01,  4.71404521e-01, 0.00000000e+00],
             [ 8.16496581e-01,  4.71404521e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[-2.22044605e-16, -5.77350269e-01, -8.16496581e-01, 0.00000000e+00],
             [-5.77350269e-01, -6.66666667e-01,  4.71404521e-01, 0.00000000e+00],
             [-8.16496581e-01,  4.71404521e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[-1.00000000e+00,  8.58760498e-18, -1.56346968e-16, 0.00000000e+00],
             [ 5.55111512e-17,  3.33333333e-01, -9.42809042e-01, 0.00000000e+00],
             [ 1.66533454e-16, -9.42809042e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[-5.00000000e-01, -8.66025404e-01,  1.21208789e-16, 0.00000000e+00],
             [ 8.66025404e-01, -5.00000000e-01, -5.55111512e-17, 0.00000000e+00],
             [ 0.00000000e+00,  1.38777878e-16,  1.00000000e+00, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[-5.00000000e-01,  8.66025404e-01,  1.01864860e-17, 0.00000000e+00],
             [-8.66025404e-01, -5.00000000e-01,  5.55111512e-17, 0.00000000e+00],
             [ 1.66533454e-16,  0.00000000e+00,  1.00000000e+00, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[-0.5,             0.28867513,     -0.81649658,     0.],
             [-0.28867513,      0.83333333,      0.47140452,     0.],
             [ 0.81649658,      0.47140452,     -0.33333333,     0.],
             [ 0.,              0.,              0.,             1.]],
            [[-0.5,            -0.28867513,      0.81649658,     0.],
             [ 0.28867513,      0.83333333,      0.47140452,     0.],
             [-0.81649658,      0.47140452,     -0.33333333,     0.],
             [ 0.,              0.,              0.,             1.]],
            [[ 5.00000000e-01, -2.88675135e-01,  8.16496581e-01, 0.00000000e+00],
             [-8.66025404e-01, -1.66666667e-01,  4.71404521e-01, 0.00000000e+00],
             [-5.55111512e-17, -9.42809042e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[ 5.00000000e-01, -8.66025404e-01,  2.32231091e-16, 0.00000000e+00],
             [-2.88675135e-01, -1.66666667e-01, -9.42809042e-01, 0.00000000e+00],
             [ 8.16496581e-01,  4.71404521e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[ 5.00000000e-01,  8.66025404e-01, -1.56346968e-16, 0.00000000e+00],
             [ 2.88675135e-01, -1.66666667e-01, -9.42809042e-01, 0.00000000e+00],
             [-8.16496581e-01,  4.71404521e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
            [[ 5.00000000e-01,  2.88675135e-01, -8.16496581e-01, 0.00000000e+00],
             [ 8.66025404e-01, -1.66666667e-01,  4.71404521e-01, 0.00000000e+00],
             [-5.55111512e-16, -9.42809042e-01, -3.33333333e-01, 0.00000000e+00],
             [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]]
        ]
        for m1, m2 in zip(matrices,refMatrices):
            self.assertArrayAlmostEqual(m1, m2)

    @classmethod
    def tearDownClass(cls):
        pass

    @classmethod
    def setUpClass(cls):
        pass