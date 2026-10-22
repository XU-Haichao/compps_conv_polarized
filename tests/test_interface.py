import unittest
import numpy as np
from polarized import blackbody_bins, build, PolarizedComppsc

class InterfaceTests(unittest.TestCase):
    def test_seed(self):
        edges = np.geomspace(1e-8,1e2,1001)
        self.assertAlmostEqual(blackbody_bins(edges.tolist()).sum(),1.,places=10)
        with self.assertRaises(ValueError):
            blackbody_bins([1,0])

    def test_native_smoke(self):
        edges = np.geomspace(1e-8,1e3,1101)
        model = PolarizedComppsc(build("build/test",(21,40,3)),edges,blackbody_bins(edges),20.,1.)
        values,audit = model.spectrum(.5)
        self.assertLessEqual(audit["last_relative_source_increment"],1e-7)
        for i,q in ((values[0],values[1]),(values[2],values[3])):
            self.assertTrue(np.isfinite(i).all() and np.isfinite(q).all())
            self.assertGreater(i.sum(),0)
            self.assertLessEqual(np.max(abs(q)-i),1e-7*np.max(i))

if __name__ == "__main__":
    unittest.main()
