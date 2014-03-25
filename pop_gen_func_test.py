#!/usr/bin/env python

import unittest
from scipy.special import hyp2f1, hyp1f1
import pop_gen_func

class TestULinkage(unittest.TestCase):
    def runTest(self):
        self.assertGreater(pop_gen_func.uLinkage(2, 0.025, 2), 0.0)
        self.assertGreater(pop_gen_func.uLinkage(0, 0.025, 2), 0.0)
        self.assertGreater(pop_gen_func.uLinkage(-2.0, 0.025, 2), 0.0)

class TestHyp2F2(unittest.TestCase):
    def runTest(self):
        self.assertAlmostEqual(hyp1f1(1,2,3), 6.361845641062556)
        self.assertAlmostEqual(hyp2f1(1,2,3,-4.), 0.2988202609457375)
        self.assertAlmostEqual(hyp2f1(1,2,3,-0.5), 0.7562791351346849)
        self.assertAlmostEqual(hyp2f1(1,2,3,0.7), 2.0570318543915747)
        self.assertAlmostEqual(hyp2f1(1, 0.3, 1.3, -3.0), 0.7113010112875268)

class TestU0(unittest.TestCase):
    def runTest(self):
        self.assertGreater(pop_gen_func.u0(2, 0.025), 0.0)
        self.assertGreater(pop_gen_func.u0(0, 0.025), 0.0)
        self.assertGreater(pop_gen_func.u0(-2, 0.025), 0.0)

class TestMGeneral(unittest.TestCase):
    def runTest(self):
        a = 0
        mu = 0.025
        s = 2.0
        V = 1.0
        m = 5
    
        th = [
            [0.9457728047394913,0.02352668668506197,0.011982536075621617,
             0.0080213671250029,0.005991268037810812,0.004705337337012394],
            [0.9523205497232898,0.02004874548100449,0.008706643282440063,
             0.004952065583054215,0.0031305893984905973,0.010841406531721899],
            [0.9366219761117129,0.02198445658489576,0.010782199884747266,
             0.0069989917350783614,0.00510646538848006,0.01850591029508657]]
    
        for k in xrange(m + 1):
            self.assertAlmostEqual(pop_gen_func.mGeneral(a, mu, 0, 0, k, m), th[0][k])
            self.assertAlmostEqual(pop_gen_func.mGeneral(a, mu, 0, V, k, m), th[1][k])
            self.assertAlmostEqual(pop_gen_func.mGeneral(a, mu, s, V, k, m), th[2][k])

class TestQ1vsM(unittest.TestCase):
    def runTest(self):
        a = 0
        mu = 0.025
        s = 2.0
        V = 1.0
        m = 5
        t = 4.0
    
        th = [
            [0.7756278272629531,0.020356283378045733,0.011266356816185446,
             0.008737546384439078,0.009161671344827057,0.17485031481355084],
            [0.782101141273109,0.01698994323216128,0.008027815716697952,
             0.005630893148796326,0.006189391647333807,0.181060814981903],
            [0.8232228035723034,0.015933218904722574,0.006787423253861372,
             0.004412407973178893,0.004737710155317756,0.1449064361406177]]
    
        for k in xrange(m + 1):
            self.assertAlmostEqual(pop_gen_func.Q1vsM(k, m, t, mu, 0, 0), th[0][k])
            self.assertAlmostEqual(pop_gen_func.Q1vsM(k, m, t, mu, 0, V), th[1][k])
            self.assertAlmostEqual(pop_gen_func.Q1vsM(k, m, t, mu, s, V), th[2][k])

if __name__ == "__main__":
    unittest.main()