import math
from scipy.special import gammaln, hyp1f1, binom, hyp2f1

EPS = 1e-6

def u0(sigma, mu):
    ret = mu
    if math.fabs(sigma) > EPS:
        ret = mu * sigma / (1.0 - math.exp(-sigma))
    return ret

def uLinkage(sigma, mu, V):
    if V < EPS:
        return u0(sigma, mu)
    ret = mu
    if sigma > EPS:
        val = hyp2f1(1.0, V / sigma, 1.0 + V / sigma, 1.0 - u0(sigma, mu) / mu)
        ret = u0(sigma, mu) * val;
    elif sigma < -EPS:
        ret = mu / (math.fabs(sigma) + V) * (u0(sigma, mu) / mu * math.fabs(sigma) + V);
    return ret;

# a is the current_allele {0,1}, e is the fitter allele {0,1}
def lambda_(a, e, uP, uM, gamma):
    return ((uM if a != e else uP) + gamma) / (uP + uM + 2.0 * gamma)

def zNeutral(mu):
    return math.exp(2 * gammaln(mu) - gammaln(2.0 * mu))

def pNeutral(a, mu, x):
    ret = 2.0 / zNeutral(mu)
    ret *= ((1.0 - a) * (1.0 - x) + a * x) * (x * (1.0 - x))**(-1.0 + mu)
    return ret

def mNeutral(a, mu, k, m):
    ret = binom(m, k)
    ret /= zNeutral(mu)
    ret *= (2.0 * (1.0 - a) * (m - k) + 2.0 * a * k + 2.0 * mu) * math.exp(gammaln(k + mu) + gammaln(-k + m + mu) - gammaln(1.0 + m + 2.0 * mu))
    return ret

def hitchhikingFraction(a, b, mu, V):
    hyperg1 = hyp1f1(b + mu, 1.0 + 2.0 * mu, (2.0 * a - 1.0) * V)
    hyperg2 = hyp1f1(b + mu, 2.0 + 2.0 * mu, (2.0 * a - 1.0) * V)
    ret = -a * math.exp(-V) + 1.0 - a
    ret /= zNeutral(mu)
    ret *= math.exp(gammaln(b + mu) + gammaln(1.0 - b + mu) - gammaln(1.0 + 2.0 * mu))
    term = -2.0 * a * hyperg1
    term += 2.0 * (1.0 - b + mu) / (1.0 + 2.0 * mu) * hyperg2
    ret *= term
    diff = mNeutral(a, mu, b, 1) - ret
    return diff

def pSweep(a, mu, V, x):
    return pNeutral(a, mu, x) * math.exp(-V * ((1.0 - a) * x + a * (1.0 - x)))

def mSweep(a, mu, V, k, m):
    ret = binom(m, k)
    hyperg1 = hyp1f1(k + mu, m + 2.0 * mu, (2.0 * a - 1.0) * V)
    hyperg2 = hyp1f1(k + mu, 1.0 + m + 2.0 * mu, (2.0 * a - 1.0) * V)
    ret /= zNeutral(mu)
    ret *= (-a * math.exp(-V) + (1.0 - a)) * math.exp(gammaln(k + mu) + gammaln(-k + m + mu) - gammaln(m + 2.0 * mu)) * ((-2.0 * a) * hyperg1 + 2.0 * (-k + m + mu) / (m + 2.0 * mu) * hyperg2)
    ret += hitchhikingFraction(a, 1, mu, V) if k == m else 0.0
    ret += hitchhikingFraction(a, 0, mu, V) if k == 0 else 0.0
    
    return ret;

def zSel(a, mu, s):
    hyperg = hyp1f1(mu, 2.0 * mu, s)
    ret = math.exp(2.0 * gammaln(mu) - gammaln(2.0 * mu)) * (1.0 - math.exp(-(1.0 - a) * s) * hyperg)
    return ret

def pSel(a, mu, s, x):
    ret = (x * (1.0 - x))**(-1.0 + mu) / zSel(a, mu, s)
    ret *= (1.0 - math.exp(s * ((1.0 - a) * (x - 1.0) + a * x)))
    return ret
    
def mSel(a, mu, s, k, m):
    hyperg = hyp1f1(k + mu, m + 2.0 * mu, s)
    ret = binom(m, k)
    ret /= zSel(a, mu, s)
    ret *= math.exp(gammaln(k + mu) + gammaln(m - k + mu) - gammaln(m + 2.0 * mu))
    ret *= (1.0 - math.exp(-(1.0 - a) * s) * hyperg)
    return ret

def hitchhikingFractionSel(a, b, mu, s, V):
    hyperg1 = hyp1f1(b + mu, 1.0 + 2.0 * mu, -(1.0 - a) * V + a * V)
    hyperg2 = hyp1f1(b + mu, 1.0 + 2.0 * mu, s - (1.0 - a) * V + a * V)
    ret = math.exp(gammaln(b + mu) + gammaln(1.0 - b + mu) - gammaln(1.0 + 2.0 * mu))
    ret *= math.exp(-(1.0 - a) * s - a * V)
    ret /= zSel(a, mu, s)
    ret *= math.exp((1.0 - a) * s) * hyperg1 - hyperg2
    diff = mSel(a, mu, s, b, 1)
    diff -= ret;
    return diff

def pSweepSel(a, mu, s, V, x):
    return pSel(a, mu, s, x) * math.exp(-V * ((1.0 - a) * x + a * (1.0 - x)))

def mSweepSel(a, mu, s, V, k, m):
    hyperg1 = hyp1f1(k + mu, m + 2.0 * mu, -(1.0 - a) * V + a * V)
    hyperg2 = hyp1f1(k + mu, m + 2.0 * mu, s - (1.0 - a) * V + a * V)
    ret = math.exp(gammaln(k + mu) + gammaln(m - k + mu) - gammaln(m + 2.0 * mu))
    ret *= binom(m, k) * math.exp(-(1.0 - a) * s - a * V) / zSel(a, mu, s)
    ret *= math.exp((1.0 - a) * s) * hyperg1 - hyperg2;
    ret += hitchhikingFractionSel(a, 1, mu, s, V) if k == m else 0.0
    ret += hitchhikingFractionSel(a, 0, mu, s, V) if k == 0 else 0.0
    return ret

def mGeneral(a, mu, s, V, k, m):
    if math.fabs(s) < EPS and V < EPS:
        ret = mNeutral(a, mu, k, m)
    elif math.fabs(s) < EPS:
        ret = mSweep(a, mu, V, k, m)
    elif math.fabs(V) < EPS:
        ret = mSel(a, mu, s, k, m)
    else:
        ret = mSweepSel(a, mu, s, V, k, m)
    return ret

def prop(aTo, aFrom, uP, uM, t):
    if aTo == 0 and aFrom == 0:
        return (uM + math.exp(-t * (uP + uM)) * uP) / (uP + uM)
    if aTo == 0 and aFrom == 1:
        return  (uM - math.exp(-t * (uP + uM)) * uM) / (uP + uM)
    if aTo == 1 and aFrom == 0:
        return (uP - math.exp(-t * (uP + uM)) * uP) / (uP + uM)
    if aTo == 1 and aFrom == 1:
        return  (uP + math.exp(-t * (uP + uM)) * uM) / (uP + uM)

def g(a1, a2, uP, uM, t):
    ret = 0.0
    for a in [0, 1]:
        ret += lambda_(a, 1, uP, uM, 0.0) * prop(a1, a, uP, uM, t) * prop(a2, a, uP, uM, t)
    return ret

def Q(k1, m1, k2, m2, t, mu, s, V, uP=None, uM=None):
    if uP is None:
        uP = uLinkage(s, mu, V)
    if uM is None:
        uM = uLinkage(-s, mu, V)
    ret = 0.0
    for a1 in [0, 1]:
        v1 = mGeneral(a1, mu, s, V, k1, m1)
        for a2 in [0, 1]:
            val = g(a1, a2, uP, uM, t) * v1 * mGeneral(a2, mu, s, V, k2, m2)
            ret += val
    return ret

def Q1vsM(k, m, t, mu, s, V):
    return Q(0, 1, k, m, t, mu, s, V) + Q(1, 1, m - k, m, t, mu, s, V)

def Qvec(k1, m1, m2, t, mu, s, V):
    ret = [0.0] * (m2 + 1)
    uP = uLinkage(s, mu, V)
    uM = uLinkage(-s, mu, V)
    for k in xrange(0, m2 + 1):
        ret[k] = Q(k1, m1, k, m2, t, mu, s, V, uP, uM)
    return ret

def Q1vsMvec(m, t, mu, s, V):
    q1 = Qvec(0, 1, m, t, mu, s, V)
    q2 = Qvec(1, 1, m, t, mu, s, V)
    for k in xrange(0, m + 1):
        q1[k] += q2[m - k]
    return q1

class IllegalParametersError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SingleSpectrumProb:
    def __init__(self, m):
        self.m = m
        self.lastP = (0.0, 0.0, 0.0, 0.0)
        self.last_prob = [0.0] * (m + 1)
  
    def __call__(t, mu, s, V):
        if self.lastP != (t, mu, s, V):
            self.last_prob = self.prob(t, mu, s, V)
            self.lastP = (t, mu, s, gamma, V)
        return last_prob
    
    def prob(t, mu, s, V):
        if t < 0:
            raise IllegalParametersError("t negative")
        if mu < 0:
            raise IllegalParametersError("mu negative")
        if s < 0:
            raise IllegalParametersError("s negative")
        if gamma < 0:
            raise IllegalParametersError("gamma negative")
        if V < 0:
            raise IllegalParametersError("V negative")
        return Q1vsM(m, t, mu, s, V)
