import numpy as np
from pop_gen_func import SpectrumProb

PENALTY = 1.0e20;

class SpectrumScore(object):    
    def __init__(self, spectrum_data, nr_params):
        self.m = len(spectrum_data) - 1
        self.neutral_spectrum_prob_func = SpectrumProb(self.m)
        self.selected_spectrum_prob_func = SpectrumProb(self.m)
        self.norm = sum(spectrum_data)
        self.spectrum_data = np.array(spectrum_data)
        self.nr_params = nr_params
    
    def check_params(self, p):
        if p["c"] < 0.0 or p["c"] > 1.0:
            raise ValueError("c out of range")
        if p["cw"] < 0.0 or p["cw"] > 1.0:
            raise ValueError("cw out of range")
        if p["gamma"] < 0.0 or 2.0 * p["gamma"] * p["t"] > 1.0:
            raise ValueError("gamma out of range")
        if p["s"] < 0.0 or p["s"] > 150.0:
            raise ValueError("selection out of range")
    
    def get_single_spectrum_probs(self, p):
        self.check_params(p)
        ret = []
        if self.norm > 0:
            self.check_params(p)
            q0 = self.neutral_spectrum_prob_func(p["t"], p["mu"], 0.0, p["V"])
            qs = self.selected_spectrum_prob_func(p["t"], p["mu"], p["s"], p["V"])
            neutral_part = q0 * p["c"]
            weak_sel_part = qs * (1.0 - p["c"]) * p["cw"]
            strong_sel_part = np.zeros(self.m + 1)
            strong_sel_part[0] = 1.0 - 2.0 * p["gamma"] * p["t"]
            strong_sel_part[-1] = 2.0 * p["gamma"] * p["t"]
            strong_sel_part *= (1.0 - p["c"]) * (1.0 - p["cw"])
            return neutral_part + weak_sel_part + strong_sel_part
        else:
            return np.ones(self.m + 1)
        return ret
    
    def get_score(self, p):
        probs = self.get_single_spectrum_probs(p)
        print "probs:", probs
        return -np.sum(np.log(probs) * self.spectrum_data)
    
    def __call__(self, x):
        assert len(x) == self.nr_params
        p = self.make_spectrum_params(x)
        return self.get_score(p)

    def make_vec_from_params(self, p):
        keys = self.param_names()
        return [p[k] for k in keys]
    
    def param_names(self):
        raise NotImplementedError()
    
    def make_spectrum_params(self):
        raise NotImplementedError()
    
    def initial_params(self):
        raise NotImplementedError()
    
    def bounds(self):
        raise NotImplementedError()

class FullScore(SpectrumScore):
    def __init__(self, spectrum):
        super(FullScore, self).__init__(spectrum, 7)
    
    def make_spectrum_params(self, x):
        keys = self.param_names()
        return dict(zip(keys, x))
    
    def param_names(self):
        return ["mu", "t", "V", "gamma", "c", "s", "cw"]
    
    def initial_params(self):
        return [0.01, 0.1, 4, 0.0, 0.5, 1.0, 0.5]
    
    def bounds(self):
        return [(0.0, None), (0.0, None), (0.0, None), (0.0, None), (0.0, 1.0), (0.0, 150.0), (0.0, 1.0)]
    
class DriverfieldNeutralScoreSel(SpectrumScore):
    def __init__(self, spectrum):
        super(DriverfieldNeutralScoreSel, self).__init__(spectrum, 4)
    
    def make_spectrum_params(self, x):
        mu, V, t, s = x
        return {"mu":mu, "t":t, "V":V, "s":s, "gamma":0, "c":0.0, "cw":1.0}
    
    def param_names(self):
        return ["mu", "t", "V", "s"]
    
    def initial_params(self):
        return [0.01, 0.1, 4.0, 0.0]

    def bounds(self):
        return [(0.0, None), (0.0, None), (0.0, None), (0.0, 150.0)]

class DriverfieldNeutralScore(SpectrumScore):
    def __init__(self, spectrum):
        super(DriverfieldNeutralScore, self).__init__(spectrum, 3)
    
    def make_spectrum_params(self, x):
        mu, V, t = x
        return {"mu":mu, "t":t, "V":V, "s":0.0, "gamma":0.0, "c":1.0, "cw":0.0}
    
    def param_names(self):
        return ["mu", "t", "V"]
    
    def initial_params(self):
        return [0.01, 0.1, 4.0]

    def bounds(self):
        return [(0.0, None), (0.0, None), (0.0, None)]

class UnlinkedNeutralScoreSel(SpectrumScore):
    def __init__(self, spectrum):
        super(UnlinkedNeutralScoreSel, self).__init__(spectrum, 3)
  
    def make_spectrum_params(self, x):
        mu, t, s = x
        return {"mu":mu, "t":t, "s":s, "V":0.0, "gamma":0.0, "c":0.0, "cw":1.0}
    
    def param_names(self):
        return ["mu", "t", "s"]
    
    def initial_params(self):
        return [0.01, 4.0, 0.0]

    def bounds(self):
        return [(0.0, None), (0.0, None), (0.0, 150.0)]
    
class UnlinkedNeutralScore(SpectrumScore):
    def __init__(self, spectrum):
        super(UnlinkedNeutralScore, self).__init__(spectrum, 2)
  
    def make_spectrum_params(self, x):
        mu, t = x
        return {"mu":mu, "t":t, "s":0.0, "V":0.0, "gamma":0.0, "c":1.0, "cw":0.0}
  
    def param_names(self):
        return ["mu", "t"]
  
    def initial_params(self):
        return [0.01, 4]

    def bounds(self):
        return [(0.0, None), (0.0, None)]

class ConstrainedNonNeutralScore(SpectrumScore):
  
    def __init__(self, spectrum, mu, V, t):
        super(ConstrainedNonNeutralScore, self).__init__(spectrum, 4)
        self.params = {"mu":mu, "V":V, "t":t}
    
    def make_spectrum_params(self, x):
        gamma, c, s, cw = x
        self.params["gamma"] = gamma
        self.params["c"] = c
        self.params["s"] = s
        self.params["cw"] = cw
        return params
    
    def param_names(self):
        return ["gamma", "c", "s", "cw"]
    
    def initial_params(self):
        return [0.0, 0.5, 2.0, 0.5]

    def bounds(self):
        return [(0.0, None), (0.0, 1.0), (0.0, 150.0), (0.0, 1.0)]
    