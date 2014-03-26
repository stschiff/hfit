#!/usr/bin/env python

import argparse
import string
import sys
import model
from scipy import optimize

models = ["driverfield_neutral", "driverfield_neutral_sel", "unlinked_neutral", "unlinked_neutral_sel", "constrained_non_neutral"]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", "-m", default="driverfield_neutral")
    parser.add_argument("--mu", type=float)
    parser.add_argument("--V", type=float)
    parser.add_argument("--t", type=float)
    parser.add_argument("--spectrum_file")
    parser.add_argument("input_filename")
    args = parser.parse_args()
    
    run(args)

def run(args):
    spectrum = read_spectrum(args.input_filename)
    score_func = get_score_func(args, spectrum)
    
    x_min = get_minimum(score_func)
    
    write_param_header()
    write_params(score_func, x_min)

def read_spectrum(input_filename):
    input_file = sys.stdin if input_filename == "-" else open(input_filename, "r")
    spectrum = []
    for line in input_file:
        fields = string.split(string.strip(line))
        spectrum.append(int(fields[1]))
    sys.stderr.write("read spectrum: {}\n".format(spectrum))
    return spectrum

def get_score_func(args, spectrum):
    if args.model == "driverfield_neutral":
        return model.DriverfieldNeutralScore(spectrum)
    elif args.model == "driverfield_neutral_sel":
        return model.DriverfieldNeutralScoreSel(spectrum)
    elif args.model == "unlinked_neutral":
        return model.UnlinkedNeutralScore(spectrum)
    elif args.model == "unlinked_neutral_sel":
        return model.UnlinkedNeutralScoreSel(spectrum)
    elif args.model == "constrained_non_neutral":
        return model.ConstrainedNonNeutralScore(spectrum, args.mu, args.V, args.t)
    else:
        raise ValueError("unknown model name: ", args.model)

def get_minimum(score_func):
    x0 = score_func.initial_params()
    bounds = score_func.bounds()
    def cb(x):
        print "minimizing, current values", x
    res = optimize.minimize(score_func, x0, method='L-BFGS-B', bounds=bounds, callback=cb)
    
    if res.success:
        print "found result:", res.x
    else:
        print "problem occurred:", res

def write_param_header():
    print "mu\tV\tt\ts\tgamma\tc\tcw\tscore"

def write_params(score_func, x):
    p = score_func.make_params(x)
    score = score_func.get_score(p)
    print "\t".join(map(str, [p["mu"], p["V"], p["t"], p["s"], p["gamma"], p["c"], p["cw"], score]))

main()

# void writeParamVar(File file, SingleSpectrumScore scoreFunc, double[] x) {
#   auto p = ["mu":0.0, "V":0.0, "t":0.0, "s":0.0, "gamma":0.0, "c":0.0, "cw":0.0];
#   auto names = scoreFunc.paramNames();
#   foreach(i, n; names)
#     p[n] = x[i];
#   file.writefln("%s\t%s\t%s\t%s\t%s\t%s\t%s", p["mu"], p["V"], p["t"], p["s"], p["gamma"], p["c"], p["cw"]);
# }
# 
# class Accumulator {
#   int nrParams;
#   double[] sums;
#   double[] sums_sq;
#   double norm;
#   
#   this(int nrParams) {
#     this.nrParams = nrParams;
#     sums = new double[nrParams];
#     sums[] = 0.0;
#     sums_sq = new double[nrParams];
#     sums_sq[] = 0.0;
#     norm = 0.0;
#   }
#   
#   void add(double[] x)
#   in {
#     assert(x.length == nrParams);
#   }
#   body {
#     foreach(i, xx; x) {
#       sums[i] += xx;
#       sums_sq[i] += xx * xx;
#     }
#     norm += 1.0;
#   }
# 
#   double[] mean() {
#     return iota(nrParams).map!(i => sums[i] / norm)().array;
#   }
#   
#   double[] stddev() {
#     return iota(nrParams).map!(i => sqrt(sums_sq[i] / norm - (sums[i] / norm)^^2))().array;
#   }
# }