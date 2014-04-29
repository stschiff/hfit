import std.stdio;
import std.algorithm;
import std.conv;
import std.getopt;
import std.exception;
import std.file;
import std.string;
import std.array;
import std.range;
import std.math;
import std.random;
import modelScores;
import popGenFunc;
import powell;

enum model_t {driverfieldNeutral, driverfieldNeutralSel, unlinkedNeutral, unlinkedNeutralSel, constrainedNonNeutral, constrainedSel}

string neutralFitFileName = "";
model_t model = model_t.driverfieldNeutral;
auto maxSteps = 500U;
double mu, V, t; // for constrained model
string spectrumFilename = "/dev/null";
string inputData;
ulong[] spectrum;
bool noBootstrap;
auto nrBootStrapSamples = 20U;

void main(string[] args) {
  
    gsl_set_error_handler_off();
    
    try {
        readCommandlineParams(args);
    }
    catch(Exception e) {
        stderr.writeln(e.msg);
        displayHelpMessage();
        return;
    }
    run();
}

void readCommandlineParams(string[] args) {
    getopt(args, std.getopt.config.passThrough, "model|m", &model, "maxSteps", &maxSteps, "mu", &mu, "V", &V, "t", &t, "spectrumFile", &spectrumFilename, "noBootstrap", &noBootstrap, "nrBootstrap|n", &nrBootStrapSamples, "inputData|i", &inputData);
    if(inputData.length > 0)
      spectrum = inputData.split(",").map!"a.to!ulong()"().array();
    else {
      enforce(args.length == 2, "need exactly one binning file");
      auto inputFilename = args[1];
      auto inputFile = inputFilename == "-" ? stdin : File(inputFilename, "r");
      spectrum = inputFile.byLine.map!"a.to!ulong"().array;
    }
    stderr.writeln("read spectrum: ", spectrum);
}

void displayHelpMessage() {
    stderr.writeln("singleSpectrumFit [options] <inputFile>
        Choose \"-\" to read input from stdin
        -m, --model <model>                 fit model for synonymous sites. Must be one of:
                                            (driverfieldNeutral, driverfieldNeutralSel, unlinkedNeutral, 
                                             unlinkedNeutralSel, constrainedNonNeutral, constrainedSel) 
                                            [default: driverfieldNeutral]
        --maxSteps=<int>                    maximum number of iterations of Powell's method for minimization
                                            [default: 200]
        --noBootstrap
        -i, --inputData                     input data given as comma-separated list of integers
        --spectrumFile                      file to write the spectrum to
        -n, --nrBootstrap                   nr of bootstrap samples
        --mu
        --V
        --t");
}

void run() {
    auto scoreFunc = getScoreFunc(spectrum, model);
    
    auto xMin = getMinimum(scoreFunc);
    writeSpectrum(spectrumFilename, scoreFunc, xMin);
    stderr.writeln("params: ", xMin);
    auto substLoadStats = getSubstAndLoadStats(scoreFunc, xMin);
    
    if(noBootstrap) {
        writeParams(scoreFunc, xMin, substLoadStats);
    }
    else {
        auto accumulator = new Accumulator(scoreFunc.nrParams);
        auto substLoadAccumulator = new Accumulator(7);
        foreach(i; 0 .. nrBootStrapSamples) {
            stderr.writefln("Bootstrap sample %s", i);
            auto bootstrapSpectrum = getBootstrap(spectrum);
            stderr.writeln("bootstrap spectrum: ", bootstrapSpectrum);
            auto bootstrapScoreFunc = getScoreFunc(bootstrapSpectrum, model);
            auto xMinBootstrap = getMinimum(bootstrapScoreFunc);
            stderr.writeln("params: ", xMinBootstrap);
            auto substLoadStatsBS = getSubstAndLoadStats(bootstrapScoreFunc, xMinBootstrap);
            accumulator.add(xMinBootstrap);
            substLoadAccumulator.add(substLoadStatsBS);
        }
        writeParamsWithBootstrap(scoreFunc, xMin, substLoadStats, accumulator, substLoadAccumulator);
    }
}  

SingleSpectrumScore getScoreFunc(ulong[] spectrum, model_t model) {
    SingleSpectrumScore scoreFunc;
    final switch(model) {
        case model_t.driverfieldNeutral:
        scoreFunc = new DriverfieldNeutralScore(spectrum);
        break;
        case model_t.driverfieldNeutralSel:
        scoreFunc = new DriverfieldNeutralScoreSel(spectrum);
        break;
        case model_t.unlinkedNeutral:
        scoreFunc = new UnlinkedNeutralScore(spectrum);
        break;
        case model_t.unlinkedNeutralSel:
        scoreFunc = new UnlinkedNeutralScoreSel(spectrum);
        break;
        case model_t.constrainedNonNeutral:
        scoreFunc = new ConstrainedNonNeutralScore(spectrum, mu, V, t);
        break;
        case model_t.constrainedSel:
        scoreFunc = new ConstrainedSelScore(spectrum, mu, V, t);
    }
    return scoreFunc;
}

double[] getMinimum(SingleSpectrumScore scoreFunc) {
    auto powell = new Powell!SingleSpectrumScore(scoreFunc);
    auto xInitial = scoreFunc.initialParams();
    powell.init(xInitial);
    
    double[] x;
    while(powell.iter < maxSteps && !powell.finished()) {
        stderr.writef("\rScore Minimization [%s/%s (max)]", powell.iter, maxSteps);
        x = powell.step();
    }
    stderr.writeln("\nScore:", scoreFunc(x));
    return x;
}

void writeParams(SingleSpectrumScore scoreFunc, in double[] x, in double[] substLoadStats) {
    auto p = scoreFunc.makeSingleSpectrumParams(x);
    auto score = scoreFunc.getScore(p);
    foreach(key; ["mu", "V", "t", "s", "gamma", "c", "cw"])
        writefln("%s\t%s", key, p[key]);
    foreach(i, key; ["neutralDrift", "weakSelDrift", "neutralHH", "weakSelHH", "adaptive", "driftLoad", "hhLoad"])
        writefln("%s\t%s", key, substLoadStats[i]);
    writefln("score\t%s", score);
}

void writeParamsWithBootstrap(SingleSpectrumScore scoreFunc, in double[] xMin, in double[] substLoadStats,
                              Accumulator accumulator, Accumulator substLoadAccumulator) {
    auto p = scoreFunc.makeSingleSpectrumParams(xMin);
    auto accMean = accumulator.mean();
    auto accStddev = accumulator.stddev();
    auto paramNames = scoreFunc.paramNames();
    auto pBSmean = scoreFunc.makeSingleSpectrumParams(accMean);
    auto pBSstddev = pBSmean.dup;
    foreach(key, ref val; pBSstddev)
        val = 0.0;
    foreach(i, param; paramNames) {
        pBSstddev[param] = accStddev[i];
    }
    auto substLoadMean = substLoadAccumulator.mean();
    auto substLoadStddev = substLoadAccumulator.stddev();
    
    auto score = scoreFunc.getScore(p);
    foreach(key; ["mu", "V", "t", "s", "gamma", "c", "cw"])
        writefln("%s\t%s\t%s\t%s", key, p[key], pBSmean[key], pBSstddev[key]);
    foreach(i, key; ["neutralDrift", "weakSelDrift", "neutralHH", "weakSelHH", "adaptive", "driftLoad", "hhLoad"])
        writefln("%s\t%s\t%s\t%s", key, substLoadStats[i], substLoadMean[i], substLoadStddev[i]);
    writefln("score\t%s", score);
}

ulong[] getBootstrap(in ulong[] spectrum) {
    auto norm = spectrum.reduce!"a+b"();
    stderr.writef("bootstrapping from %s data points...", norm);
    auto ret = new ulong[spectrum.length];
    foreach(i; 0 .. norm) {
        auto k = dice(spectrum);
        ret[k] += 1;
    }
    stderr.writeln("done");
    return ret;
}

double[7] getSubstAndLoadStats(SingleSpectrumScore scoreFunc, in double[] xMin) {
    auto p = scoreFunc.makeSingleSpectrumParams(xMin);
    auto m = cast(int)scoreFunc.m;
    auto mu = p["mu"];
    auto t = p["t"];
    auto V = p["V"];
    auto s = p["s"];
    auto gamma = p["gamma"];
    auto c = p["c"];
    auto cw = p["cw"];
    
    // definition: double Q1vsM(int k, int m, double t, double mu, double s, double gamma, double V)
    auto neutral = c * Q1vsM(m, m, t, mu, 0, 0, 0);
    auto neutralAndHH = c * Q1vsM(m, m, t, mu, 0, 0, V);
    auto driftOnly = neutral * exp(-V);
    auto hh = neutralAndHH - driftOnly;

    auto selDrift = (1.0 - c) * cw * Q1vsM(m, m, t, mu, s, 0, 0);
    auto selDriftAndHH = (1.0 - c) * cw * Q1vsM(m, m, t, mu, s, 0, V);
    auto selHH = selDriftAndHH - selDrift;
  
    auto adaptive = (1.0 - c) * (1.0 - cw) * 2.0 * gamma * t;
  
    auto totalSubst = neutralAndHH + selDriftAndHH + adaptive;
  
    auto driftLoad = selDrift * s;
    auto hhLoad = selHH * s;
    
    return [driftOnly, selDrift, hh, selHH, adaptive, driftLoad, hhLoad];
}

void writeSpectrum(string spectrumFile, SingleSpectrumScore scoreFunc, double[] xmin) {
  auto p = scoreFunc.makeSingleSpectrumParams(xmin);
  auto spectrum = scoreFunc.getSingleSpectrumProbs(p);
  auto f = File(spectrumFile, "w");
  foreach(x; spectrum)
    f.writeln(x);
}


class Accumulator {
    ubyte nrParams;
    double[] sums;
    double[] sums_sq;
    double norm;
    double[][] inputs;
    
    this(ubyte nrParams) {
        this.nrParams = nrParams;
        sums = new double[nrParams];
        sums[] = 0.0;
        sums_sq = new double[nrParams];
        sums_sq[] = 0.0;
        norm = 0.0;
    }
    
    void add(in double[] x)
    in {
        assert(x.length == nrParams);
    }
    body {
        inputs ~= x.dup;
        foreach(i, xx; x) {
            sums[i] += xx;
            sums_sq[i] += xx * xx;
        }
        norm += 1.0;
    }
  
    double[] mean() const {
        return iota(nrParams).map!(i => sums[i] / norm)().array.dup;
    }
    
    double[] stddev() const {
        return iota(nrParams).map!(i => sqrt(max(sums_sq[i] / norm - (sums[i] / norm)^^2, 0.0)))().array.dup;
    }
}