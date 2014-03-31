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

enum model_t {driverfieldNeutral, driverfieldNeutralSel, unlinkedNeutral, unlinkedNeutralSel, constrainedNonNeutral}

string neutralFitFileName = "";
model_t model = model_t.driverfieldNeutral;
auto maxSteps = 500U;
double mu, V, t; // for constrained model
ulong[] spectrum;
string spectrumFilename;
bool noBootstrap;
auto nrBootStrapSamples = 4U;

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
    getopt(args, std.getopt.config.passThrough, "model|m", &model, "maxSteps", &maxSteps, "mu", &mu, "V", &V, "t", &t, "spectrumFile", &spectrumFilename, "noBootstrap", &noBootstrap);
    enforce(args.length == 2, "need exactly one binning file");
    auto inputFilename = args[1];
    auto inputFile = inputFilename == "-" ? stdin : File(inputFilename, "r");
    spectrum = inputFile.byLine.map!"a.to!ulong"().array;
    stderr.writeln("read spectrum: ", spectrum);
    enforce(inputFilename == "-" || isFile(inputFilename), "illegal input file");
}

void displayHelpMessage() {
    stderr.writeln("singleSpectrumFit [options] <inputFile>
        Choose \"-\" to read input from stdin
        -m, --model <model> : fit model for synonymous sites. Must be one of:
              (driverfieldNeutral, driverfieldNeutralSel, unlinkedNeutral, unlinkedNeutralSel, constrainedNonNeutral) 
              [default: driverfieldNeutral]
        --maxSteps=<int> maximum number of iterations of Powell's method for minimization [default: 200]
        --noBootstrap
        --mu
        --V
        --t");
}

void run() {
    auto scoreFunc = getScoreFunc(spectrum, model);
    
    auto xMin = getMinimum(scoreFunc);
    
    if(noBootstrap) {
        writeParams(scoreFunc, xMin);
    }
    else {
        auto accumulator = new Accumulator(scoreFunc.nrParams);
        foreach(i; 0 .. nrBootStrapSamples) {
            stderr.writefln("Bootstrap sample %s", i);
            auto bootstrapSpectrum = getBootstrap(spectrum);
            writeln("bootstrap spectrum: ", bootstrapSpectrum);
            auto bootstrapScoreFunc = getScoreFunc(bootstrapSpectrum, model);
            auto xMinBootstrap = getMinimum(bootstrapScoreFunc);
            accumulator.add(xMinBootstrap);
        }
        writeParamsWithBootstrap(scoreFunc, xMin, accumulator);
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
    stderr.write("\n");
    return x;
}

void writeParams(SingleSpectrumScore scoreFunc, double[] x) {
    auto p = scoreFunc.makeSingleSpectrumParams(x);
    auto score = scoreFunc.getScore(p);
    foreach(key; ["mu", "V", "t", "s", "gamma", "c", "cw"])
        writefln("%s\t%s", key, p[key]);
    writefln("score\t%s", score);
}

void writeParamsWithBootstrap(SingleSpectrumScore scoreFunc, double[] xMin, Accumulator accumulator) {
    auto p = scoreFunc.makeSingleSpectrumParams(xMin);
    auto pBSmean = scoreFunc.makeSingleSpectrumParams(accumulator.mean());
    auto pBSstddev = scoreFunc.makeSingleSpectrumParams(accumulator.stddev());
    
    auto score = scoreFunc.getScore(p);
    foreach(key; ["mu", "V", "t", "s", "gamma", "c", "cw"])
        writefln("%s\t%s\t%s\t%s", key, p[key], pBSmean[key], pBSstddev[key]);
    writefln("score\t%s", score);
}

ulong[] getBootstrap(ulong[] spectrum) {
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

class Accumulator {
    ubyte nrParams;
    double[] sums;
    double[] sums_sq;
    double norm;
    
    this(ubyte nrParams) {
        this.nrParams = nrParams;
        sums = new double[nrParams];
        sums[] = 0.0;
        sums_sq = new double[nrParams];
        sums_sq[] = 0.0;
        norm = 0.0;
    }
    
    void add(double[] x)
    in {
        assert(x.length == nrParams);
    }
    body {
        foreach(i, xx; x) {
            sums[i] += xx;
            sums_sq[i] += xx * xx;
        }
        norm += 1.0;
    }
  
    double[] mean() {
        return iota(nrParams).map!(i => sums[i] / norm)().array;
    }
    
    double[] stddev() {
        return iota(nrParams).map!(i => sqrt(sums_sq[i] / norm - (sums[i] / norm)^^2))().array;
    }
}