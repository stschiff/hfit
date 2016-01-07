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
import mcmc;

enum model_t {LINKED_SELECTION, LINKED_SELECTION_WITHS, REDUCED_NE, REDUCED_NE_WITHS, MIXED, MIXED_SIMPLE}

string neutralFitFileName = "";
model_t model = model_t.LINKED_SELECTION;
auto maxSteps = 500U;
double theta, nu, tau;
string spectrumFilename = "/dev/null";
string inputData;
ulong[] spectrum;
bool noMCMC;
auto nrMCMCcycles = 1000;
string mcmcTraceFile = "/dev/null";

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
    getopt(args, std.getopt.config.passThrough,
           "model|m", &model,
           "maxSteps", &maxSteps,
           "theta", &theta,
           "nu", &nu,
           "tau", &tau,
           "spectrumFile", &spectrumFilename,
           "noMCMC", &noMCMC,
           "nrMCMCcycles|n", &nrMCMCcycles,
           "mcmcTraceFile", &mcmcTraceFile,
           "inputData|i", &inputData);
   
    if(inputData.length > 0)
        spectrum = inputData.split(",").map!"a.to!ulong()"().array();
    else {
        enforce(args.length == 2, "need exactly one binning file");
        auto inputFilename = args[1];
        auto inputFile = inputFilename == "-" ? stdin : File(inputFilename, "r");
        spectrum = inputFile.byLine.map!"a.to!ulong"().array;
    }
    
    if(model == model_t.MIXED || model == model_t.MIXED_SIMPLE) {
        enforce(tau > 0.0, "need --tau");
        enforce(theta > 0.0, "need --theta");
    }
    
    stderr.writeln("read spectrum: ", spectrum);
}

void displayHelpMessage() {
    stderr.writeln("singleSpectrumFit [options] <inputFile>
        Choose \"-\" to read input from stdin
        -m, --model <model>                 fit model. Must be one of the following:
                                            * LINKED_SELECTION:
                                                Free parameters: theta, nu, tau
                                            * LINKED_SELECTION_WITHS:
                                                Free parameters: theta, nu, tau, sigma
                                            * REDUCED_NE:
                                                Free parameters: theta, tau
                                            * REDUCED_NE_WITHS:
                                                Free parameters: theta, tau, sigma
                                            * MIXED:
                                                Free parameters: cn, cw, ca, sigma
                                                Additional Input: theta, nu, tau
                                            * MIXED_SIMPLE:
                                                Free parameters: cn, ca
                                                AdditionalInput: theta, nu, tau
                                               
                                            [default: LINKED_SELECTION]
        --maxSteps=<int>                    maximum number of iterations of Powell's method for minimization
                                            [default: 200]
        --noMCMC                            no MCMC
        -i, --inputData                     input data given as comma-separated list of integers
        --spectrumFile                      file to write the spectrum to
        -n, --nrMCMCcycles                  nr of MCMC main cycles [default: 1000]
        --tau                               input for fixed model parameters
        --theta                             input for fixed model parameters
        --nu                                input for fixed model parameters
        --mcmcTraceFile                     file to write MCMC trace to

        The output of the program contains columns with model estimates. The first column is the name of the parameter. The second column is the maximum likelihood estimate. The third column is the standard deviation for the parameter, obtained from boostrapping (omitted if --noMCMC is given)\n");
}

void run() {
    auto scoreFunc = getScoreFunc(spectrum, model);
    
    auto xMin = findMinimum(scoreFunc);
    stderr.writeln("found minimum: ", xMin);
    writeSpectrum(spectrumFilename, scoreFunc, xMin);
    auto pMin = getNamedParams(scoreFunc, xMin);

    if(noMCMC) {
        writeParams(pMin);
    }
    else {
        stderr.writeln("starting MCMC");
        auto mcmc = new MCMC!SingleSpectrumScore(scoreFunc, xMin);
        auto accumulator = new Accumulator(pMin.keys());
        mcmc.run(nrMCMCcycles);
        auto nrBurninSteps = mcmc.getNrBurninCycles();
        stderr.writefln("done MCMC, with %s burnin cycles and %s main cycles", nrBurninSteps, mcmc.trace.length - nrBurninSteps);
        foreach(i; nrBurninSteps .. mcmc.trace.length) {
            auto x = mcmc.trace[i];
            auto pMinCycle = getNamedParams(scoreFunc, x);
            accumulator.add(pMinCycle);
        }
        writeParamsWithBootstrap(pMin, accumulator);
        mcmc.reportTraces(mcmcTraceFile);
    }
}  

SingleSpectrumScore getScoreFunc(ulong[] spectrum, model_t model) {
    SingleSpectrumScore scoreFunc;
    final switch(model) {
        case model_t.LINKED_SELECTION:
        scoreFunc = new LinkedSelection_Score(spectrum);
        break;
        case model_t.LINKED_SELECTION_WITHS:
        scoreFunc = new LinkedSelection_withS_Score(spectrum);
        break;
        case model_t.REDUCED_NE:
        scoreFunc = new ReducedNe_Score(spectrum);
        break;
        case model_t.REDUCED_NE_WITHS:
        scoreFunc = new ReducedNe_withS_Score(spectrum);
        break;
        case model_t.MIXED:
        scoreFunc = new MixedScore(spectrum, theta, nu, tau);
        break;
        case model_t.MIXED_SIMPLE:
        scoreFunc = new MixedSimpleScore(spectrum, theta, nu, tau);
        break;
    }
    return scoreFunc;
}

double[] findMinimum(SingleSpectrumScore scoreFunc) {
    
    auto powell = new Powell!SingleSpectrumScore(scoreFunc);
    auto xInitial = scoreFunc.initialParams();
    powell.init(xInitial);
    
    double[] x;
    while(powell.iter < maxSteps && !powell.finished()) {
        stderr.writef("\rScore Minimization [%s/%s (max)]", powell.iter, maxSteps);
        x = powell.step();
    }
    return x;
}

double[string] getNamedParams(SingleSpectrumScore scoreFunc, double[] x) {
    auto p = scoreFunc.makeSingleSpectrumParams(x);
    auto score = scoreFunc(x);
    
    auto mu = p["mu"];
    auto t = p["t"];
    auto V = p["V"];
    auto s = p["s"];
    auto gamma = p["gamma"];
    auto c = p["c"];
    auto cw = p["cw"];
    
    // the internals of this software use a different parameterization, using 2N, not 2N_0 as the scaling factor. Also there is a bit different usage of cn, cw and ca (see modelScores.d). That's why we scale everything to match the notation in the preprint (Schiffels, Mustonen, Laessig, 2014)
    double[string] ret;
    ret["tau"] = tau;
    ret["score"] = score;
    ret["theta"] = mu;
    ret["nu"] = V;
    ret["sigma"] = s;
    ret["cn"] = c;
    ret["cw"] = (1.0 - c) * cw;
    ret["ca"] = (1.0 - c) * (1.0 - cw) * 2.0 * gamma * t;
    
    // definition: double Q1vsM(int k, int m, double t, double mu, double s, double gamma, double V)
    auto m = cast(int)scoreFunc.m;
    auto neutral = c * Q1vsM(m, m, t, mu, 0, 0, 0);
    auto neutralAndHH = c * Q1vsM(m, m, t, mu, 0, 0, V);
    auto driftOnly = neutral * exp(-V);
    auto hh = neutralAndHH - driftOnly;

    auto selDrift = (1.0 - c) * cw * Q1vsM(m, m, t, mu, s, 0, 0);
    auto selDriftAndHH = (1.0 - c) * cw * Q1vsM(m, m, t, mu, s, 0, V);
    auto selDriftOnly = selDrift * exp(-V);
    auto selHH = selDriftAndHH - selDriftOnly;
  
    auto adaptive = (1.0 - c) * (1.0 - cw) * 2.0 * gamma * t;

    auto driftLoad = (1.0 - c) * cw * 2.0 * u0(-s, mu) / (u0(s, mu) + u0(-s, mu));
    auto hhAndDriftLoad = (1.0 - c) * cw * 2.0 * uLinkage(-s, mu, V) / (uLinkage(s, mu, V) + uLinkage(-s, mu, V));
    auto hhLoad = hhAndDriftLoad - driftLoad;
    ret["neutralDrift"] = driftOnly;
    ret["weakSelDrift"] = selDrift;
    ret["neutralHH"] = hh;
    ret["weakSelHH"] = selHH;
    ret["adaptive"] = adaptive;
    ret["driftLoad"] = driftLoad;
    ret["hhLoad"] = hhLoad;
    
    return ret;
}


void writeParams(in double[string] p) {
    writefln("Parameter\tMaximumLL");
    foreach(k, v; p)
        writefln("%s\t%s", k, v);
}

void writeParamsWithBootstrap(in double[string] p, in Accumulator accumulator) {
    writefln("Parameter\tMaximumLL\tMCMCmedian\tMCMClowerCI\tMCMCupperCI");
    foreach(k, v; p) {
        auto orderStats = accumulator.orderStats(k);
        writefln("%s\t%s\t%s\t%s\t%s", k, v, orderStats[1], orderStats[0], orderStats[2]);
    }
}

void writeSpectrum(string spectrumFile, SingleSpectrumScore scoreFunc, double[] xMin) {
  
    auto p = scoreFunc.makeSingleSpectrumParams(xMin);

    auto spectrum = scoreFunc.getSingleSpectrumProbs(p);
    auto f = File(spectrumFile, "w");
    foreach(x; spectrum)
        f.writeln(x);
}

class Accumulator {
    double[][string] values;

    this(string[] keys) {
        foreach(k; keys)
            values[k] = [];
    }

    void add(in double[string] p)
    {
        foreach(k, v; p)
            values[k] ~= v;
    }
    
    double[3] orderStats(string key) const {
        auto sorted = values[key].dup;
        sorted.sort();
        auto midIndex = to!int(sorted.length / 2);
        auto lowerCIindex = to!int(sorted.length * 0.025);
        auto upperCIindex = to!int(sorted.length * 0.975);
        return [sorted[lowerCIindex], sorted[midIndex], sorted[upperCIindex]];
    }
}