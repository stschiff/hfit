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
import modelScores;
import popGenFunc;
import powell;

enum model_t {driverfieldNeutral, driverfieldNeutralSel, unlinkedNeutral, unlinkedNeutralSel, constrainedNonNeutral}

string neutralFitFileName = "";
model_t model = model_t.driverfieldNeutral;
size_t maxSteps = 500;
double mu, V, t; // for constrained model
int[] spectrum;
string spectrumFilename;

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
  getopt(args, std.getopt.config.passThrough, "model|m", &model, "maxSteps", &maxSteps, "mu", &mu, "V", &V, "t", &t, "spectrumFile", &spectrumFilename);
  enforce(args.length == 2, "need exactly one binning file");
  auto inputFilename = args[1];
  auto inputFile = inputFilename == "-" ? stdin : File(inputFilename, "r");
  spectrum = inputFile.byLine.map!"a.to!int"().array;
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
      --mu
      --V
      --t");
}

void run() {
  auto scoreFunc = getScoreFunc(spectrum, model);
  
  auto xMin = getMinimum(scoreFunc);
  
  stdout.writeParamHeader();
  stdout.writeParams(scoreFunc, xMin);
  // if(!no_mcmc) {
  //   auto mcmcFile = File(mcmcStepFilename, "w");
  //   mcmcFile.writeParamHeader();
  //   auto sampler = new MCMCsampler(scoreFunc, xMin);
  //   auto accumulator = new Accumulator(scoreFunc.nrParams);
  //   foreach(loop; 0 .. mcmc_burnin + mcmc_steps) {
  //     auto x = sampler.next();
  //     mcmcFile.writeParams(scoreFunc, x);
  //     if(loop > mcmc_burnin)
  //       accumulator.add(x);
  //   }
  //   stdout.writeParams(scoreFunc, accumulator.mean());
  //   stdout.writeParamVar(scoreFunc, accumulator.stddev());
  // }
}  

SingleSpectrumScore getScoreFunc(int[] spectrum, model_t model) {
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

void writeParamHeader(File file) {
  file.writeln();
}

void writeParams(File file, SingleSpectrumScore scoreFunc, double[] x) {
  auto p = scoreFunc.makeSingleSpectrumParams(x);
  auto score = scoreFunc.getScore(p);
  file.writeln("Parameter\tMean(BS)\tStddev(BS)");
  foreach(key; ["mu", "V", "t", "s", "gamma", "c", "cw"])
    file.writefln("%s\t%s", key, p[key]);
  file.writefln("score\t%s", score);
}

void writeParamVar(File file, SingleSpectrumScore scoreFunc, double[] x) {
  auto p = ["mu":0.0, "V":0.0, "t":0.0, "s":0.0, "gamma":0.0, "c":0.0, "cw":0.0];
  auto names = scoreFunc.paramNames();
  foreach(i, n; names)
    p[n] = x[i];
  file.writefln("%s\t%s\t%s\t%s\t%s\t%s\t%s", p["mu"], p["V"], p["t"], p["s"], p["gamma"], p["c"], p["cw"]);
}

class Accumulator {
  int nrParams;
  double[] sums;
  double[] sums_sq;
  double norm;
  
  this(int nrParams) {
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