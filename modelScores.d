import std.math;
import std.algorithm;
import std.stdio;
import std.json;
import std.string;
import std.exception;
import popGenFunc;

immutable static double penalty = 1.0e20;

// c is the fraction of neutral sites
// (1-c)*cw is the fraction of weakly selected sites
// (1-c)(1-cw)*2*gamma*t is the fraction of adaptive sites.
// these definitions are different than in the paper, due to easier contraints in the optimization: c, cw and 2*gamma*t are all free between 0 and 1. The usage of gamma is historical, see Mustonen and Laessig, 2007.

class SingleSpectrumScore {
  SingleSpectrumProb neutralSpectrumProbFunc;
  SingleSpectrumProb selectedSpectrumProbFunc;
  const(ulong[]) spectrumData;
  ulong m;
  ulong norm;
  ubyte nrParams;
  int minAc;
  
  this(in ulong[] spectrumData, ubyte nrParams, int minAc=1) {
    enforce(spectrumData.length > 0, "need non-empty data");
    m = spectrumData.length - 1;
    neutralSpectrumProbFunc = new SingleSpectrumProb(cast(int)m);
    selectedSpectrumProbFunc = new SingleSpectrumProb(cast(int)m);
    norm = reduce!"a+b"(spectrumData.dup);
    this.spectrumData = spectrumData;
    this.nrParams = nrParams;
    this.minAc = minAc;
  }
  
  abstract double[string] makeSingleSpectrumParams(in double[] x);
  abstract double[] makeVecFromParams(double[string]);
  abstract double[] initialParams();
  abstract string[] paramNames();
  
  void checkParams(double[string] p) {
    if(p["c"] < 0.0 || p["c"] > 1.0)
      throw new IllegalParametersException("c out of range");
    if(p["cw"] < 0.0 || p["cw"] > 1.0)
      throw new IllegalParametersException("cw out of range");
    if(p["gamma"] < 0.0 || 2.0 * p["gamma"] * p["t"] > 1.0 + 1.0e-8)
      throw new IllegalParametersException(format("2*gamma*t=%s out of range", 2.0 * p["gamma"] * p["t"]));
    if(p["s"] < 0.0 || p["s"] > 150.0)
      throw new IllegalParametersException(format("selection out of range:%s", p["s"]));
  }
  
  immutable(double)[] getSingleSpectrumProbs(double[string] p) {
    immutable(double)[] ret;
    if(norm > 0) {
      checkParams(p);
      auto Q0 = neutralSpectrumProbFunc(p["t"], p["mu"], 0.0, 0.0, p["V"]);
      auto Qs = selectedSpectrumProbFunc(p["t"], p["mu"], p["s"], 0.0, p["V"]);
      auto neutral_part = Q0.dup;
      neutral_part[] *= p["c"];
      auto weak_sel_part = Qs.dup;
      weak_sel_part[] *= (1.0 - p["c"]) * p["cw"];
      auto strong_sel_part = new double[m + 1];
      strong_sel_part[] = 0.0;
      strong_sel_part[0] = 1.0 - 2.0 * p["gamma"] * p["t"];
      strong_sel_part[$ - 1] = 2.0 * p["gamma"] * p["t"];
      strong_sel_part[] *= (1.0 - p["c"]) * (1.0 - p["cw"]);
      auto full = neutral_part.dup;
      full[] += weak_sel_part[];
      full[] += strong_sel_part[];
      ret = full.idup;
    }
    else {
      auto dummy = new double[m + 1];
      dummy[] = 1.0; // dummy value
      ret = dummy.idup;
    }
    return ret;
  }
  
  double getScore(double[string] p) {
    auto probs = getSingleSpectrumProbs(p);
    auto score = 0.0;
    auto norm = spectrumData.reduce!"a+b"();
    auto probTotal = 0.0;
    auto countTotal = 0;
    foreach(k; minAc .. m + 1) {
      score += -log(probs[k]) * spectrumData[k];
      probTotal += probs[k];
      countTotal += spectrumData[k];
    }
    score += -log(1.0 - probTotal) * (norm - countTotal);
    return score;
  }
  
  double opCall(double[string] p) {
    double score;
    try{
      score = getScore(p);
    }
    catch(IllegalParametersException e) {
      return penalty;
    }
    catch(GSLNumericsException e) {
      return penalty;
    }
    if(isNaN(score)) {
      return penalty;
    }
    return score;
  }
  
  double opCall(in double[] x) {
    assert(x.length == nrParams);
    auto p = makeSingleSpectrumParams(x);
    return opCall(p);
  }
}

class FullScore : SingleSpectrumScore {
  
  this(in ulong[] spectrum, int minAc=1) {
    super(spectrum, 7, minAc);
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto mu = x[0];
    auto V = x[1];
    auto t = x[2];
    auto gamma = x[3];
    auto c = x[4];
    auto s = x[5];
    auto cw = x[6];
    double[string] params;
    params["mu"] = mu;
    params["t"] = t;
    params["V"] = V;
    params["gamma"] = gamma;
    params["c"] = c;
    params["s"] = s;
    params["cw"] = cw;
    return params;
  }
  
  override string[] paramNames() {
    return ["mu", "t", "V", "gamma", "c", "s", "cw"];
  }
  
  override double[] initialParams() {
    auto ret = new double[7];
    ret[0] = 0.01;
    ret[1] = 0.1;
    ret[2] = 4;
    ret[3] = 0.0;
    ret[4] = 0.5;
    ret[5] = 1.0;
    ret[6] = 0.5;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[7];
    ret[0] = p["mu"];
    ret[1] = p["V"];
    ret[2] = p["t"];
    ret[3] = p["gamma"];
    ret[4] = p["c"];
    ret[5] = p["s"];
    ret[6] = p["cw"];
    return ret;
  } 
}


class LinkedSelection_withS_Score : SingleSpectrumScore {

  this(in ulong[] spectrum, int minAc=1) {
    super(spectrum, 4, minAc);
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto mu = x[0];
    auto V = x[1];
    auto t = x[2];
    auto s = x[3];
    double[string] params;
    params = ["mu":mu, "t":t, "V":V, "s":s, "gamma":0, "c":0.0, "cw":1.0];
    return params;
  }
  
  override string[] paramNames() {
    return ["mu", "t", "V", "s"];
  }
  
  override double[] initialParams() {
    auto ret = new double[4];
    ret[0] = 0.01;
    ret[1] = 0.1;
    ret[2] = 4.0;
    ret[3] = 0.0;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[4];
    ret[0] = p["mu"];
    ret[1] = p["V"];
    ret[2] = p["t"];
    ret[3] = p["s"];
    return ret;
  }
  
}

class LinkedSelection_Score : SingleSpectrumScore {

  this(in ulong[] spectrum, int minAc=1) {
    super(spectrum, 3, minAc);
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto mu = x[0];
    auto V = x[1];
    auto t = x[2];
    double[string] params;
    params = ["mu":mu, "t":t, "V":V, "s":0.0, "gamma":0.0, "c":1.0, "cw":0.0];
    return params;
  }
  
  override string[] paramNames() {
    return ["mu", "t", "V"];
  }
  
  override double[] initialParams() {
    auto ret = new double[3];
    ret[0] = 0.01;
    ret[1] = 0.1;
    ret[2] = 4;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[3];
    ret[0] = p["mu"];
    ret[1] = p["V"];
    ret[2] = p["t"];
    return ret;
  }
}

class ReducedNe_withS_Score : SingleSpectrumScore {
  
  this(in ulong[] spectrum, int minAc=1) {
    super(spectrum, 3, minAc);
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto mu = x[0];
    auto t = x[1];
    auto s = x[2];
    double[string] params;
    params = ["mu":mu, "t":t, "s":s, "V":0.0, "gamma":0.0, "c":0.0, "cw":1.0];
    return params;
  }
  
  override string[] paramNames() {
    return ["mu", "t", "s"];
  }
  
  override double[] initialParams() {
    auto ret = new double[3];
    ret[0] = 0.01;
    ret[1] = 4.0;
    ret[2] = 0.0;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[3];
    ret[0] = p["mu"];
    ret[1] = p["t"];
    ret[2] = p["s"];
    return ret;
  }
  
}

class ReducedNe_Score : SingleSpectrumScore {
  
  this(in ulong[] spectrum, int minAc=1) {
    super(spectrum, 2, minAc);
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto mu = x[0];
    auto t = x[1];
    double[string] params;
    params = ["mu":mu, "t":t, "s":0.0, "V":0.0, "gamma":0.0, "c":1.0, "cw":0.0];
    return params;
  }
  
  override string[] paramNames() {
    return ["mu", "t"];
  }
  
  override double[] initialParams() {
    auto ret = new double[2];
    ret[0] = 0.01;
    ret[1] = 4;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[2];
    ret[0] = p["mu"];
    ret[1] = p["t"];
    return ret;
  }
  
}

class MixedScore : SingleSpectrumScore {
  
  double mu, V, t;
  
  this(in ulong[] spectrum, double mu, double V, double t, int minAc=1) {
    super(spectrum, 4, minAc);
    this.mu = mu;
    this.V = V;
    this.t = t;
  }
  
  override void checkParams(double[string] p) {
    super.checkParams(p);
    if(p["s"] < 1.0)
      throw new IllegalParametersException(format("selection out of range:%s", p["s"]));
  }
  
  override string[] paramNames() {
    return ["gamma", "c", "s", "cw"];
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto gamma = x[0];
    auto c = x[1];
    auto s = x[2];
    auto cw = x[3];
    double[string] params;
    params["mu"] = mu;
    params["V"] = V;
    params["t"] = t;
    params["gamma"] = gamma;
    params["c"] = c;
    params["s"] = s;
    params["cw"] = cw;
    return params;
  }
  
  override double[] initialParams() {
    auto ret = new double[4];
    ret[0] = 0.0;
    ret[1] = 0.5;
    ret[2] = 2.0;
    ret[3] = 0.5;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[4];
    ret[0] = p["gamma"];
    ret[1] = p["c"];
    ret[2] = p["s"];
    ret[3] = p["cw"];
    return ret;
  }
}

class MixedSimpleScore : SingleSpectrumScore {
  
  double mu, V, t;
  
  this(in ulong[] spectrum, double mu, double V, double t, int minAc=1) {
    super(spectrum, 2, minAc);
    this.mu = mu;
    this.V = V;
    this.t = t;
  }
  
  override string[] paramNames() {
    return ["gamma", "c"];
  }
  
  override double[string] makeSingleSpectrumParams(in double[] x) {
    auto gamma = x[0];
    auto c = x[1];
    double[string] params;
    params["mu"] = mu;
    params["V"] = V;
    params["t"] = t;
    params["gamma"] = gamma;
    params["c"] = c;
    params["s"] = 0.0;
    params["cw"] = 0.0;
    return params;
  }
  
  override double[] initialParams() {
    auto ret = new double[2];
    ret[0] = 0.0;
    ret[1] = 0.5;
    return ret;
  }
  
  override double[] makeVecFromParams(double[string] p) {
    auto ret = new double[2];
    ret[0] = p["gamma"];
    ret[1] = p["c"];
    return ret;
  }
}


