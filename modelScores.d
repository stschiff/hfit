import std.math;
import std.algorithm;
import std.stdio;
import std.json;
import std.string;
import std.exception;
import popGenFunc;

immutable static double penalty = 1.0e20;

class SingleSpectrumScore {
  SingleSpectrumProb neutralSpectrumProbFunc;
  SingleSpectrumProb selectedSpectrumProbFunc;
  const(int[]) spectrumData;
  int m;
  int norm;
  int nrParams;
  
  this(in int[] spectrumData, int nrParams) {
    enforce(spectrumData.length > 0, "need non-empty data");
    m = cast(int)spectrumData.length - 1;
    neutralSpectrumProbFunc = new SingleSpectrumProb(m);
    selectedSpectrumProbFunc = new SingleSpectrumProb(m);
    norm = reduce!"a+b"(spectrumData.dup);
    this.spectrumData = spectrumData;
    this.nrParams = nrParams;
  }
  
  abstract double[string] makeSingleSpectrumParams(double[] x);
  abstract double[] makeVecFromParams(double[string]);
  abstract double[] initialParams();
  abstract string[] paramNames();
  
  void checkParams(double[string] p) {
    if(p["c"] < 0.0 || p["c"] > 1.0)
      throw new IllegalParametersException("c out of range");
    if(p["cw"] < 0.0 || p["cw"] > 1.0)
      throw new IllegalParametersException("cw out of range");
    if(p["gamma"] < 0.0 || 2.0 * p["gamma"] * p["t"] > 1.0)
      throw new IllegalParametersException("gamma out of range");
    if(p["s"] < 0.0 || p["s"] > 150.0)
      throw new IllegalParametersException("selection out of range");
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
    foreach(k; 0 .. m + 1) {
      score += -log(probs[k]) * spectrumData[k];
    }
    return score;
  }
  
  double opCall(double[string] p) {
    double score;
    try{
      score = getScore(p);
    }
    catch(IllegalParametersException e) {
      debug stderr.writeln(e.msg, ", returning penalty");
      return penalty;
    }
    catch(GSLNumericsException e) {
      debug stderr.writeln(e.msg, ", returning penalty");
      return penalty;
    }
    if(isnan(score)) {
      return penalty;
    }
    return score;
  }
  
  double opCall(double[] x) {
    assert(x.length == nrParams);
    auto p = makeSingleSpectrumParams(x);
    return opCall(p);
  }
}

class FullScore : SingleSpectrumScore {
  
  this(in int[] spectrum) {
    super(spectrum, 7);
  }
  
  override double[string] makeSingleSpectrumParams(double[] x) {
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


class DriverfieldNeutralScoreSel : SingleSpectrumScore {

  this(in int[] spectrum) {
    super(spectrum, 4);
  }
  
  override double[string] makeSingleSpectrumParams(double[] x) {
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

class DriverfieldNeutralScore : SingleSpectrumScore {

  this(in int[] spectrum) {
    super(spectrum, 3);
  }
  
  override double[string] makeSingleSpectrumParams(double[] x) {
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


class UnlinkedNeutralScoreSel : SingleSpectrumScore {
  
  this(in int[] spectrum) {
    super(spectrum, 3);
  }
  
  override double[string] makeSingleSpectrumParams(double[] x) {
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

class UnlinkedNeutralScore : SingleSpectrumScore {
  
  this(in int[] spectrum) {
    super(spectrum, 2);
  }
  
  override double[string] makeSingleSpectrumParams(double[] x) {
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

class ConstrainedNonNeutralScore : SingleSpectrumScore {
  
  double[string] params;
  
  this(in int[] spectrum, double mu, double V, double t) {
    super(spectrum, 4);
    params["mu"] = mu;
    params["V"] = V;
    params["t"] = t;
  }
  
  override void checkParams(double[string] p) {
    super.checkParams(p);
    if(p["s"] < 1.0)
      throw new IllegalParametersException("selection out of range");
  }
  
  override string[] paramNames() {
    return ["gamma", "c", "s", "cw"];
  }
  
  override double[string] makeSingleSpectrumParams(double[] x) {
    auto gamma = x[0];
    auto c = x[1];
    auto s = x[2];
    auto cw = x[3];
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

