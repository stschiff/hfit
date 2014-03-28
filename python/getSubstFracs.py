import std.getopt;
import std.stdio;
import std.math;
import popGenFunc;

void main(string[] args) {
  double mu, t, V, gamma, c, cw, s;
  int m;
  bool normalized = false;

  getopt(args, "mu|u", &mu, "t|t", &t, "gamma|g", &gamma, "c|c", &c, "s|s", &s, 
         "V|V", &V, "m|m", &m, "cw|w", &cw, "normalized", &normalized);
  
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
  
  // definition: double uLinkage(double sigma, double mu, double V)
  auto uP = uLinkage(s, mu, V);
  auto uM = uLinkage(-s, mu, V);
  
  writeln("neutralDrift\twSelDrift\tneutralHH\twSelHH\tadaptive\tdriftLoad\thhLoad");
  if(normalized)
    writefln("%s\t%s\t%s\t%s\t%s\t%s\t%s", driftOnly / totalSubst, selDrift / totalSubst, hh / totalSubst, selHH / totalSubst, adaptive / totalSubst, driftLoad, hhLoad);
  else
    writefln("%s\t%s\t%s\t%s\t%s\t%s\t%s", driftOnly, selDrift, hh, selHH, adaptive, driftLoad, hhLoad);
}

