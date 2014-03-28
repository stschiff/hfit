import std.math;
import std.mathspecial;
import std.string;
import std.stdio;
import std.conv;

static this() {
  gsl_set_error_handler_off();
}

void assertApproxE(double val1, double val2, double tol1=0.0, double tol2=1.0e-12) {
  assertApproxE(val1, val2, tol1, tol2), text([val1, val2]);
}


class GSLNumericsException : Exception {
  this(string msg) {
    super(msg);
  }
}


extern (C) {
  struct gsl_sf_result {
    double val;
    double err;
  }
  
  int gsl_sf_hyperg_1F1_e(double a, double b, double x, gsl_sf_result* res);
  int gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, gsl_sf_result* res);

  void gsl_set_error_handler_off();
}

double hyperg_1F1(double a, double b, double x) {
  int err;
  gsl_sf_result res;
  err = gsl_sf_hyperg_1F1_e(a, b, x, &res);
  if(err) {
    throw new GSLNumericsException("hyperg_1F1 overflow");
  }
  return res.val;
}
double hyperg_2F1(double a, double b, double c, double x) {
  int err;
  gsl_sf_result res;
  err = gsl_sf_hyperg_2F1_e(a, b, c, x, &res);
  if(err) {
    throw new GSLNumericsException("hyperg_2F1 overflow");
  }
  return res.val;
}


// CODE FROM http://lists.gnu.org/archive/html/bug-gsl/2010-07/msg00000.html
// SEE Equation 32 in http://mathworld.wolfram.com/HypergeometricFunction.html
double hyperg_2F1_full(double a, double b, double c, double x) {
  double ret;
  double hyperg;
  if (-1.<=x && x<1.) {
    hyperg = hyperg_2F1(a, b, c, x);
    return hyperg;
  }
  else if (x<-1.) {
    // gsl_sf_hyperg_2F1 may have problems with negative arguments
    if (c-a<0) {
      hyperg = hyperg_2F1(a, c-b, c, x/(x-1.));
      hyperg *= pow(1.-x, -a);
      return hyperg;
    }
    else if (c-b<0) {
      hyperg = hyperg_2F1(c-a, b, c, x/(x-1.));
      hyperg *= pow(1.-x, -b);
      return hyperg;
    }
    // choose one of two equivalent formulas which is expected to be
    // more accurate
    else if (a*(c-b)<(c-a)*b) {
      hyperg = hyperg_2F1(a, c-b, c, x/(x-1.));
      hyperg *= pow(1.-x, -a);
      return hyperg;
    }
    else {
      hyperg = hyperg_2F1(c-a, b, c, x/(x-1.));
      hyperg *= pow(1.-x, -b);
      return hyperg;
    }
  }
  else {
    throw new GSLNumericsException("ERROR: x>=1 in hypgeo");
  }
  // return ret;
}

unittest {
  auto h = hyperg_1F1(1,2,3);
  assert(approxEqual(h, 6.361845641062556, 0.0, 1e-12));
  assert(approxEqual(hyperg_2F1_full(1,2,3,-4.), 0.2988202609457375, 0.0, 1e-12));
  assert(approxEqual(hyperg_2F1_full(1,2,3,-0.5), 0.7562791351346849, 0.0, 1e-12));
  assert(approxEqual(hyperg_2F1_full(1,2,3,0.7), 2.0570318543915747, 0.0, 1e-12));
  assert(approxEqual(hyperg_2F1_full(1, 0.3, 1.3, -3.0), 0.7113010112875268, 0.0, 1e-12));
}

immutable double eps = 1e-6;

double factln(int n) {
  static double a[100];
  if(n <= 1) return 0.0;
  if(n < 100) return !isNaN(a[n]) ? a[n] : (a[n] = logGamma(n + 1.0));
  return logGamma(n + 1.0);
}

double bico(int m, int k) {
  return floor(0.5 + exp(factln(m) - factln(k) - factln(m - k)));
}

double u0(double sigma, double mu) {
  double ret = mu;
  if(fabs(sigma) > eps) ret = mu * sigma / (1.0 - exp(-sigma));
  return ret;
}

unittest {
  assert(isNormal(u0(2, 0.025)));
  assert(isNormal(u0(0, 0.025)));
  assert(isNormal(u0(-2, 0.025)));
}

double uLinkage(double sigma, double mu, double V) {
  
  if(V < eps)
      return u0(sigma, mu);
  double ret = mu;
  if(sigma > eps) {
    double val = hyperg_2F1_full(1.0, V / sigma, 1.0 + V / sigma,
                                 1.0 - u0(sigma, mu) / mu);
    ret = u0(sigma, mu) * val;
  }
  else if(sigma < -eps)
      ret = mu / (fabs(sigma) + V) * (u0(sigma, mu) / mu * fabs(sigma) + V);
  return ret;
}

unittest {
  assert(isNormal(uLinkage(2, 0.025, 2)));
  assert(isNormal(uLinkage(0, 0.025, 2)));
  assert(isNormal(uLinkage(-2.0, 0.025, 2)));
}

// a is the current_allele {0,1}, e is the fitter allele {0,1}
double lambda(int a, int e, double uP, double uM, double gamma) {
  return (((a != e ? uM : uP) + gamma) / (uP + uM + 2.0 * gamma));
}

double zNeutral(double mu) {
  return exp(2 * logGamma(mu) - logGamma(2.0 * mu));
}

double mNeutral(double a, double mu, int k, int m) {
  double ret = bico(m, k);
  ret /= zNeutral(mu);
  ret *= (2.0 * (1.0 - a) * (m - k) + 2.0 * a * k + 2.0 * mu) * exp(logGamma(k + mu) + logGamma(-k + m + mu) - logGamma(1.0 + m + 2.0 * mu));
  return ret;
}

double hitchhikingFraction(int a_i, int b_i, double mu, double V) {
  auto a = cast(double)a_i;
  auto b = cast(double)b_i;
  auto hyperg1 = hyperg_1F1(b + mu, 1.0 + 2.0 * mu, (2.0 * a - 1.0) * V);
  auto hyperg2 = hyperg_1F1(b + mu, 2.0 + 2.0 * mu, (2.0 * a - 1.0) * V);
  auto ret = -a * exp(-V) + 1.0 - a;
  ret /= zNeutral(mu);
  ret *= exp(logGamma(b + mu) + logGamma(1.0 - b + mu) - logGamma(1.0 + 2.0 * mu));
  auto term = -2.0 * a * hyperg1;
  term += 2.0 * (1.0 - b + mu) /
     (1.0 + 2.0 * mu) * hyperg2;
  ret *= term;
  auto diff = mNeutral(a, mu, cast(int)b, 1);
  diff -= ret;
  return diff;
}

double mSweep(int a_i, double mu, double V, int k_i, int m_i) {
  auto a = cast(double)a_i;
  auto k = cast(double)k_i;
  auto m = cast(double)m_i;
  double ret = bico(m_i, k_i);
  auto hyperg1 = hyperg_1F1(k + mu, m + 2.0 * mu, (2.0 * a - 1.0) * V);
  auto hyperg2 = hyperg_1F1(k + mu, 1.0 + m + 2.0 * mu, (2.0 * a - 1.0) * V);
  ret /= zNeutral(mu);
  ret *= (-a * exp(-V) + (1.0 - a)) * exp(logGamma(k + mu) + logGamma(-k + m + mu) - logGamma(m + 2.0 * mu)) * ((-2.0 * a) * hyperg1 + 2.0 * (-k + m + mu) / (m + 2.0 * mu) * hyperg2);
  ret += k_i == m_i ? hitchhikingFraction(a_i, 1, mu, V) : 0.0;
  ret += k_i == 0 ? hitchhikingFraction(a_i, 0, mu, V) : 0.0;
  
  return ret;
}

double zSel(int a_i, double mu, double s) {
  auto a = cast(double)a_i;
  auto hyperg = hyperg_1F1(mu, 2.0 * mu, s);
  auto ret = exp(2.0 * logGamma(mu) - logGamma(2.0 * mu)) *
      (1.0 - exp(-(1.0 - a) * s) * hyperg);
  return ret;
}

double mSel(int a_i, double mu, double s, int k_i, int m_i) {
  auto a = cast(double)a_i;
  auto k = cast(double)k_i;
  auto m = cast(double)m_i;
  auto hyperg = hyperg_1F1(k + mu, m + 2.0 * mu, s);
  double ret = bico(m_i, k_i);
  ret /= zSel(a_i, mu, s);
  ret *= exp(logGamma(k + mu) + logGamma(m - k + mu) - logGamma(m + 2.0 * mu));
  ret *= (1.0 - exp(-(1.0 - a) * s) * hyperg);
  return ret;
}

double hitchhikingFractionSel(int a_i, int b_i, double mu, double s, double V) {
  auto a = cast(double)a_i;
  auto b = cast(double)b_i;
  auto hyperg1 = hyperg_1F1(b + mu, 1.0 + 2.0 * mu, -(1.0 - a) * V + a * V);
  auto hyperg2 = hyperg_1F1(b + mu, 1.0 + 2.0 * mu, s - (1.0 - a) * V + a * V);
  auto ret = exp(logGamma(b + mu) + logGamma(1.0 - b + mu) - logGamma(1.0 + 2.0 * mu));
  ret *= exp(-(1.0 - a) * s - a * V);
  ret /= zSel(a_i, mu, s);
  ret *= exp((1.0 - a) * s) * hyperg1 - hyperg2;
  auto diff = mSel(a_i, mu, s, b_i, 1);
  diff -= ret;
  return diff;
}

double mSweepSel(int a_i, double mu, double s, double V, int k_i, int m_i) {
  auto a = cast(double)a_i;
  auto k = cast(double)k_i;
  auto m = cast(double)m_i;
  auto hyperg1 = hyperg_1F1(k + mu, m + 2.0 * mu, -(1.0 - a) * V + a * V);
  auto hyperg2 = hyperg_1F1(k + mu, m + 2.0 * mu, s - (1.0 - a) * V + a * V);
  auto ret = exp(logGamma(k + mu) + logGamma(m - k + mu) - logGamma(m + 2.0 * mu));
  ret *= bico(m_i, k_i);
  ret *= exp(-(1.0 - a) * s - a * V);
  ret /= zSel(a_i, mu, s);
  ret *= exp((1.0 - a) * s) * hyperg1 - hyperg2;
  ret += k_i == m_i ? hitchhikingFractionSel(a_i, 1, mu, s, V) : 0.0;
  ret += k_i == 0 ? hitchhikingFractionSel(a_i, 0, mu, s, V) : 0.0;
  return ret;
}

double mGeneral(int a, double mu, double s, double V, int k, int m) {
  double ret;
  if(fabs(s) < eps && V < eps) ret = mNeutral(a, mu, k, m);
  else if(fabs(s) < eps) ret = mSweep(a, mu, V, k, m);
  else if(fabs(V) < eps) ret = mSel(a, mu, s, k, m);
  else ret = mSweepSel(a, mu, s, V, k, m);
  return ret;
}

unittest {
  auto a = 0;
  auto mu = 0.025;
  auto s = 2.0;
  auto V = 1.0;
  auto m = 5U;
  auto gamma = 0.1;
  
  auto th = [
      [0.9457728047394913,0.02352668668506197,0.011982536075621617,
       0.0080213671250029,0.005991268037810812,0.004705337337012394],
      [0.9523205497232898,0.02004874548100449,0.008706643282440063,
       0.004952065583054215,0.0031305893984905973,0.010841406531721899],
      [0.9366219761117129,0.02198445658489576,0.010782199884747266,
       0.0069989917350783614,0.00510646538848006,0.01850591029508657]];
  
  foreach(k; 0 .. m + 1) {
    assert(approxEqual(mGeneral(a, mu, 0, 0, k, m), th[0][k], 0.0, 1e-12));
    assert(approxEqual(mGeneral(a, mu, 0, V, k, m), th[1][k], 0.0, 1e-12));
    assert(approxEqual(mGeneral(a, mu, s, V, k, m), th[2][k], 0.0, 1e-12));
  }
}

double prop(int aTo, int aFrom, double uP, double uM, double t) {
  double ret = 0.0;
  if(aTo == 0 && aFrom == 0)
    ret = (uM + exp(-t * (uP + uM)) * uP) / (uP + uM);
  if(aTo == 0 && aFrom == 1)
    ret =  (uM - exp(-t * (uP + uM)) * uM) / (uP + uM);
  if(aTo == 1 && aFrom == 0)
    ret = (uP - exp(-t * (uP + uM)) * uP) / (uP + uM);
  if(aTo == 1 && aFrom == 1)
    ret =  (uP + exp(-t * (uP + uM)) * uM) / (uP + uM);
  return ret;
}

double prop_gamma(int aTo, int aFrom, int eTo, int eFrom,
                  double uP, double uM, double gamma, double t)
{

  double ret = 0.0;

  if(aTo == 0 && eTo == 1) {
    if(aFrom == 0 && eFrom == 1) {
      ret = -(-(exp(t*(uM + uP + 2*gamma))*(uM + uP - 2*gamma)*(uM + 
            gamma)) - (uM + uP - 2*gamma)*(uP + gamma) - 
            exp(t*(uM + uP))*(uM - gamma)*(uM + uP + 2*gamma) +         
            exp(2*t*gamma)*(-uP + gamma)*(uM + uP + 2*gamma))/
            (2.*exp(t*(uM + uP + 2*gamma))*(pow(uM + uP,2) -    
            4*pow(gamma,2)));
    }
    if(aFrom == 0 && eFrom == 0) {    
      ret = (2*((-uM + uP)*gamma*cosh(t*gamma)*sinh((t*(uM + uP))/2.) + 
            (uM*(uM + uP) - 2*pow(gamma,2))*cosh((t*(uM + 
            uP))/2.)*sinh(t*gamma)))/
            (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
            4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 1) {
      ret = (2*((uM*(uM + uP) - 2*pow(gamma,2))*cosh(t*gamma)*sinh((t*(uM 
            + uP))/2.) + 
            (-uM + uP)*gamma*cosh((t*(uM + uP))/2.)*sinh(t*gamma)))/
            (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
            4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 0) {
      ret = -(-(exp(t*(uM + uP + 2*gamma))*(uM + uP - 2*gamma)*(uM + 
            gamma)) - (uM + uP - 2*gamma)*(uP + gamma) + 
            exp(t*(uM + uP))*(uM - gamma)*(uM + uP + 2*gamma) -     
            exp(2*t*gamma)*(-uP + gamma)*(uM + uP + 2*gamma))/
            (2.*exp(t*(uM + uP + 2*gamma))*(pow(uM + uP,2) - 
            4*pow(gamma,2)));
    }
  }
         
  if(aTo == 0 && eTo == 0) {
    if(aFrom == 0 && eFrom == 1) {     
      ret = (2*((uM - uP)*gamma*cosh(t*gamma)*sinh((t*(uM + uP))/2.) + 
            (uP*(uM + uP) - 2*pow(gamma,2))*cosh((t*(uM +   
            uP))/2.)*sinh(t*gamma)))/
            (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
            4*pow(gamma,2)));
    }
    if(aFrom == 0 && eFrom == 0) {    
      ret = -(-(((uM + uP - 2*gamma)*(uM + gamma))/exp(t*(uM + uP))) + 
            exp(2*t*gamma)*(uP + gamma)*(-uM - uP + 2*gamma) - 
            ((uM - gamma)*(uM + uP + 2*gamma))/exp(t*(uM + uP - 
            2*gamma)) - (uP - gamma)*(uM + uP + 2*gamma))/
            (2.*exp(2*t*gamma)*(pow(uM + uP,2) - 4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 1) {
      ret = -(-(((uM + uP - 2*gamma)*(uM + gamma))/exp(t*(uM + uP))) + 
            exp(2*t*gamma)*(uP + gamma)*(-uM - uP + 2*gamma) + 
            ((uM - gamma)*(uM + uP + 2*gamma))/exp(t*(uM + uP - 
            2*gamma)) + (uP - gamma)*(uM + uP + 2*gamma))/
            (2.*exp(2*t*gamma)*(pow(uM + uP,2) - 4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 0) {
      ret = (2*((uP*(uM + uP) - 2*pow(gamma,2))*cosh(t*gamma)*sinh((t*(uM     
            + uP))/2.) + 
            (uM - uP)*gamma*cosh((t*(uM + uP))/2.)*sinh(t*gamma)))/
            (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
            4*pow(gamma,2)));
    }
  }
         
  if(aTo == 1 && eTo == 1) {
    if(aFrom == 0 && eFrom == 1) {    
        ret = (2*((uP*(uM + uP) - 2*pow(gamma,2))*cosh(t*gamma)*sinh((t*(uM 
              + uP))/2.) + 
              (uM - uP)*gamma*cosh((t*(uM + uP))/2.)*sinh(t*gamma)))/
              (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
              4*pow(gamma,2)));
    }
    if(aFrom == 0 && eFrom == 0) {
        ret = -(-(((uM + uP - 2*gamma)*(uM + gamma))/exp(t*(uM + uP))) + 
              exp(2*t*gamma)*(uP + gamma)*(-uM - uP + 2*gamma) + 
              ((uM - gamma)*(uM + uP + 2*gamma))/exp(t*(uM + uP - 
              2*gamma)) + (uP - gamma)*(uM + uP + 2*gamma))/
              (2.*exp(2*t*gamma)*(pow(uM + uP,2) - 4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 1) {
        ret = -(-(((uM + uP - 2*gamma)*(uM + gamma))/exp(t*(uM + uP))) + 
              exp(2*t*gamma)*(uP + gamma)*(-uM - uP + 2*gamma) - 
              ((uM - gamma)*(uM + uP + 2*gamma))/exp(t*(uM + uP -     
              2*gamma)) - (uP - gamma)*(uM + uP + 2*gamma))/
              (2.*exp(2*t*gamma)*(pow(uM + uP,2) - 4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 0) {
        ret = (2*((uM - uP)*gamma*cosh(t*gamma)*sinh((t*(uM + uP))/2.) + 
              (uP*(uM + uP) - 2*pow(gamma,2))*cosh((t*(uM + 
              uP))/2.)*sinh(t*gamma)))/
              (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
              4*pow(gamma,2)));
    }
  }

  if(aTo == 1 && eTo == 0) {    
    if(aFrom == 0 && eFrom == 1) {   
        ret = -(-(exp(t*(uM + uP + 2*gamma))*(uM + uP - 2*gamma)*(uM +  
              gamma)) - (uM + uP - 2*gamma)*(uP + gamma) + 
              exp(t*(uM + uP))*(uM - gamma)*(uM + uP + 2*gamma) - 
              exp(2*t*gamma)*(-uP + gamma)*(uM + uP + 2*gamma))/
              (2.*exp(t*(uM + uP + 2*gamma))*(pow(uM + uP,2) - 
              4*pow(gamma,2)));
    }
    if(aFrom == 0 && eFrom == 0) { 
        ret = (2*((uM*(uM + uP) - 2*pow(gamma,2))*cosh(t*gamma)*sinh((t*(uM     
              + uP))/2.) + 
              (-uM + uP)*gamma*cosh((t*(uM + uP))/2.)*sinh(t*gamma)))/
              (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
              4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 1) {
        ret = (2*((-uM + uP)*gamma*cosh(t*gamma)*sinh((t*(uM + uP))/2.) + 
              (uM*(uM + uP) - 2*pow(gamma,2))*cosh((t*(uM + 
              uP))/2.)*sinh(t*gamma)))/
              (exp((t*(uM + uP + 2*gamma))/2.)*(pow(uM + uP,2) - 
              4*pow(gamma,2)));
    }
    if(aFrom == 1 && eFrom == 0) {
        ret = -(-(exp(t*(uM + uP + 2*gamma))*(uM + uP - 2*gamma)*(uM +  
              gamma)) - (uM + uP - 2*gamma)*(uP + gamma) - 
              exp(t*(uM + uP))*(uM - gamma)*(uM + uP + 2*gamma) +     
              exp(2*t*gamma)*(-uP + gamma)*(uM + uP + 2*gamma))/
              (2.*exp(t*(uM + uP + 2*gamma))*(pow(uM + uP,2) - 
              4*pow(gamma,2)));
    }
  }
  return ret;
}

double g(int a1, int a2, double uP, double uM, double t) {
  double ret = 0.0;
  foreach(a; 0 .. 2)
    ret += lambda(a, 1, uP, uM, 0.0) *
           prop(a1, a, uP, uM, t) * prop(a2, a, uP, uM, t);
  return ret;
}

double g_gamma(int a1, int e1, int a2, int e2, double uP, double uM,
               double gamma, double t)
{
  double ret = 0.0;
  foreach(a; 0 .. 2) foreach(e; 0 .. 2)
    ret += 0.5 * lambda(a, e, uP, uM, gamma) *
           prop_gamma(a1, a, e1, e, uP, uM, gamma, t) * 
           prop_gamma(a2, a, e2, e, uP, uM, gamma, t);
  return ret;
}

double Q(int k1, int m1, int k2, int m2, double t, double mu, double s, 
         double gamma, double V, double uP, double uM)
{
  double ret = 0.0;
  if(fabs(gamma) < 1.0e-12) {
    foreach(a1; 0 .. 2) {
      auto v1 = mGeneral(a1, mu, s, V, k1, m1);
      foreach(a2; 0 .. 2) {
        double val = g(a1, a2, uP, uM, t);
        val *= v1;
        double m = mGeneral(a2, mu, s, V, k2, m2);
        val *= m;
        ret += val;
      }
    }
  }
  else {
    foreach(a1; 0 .. 2) foreach(e1_i; 0 .. 2) {
      auto e1 = cast(double)e1_i;
      auto v1 = mGeneral(a1, mu, e1 * s - (1.0 - e1) * s, V, k1, m1);
      foreach(a2; 0 .. 2) foreach(e2; 0 .. 2) {
        double val = g_gamma(a1, e1_i, a2, e2, uP, uM, gamma, t);
        val *= v1;
        double m = mGeneral(a2, mu, e2 * s - (1.0 - e2) * s, V, k2, m2);
        val *= m;
        ret += val;
      }
    }
  }
  return ret;                        
}

double Q(int k1, int m1, int k2, int m2, double t, double mu, double s, 
         double gamma, double V)
{
  double uP = uLinkage(s, mu, V);
  double uM = uLinkage(-s, mu, V);
  return Q(k1, m1, k2, m2, t, mu, s, gamma, V, uP, uM);
}

double Q1vsM(int k, int m, double t, double mu, double s, double gamma, 
             double V)
{
  double ret = Q(0, 1, k, m, t, mu, s, gamma, V);
  ret += Q(1, 1, m - k, m, t, mu, s, gamma, V);
  return ret;
}

unittest {
  auto a = 0;
  auto mu = 0.025;
  auto s = 2.0;
  auto V = 1.0;
  auto m = 5U;
  auto gamma = 0.1;
  auto t = 4.0;
  
  auto th = [
      [0.7756278272629531,0.020356283378045733,0.011266356816185446,
       0.008737546384439078,0.009161671344827057,0.17485031481355084],
      [0.782101141273109,0.01698994323216128,0.008027815716697952,
       0.005630893148796326,0.006189391647333807,0.181060814981903],
      [0.8232228035723034,0.015933218904722574,0.006787423253861372,
       0.004412407973178893,0.004737710155317756,0.1449064361406177],
      [0.7756278272629531,0.020356283378045733,0.011266356816185446,
       0.008737546384439078,0.009161671344827057,0.17485031481355084],
      [0.782101141273109,0.01698994323216128,0.008027815716697952,
       0.005630893148796326,0.006189391647333807,0.181060814981903],
      [0.8232228035723034,0.015933218904722574,0.006787423253861372,
       0.004412407973178893,0.004737710155317756,0.1449064361406177]];
  
  foreach(k; 0 .. m  + 1) {
    assert(approxEqual(Q1vsM(k, m, t, mu, 0, 0, 0), th[0][k], 0.0, 1e-12));
    auto v = Q1vsM(k, m, t, mu, 0, 0, V);
    assert(approxEqual(v, th[1][k], 0.0, 1e-12), text([v, th[1][k]]));
    assert(approxEqual(Q1vsM(k, m, t, mu, s, 0, V), th[2][k], 0.0, 1e-12));
    assert(approxEqual(Q1vsM(k, m, t, mu, 0, gamma, 0), th[3][k], 0.0, 1e-12));
    assert(approxEqual(Q1vsM(k, m, t, mu, 0, gamma, V), th[4][k], 0.0, 1e-12));
    assert(approxEqual(Q1vsM(k, m, t, mu, s, gamma, V), th[5][k], 0.0, 1e-12));
  }
}

double[] Q(int k1, int m1, int m2, double t, double mu, double s, 
         double gamma, double V)
{
  auto ret = new double[m2 + 1];
  double uP = uLinkage(s, mu, V);
  double uM = uLinkage(-s, mu, V);
  foreach(k; 0 .. m2 + 1) {
    ret[k] = Q(k1, m1, k, m2, t, mu, s, gamma, V, uP, uM);
  }
  return ret;
}


double[] Q1vsM(int m, double t, double mu, double s, double gamma, double V) {
  auto q1 = Q(0, 1, m, t, mu, s, gamma, V);
  auto q2 = Q(1, 1, m, t, mu, s, gamma, V);
  foreach(k; 0 .. m + 1) {
    q1[k] += q2[m - k];
  }
  return q1;
}

class IllegalParametersException : Exception {
  this(string msg) {
    super(msg);
  }
}

class SingleSpectrumProb {
  int m;
  double[5] lastP;
  double[] last_prob;
  bool first;
  
  this(int m) {
    this.m = m;
  }

  immutable(double)[] opCall(double t, double mu, double s, double gamma, 
                             double V)
  {
    if(lastP != [t, mu, s, gamma, V]) {
      last_prob = prob(t, mu, s, gamma, V);
      lastP = [t, mu, s, gamma, V];
    }
    return last_prob.idup;
  }
  
  double[] prob(double t, double mu, double s, double gamma, double V) {
    if(t < 0) {
      throw new IllegalParametersException("t negative");
    }
    if(mu < 0) {
      throw new IllegalParametersException("mu negative");
    }
    if(s < 0) {
      throw new IllegalParametersException("s negative");
    }
    if(gamma < 0) {
      throw new IllegalParametersException("gamma negative");
    }
    if(V < 0) {
      throw new IllegalParametersException("V negative");
    }

    return Q1vsM(m, t, mu, s, gamma, V);
  }
}
