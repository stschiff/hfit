import std.math;
import std.random;
import std.range;
import std.mathspecial;
import std.exception;
import std.stdio;
import std.conv;
import std.algorithm;

class MCMC(MinfuncType) {
    MinfuncType minfunc;
    int N;
    double[][] trace;
    double[] fvalues;
    double[] currentPoint;
    double[] successRate;
    double[] stepWidths;
    double f;
    double fnext;
    
    this(MinfuncType minfunc, double[] init) {
        this.minfunc = minfunc;
        N = minfunc.nrParams;
        enforce(init.length == N);
        successRate = new double[N];
        successRate[] = 0.44;
        currentPoint = init.dup;
        f = minfunc(currentPoint);
        stepWidths = init.map!"a*0.01"().array();
    }
    
    void run(int mainCycles, int maxCycles=100000) {
        auto cycle = 0;
        while(cycle < maxCycles) {
            cycle += 1;
            foreach(i; iota(N).randomCover()) {
                auto prop = proposal(i);
                currentPoint[i] += prop;
                auto fnext = minfunc(currentPoint);
                if(accept(fnext, f)) {
                    f = fnext;
                    successRate[i] = 0.99 * successRate[i] + 0.01;
                }
                else {
                    currentPoint[i] -= prop;
                    successRate[i] *= 0.99;
                }
            }
            if(cycle == 10 || cycle == 20 || cycle == 50 || cycle % 100 == 99) {
                stderr.writefln("cycle %s, Success rates: %s, adapting step widths", cycle, successRate);
                stderr.writefln("Current point: %s", currentPoint);
                stderr.writefln("Current function value: %s", f);
                adaptStepWidths();
                if(trace.length - getNrBurninCycles() > mainCycles) {
                    return;
                }
            }
            trace ~= currentPoint.dup;
            fvalues ~= f;
        }
    }
    
    int getNrBurninCycles() {
        auto scoreBound = fvalues.minCount()[0] + 10.0;
        foreach(int i, val; fvalues) {
            if(val<scoreBound)
                return i;
        }
        return fvalues.length.to!int() - 1;
    }
    
    bool accept(double fnext, double f) {
        if(fnext < f)
            return true;
        auto rnd = uniform(0.0, 1.0);
        if(rnd < exp(f - fnext))
            return true;
        return false; 
    }
    
    double proposal(int i) {
        auto rnd = uniform(0.0, 1.0);
        return normalDistributionInverse(rnd) * stepWidths[i];
    }
    
    void adaptStepWidths() {
        stderr.writeln("Sucess rates:, ", successRate);
        foreach(i; 0 .. N) {
            if(successRate[i] < 0.29)
                stepWidths[i] /= 1.5;
            if(successRate[i] > 0.59)
                stepWidths[i] *= 1.5;
        }
        stderr.writeln("new step widths: ", stepWidths);
    }

    void reportTraces(string tracefile) {
        auto f = File(tracefile, "w");
        f.write(minfunc.paramNames().join("\t"));
        f.writeln("\tscore");
        foreach(i, point; trace) {
            f.write(point.map!"text(a)"().join("\t"));
            f.writefln("\t%20.2f", fvalues[i]);
        }
        f.close();
    }
}
