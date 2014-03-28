import std.algorithm;
import std.math;

T SIGN(T)(T a, T b)	{
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

class Bracketmethod {
	double ax,bx,cx,fa,fb,fc;

	void bracket(T)(double a, double b, T func) {
		immutable double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
		ax=a; bx=b;
		double fu;
		fa=func(ax);
		fb=func(bx);
		if (fb > fa) {
      double dummy = ax;
      ax = bx;
      bx = dummy;
      dummy = fb;
      fb = fa;
      fa = dummy;
		}
		cx=bx+GOLD*(bx-ax);
		fc=func(cx);
		while(fb > fc) {
			double r=(bx-ax)*(fb-fc);
			double q=(bx-cx)*(fb-fa);
			double u=bx-((bx-cx)*q-(bx-ax)*r)/
				(2.0*SIGN(max(abs(q-r),TINY),q-r));
			double ulim=bx+GLIMIT*(cx-bx);
			if ((bx-u)*(u-cx) > 0.0) {
				fu=func(u);
				if (fu < fc) {
					ax=bx;
					bx=u;
					fa=fb;
					fb=fu;
					return;
				} else if (fu > fb) {
					cx=u;
					fc=fu;
					return;
				}
				u=cx+GOLD*(cx-bx);
				fu=func(u);
			} else if ((cx-u)*(u-ulim) > 0.0) {
				fu=func(u);
				if (fu < fc) {
					shft3(bx,cx,u,u+GOLD*(u-cx));
					shft3(fb,fc,fu,func(u));
				}
			} else if ((u-ulim)*(ulim-cx) >= 0.0) {
				u=ulim;
				fu=func(u);
			} else {
				u=cx+GOLD*(cx-bx);
				fu=func(u);
			}
			shft3(ax,bx,cx,u);
			shft3(fa,fb,fc,fu);
		}
	}

	void shft2(ref double a, ref double b, double c) {
		a=b;
		b=c;
	}
  
	void shft3(ref double a, ref double b, ref double c, double d) {
		a=b;
		b=c;
		c=d;
	}
  
	void mov3(ref double a, ref double b, ref double c, double d, double e, double f)	{
		a=d; b=e; c=f;
	}
};

class Brent : Bracketmethod {
	double xmin,fmin;
	double tol;

	this(double tol=3.0e-8) {
    this.tol = tol;
  }
  
	double minimize(T)(T func)	{
		immutable int ITMAX=100;
		immutable double CGOLD=0.3819660;
		immutable double ZEPS=double.epsilon*1.0e-3;
		double a,b,d=0.0,etemp,fu,fv,fw,fx;
		double p,q,r,tol1,tol2,u,v,w,x,xm;
		double e=0.0;
	
		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=func(x);
		for (int iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*abs(x)+ZEPS);
			if (abs(x-xm) <= (tol2-0.5*(b-a))) {
				fmin=fx;
				return xmin=x;
			}
			if (abs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=abs(q);
				etemp=e;
				e=d;
				if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x)
						|| p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(abs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=func(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				shft3(v,w,x,u);
				shft3(fv,fw,fx,fu);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		throw new Exception("Too many iterations in brent");
	}
};
