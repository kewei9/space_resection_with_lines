
#include "StdAfx.h"
#include "powell.h"

#include <math.h>
#include <stdlib.h>

//#include <fstream>

#define NR_END		1
#define FREE_ARG	char*

static double sqrarg;
#define SQR(a)	((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

static double dsqrarg;
#define DSQR(a)	((dsqrarg = (a)) == 0.0 ? 0.0 : dsqrarg * dsqrarg)

static double dmaxarg1, dmaxarg2;
#define DMAX(a, b) (dmaxarg1 = (a), dmaxarg2 = (b), (dmaxarg1) > (dmaxarg2) ? \
(dmaxarg1) : (dmaxarg2))

static double dminarg1, dminarg2;
#define DMIN(a, b) (dminarg1 = (a), dminarg2 = (b), (dminarg1) < (dminarg2) ? \
(dminarg1) : (dminarg2))

static double maxarg1, maxarg2;
#define FMAX(a, b) (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? \
(maxarg1) : (maxarg2))

static double minarg1, minarg2;
#define FMIN(a, b) (minarg1 = (a), minarg2 = (b), (minarg1) < (minarg2) ? \
(minarg1) : (minarg2))

static int imaxarg1, imaxarg2;
#define IMAX(a, b) (imaxarg1 = (a), imaxarg2 = (b), (imaxarg1) > (imaxarg2) ? \
(imaxarg1) : (imaxarg2))

static int iminarg1, iminarg2;
#define IMIN(a, b) (iminarg1 = (a), iminarg2 = (b), (iminarg1) < (iminarg2) ? \
(iminarg1) : (iminarg2))

#define SIGN(a, b)	((b) >= 0.0 ? fabs(a) : -fabs(a))
static double swap;
#define SWAP(a, b)  {swap = (a); (a) = (b); (b) = swap;}

#define ITMAX		100
#define CGOLD		0.3819660
#define ZEPS		1.0e-10

#define SHFT(a, b, c, d) (a) = (b); (b) = (c); (c) = (d);

//////////////////////////////////////////////////////////////////////////

class engine_handle
{
public:
	static engine* get_engine_instance();
	static void release_engine();
	static void register_engine(engine* energy_engine_);

private:
	static engine* energy_engine;
};

engine* engine_handle::energy_engine = 0;

engine* engine_handle::get_engine_instance()
{
	return energy_engine;
}

void engine_handle::register_engine(engine* energy_engine_)
{
	energy_engine = energy_engine_;
}

void engine_handle::release_engine()
{
	energy_engine = 0;
}

double compute_energy(double* value)
{
	engine* energy_engine = engine_handle::get_engine_instance();
	return energy_engine->compute_energy(value);
}

engine::engine()
{

}

engine::~engine()
{

}

void engine::register_engine()
{
	engine_handle::register_engine(this);
}

void engine::release_engine()
{
	engine_handle::release_engine();
}

optimization::optimization()
{
}

optimization::~optimization()
{
}

optimization_powell::optimization_powell()
{

}

optimization_powell::~optimization_powell()
{

}

double **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, double (*func)(double[]));

double optimization_powell::optimize(double* value, int value_count, double ftol /* = 0.000001 */)
{
	double **xi = matrix(1, value_count, 1, value_count);
	for(int i = 1; i <= value_count; i++){
		for(int j = 1; j <= value_count; j++)
			xi[i][j] = 0.0;
		xi[i][i] = 1.0; 
	}
	int iter;
	double fret;
	powell(value, xi, value_count, (float)ftol, &iter, &fret, compute_energy);
	free_matrix(xi, 1, value_count, 1, value_count);

	return fret;
}

//////////////////////////////////////////////////////////////////////////

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
{
	double ulim, u, r, q, fu, dum;

	*fa = (*func)(*ax);
	*fb = (*func)(*bx);
	if(*fb > *fa){
		SHFT(dum, *ax, *bx, dum);
		SHFT(dum, *fb, *fa, dum);
	}
	*cx = (*bx) + 1.618034 * (*bx - *ax);
	*fc = (*func)(*cx);
	while(*fb > *fc){
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * SIGN(FMAX(fabs(q - r), 1.0e-20), q - r));
		ulim = (*bx) + 100.0 * (*cx - *bx);
		if((*bx - u ) * ( u - *cx) > 0.0){
			fu = (*func)(u);
			if(fu < *fc){
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}else if(fu > *fb){
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + 1.618034 * (*cx - *bx);
			fu = (*func)(u);
		}else if((*cx - u) * (u - ulim) > 0.0){
			fu = (*func)(u);
			if(fu < *fc){
				SHFT(*bx, *cx, u, *cx + 1.618034 * (*cx - *bx))
				SHFT(*fb, *fc, fu, (*func)(u))
			}
		}else if ((u - ulim) * (ulim - *cx) >= 0.0){
			u = ulim;
			fu = (*func)(u);
		}else{
			u = (*cx) + 1.618034 * (*cx - *bx);
			fu = (*func)(u);
		}
		SHFT(*ax, *bx, *cx, u)
		SHFT(*fa, *fb, *fc, fu)
	}
}

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
{
	int iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*f)(x);
	for(iter = 1; iter <= ITMAX; iter++){
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if(fabs(x - xm) <= (tol2 - 0.5 * (b - a))){
			*xmin = x;
			return fx;
		}
		if(fabs(e) > tol1){
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if(q > 0.0)
				p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if(fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else{
				d = p / q;
				u = x + d;
				if(u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		}else{
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);
		if(fu <= fx){
			if( u >= x)
				a = x;
			else
				b = x;
			SHFT(v, w, x, u);
			SHFT(fv, fw, fx, fu);
		}else{
			if(u < x)
				a = u; 
			else
				b = u;
			if(fu <= fw || w == x){
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}else if(fu <= fv || v == x || v == w){
				v = u;
				fv = fu;
			}
		}
	}
	*xmin = x;
	return fx;
}

void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, double (*func)(double[]))
{
//	ofstream output("D:\\hbw_serial_values.txt", ios::out);

	void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));
	int i, ibig, j;
	double del, fp, fptt, t, pt[51], ptt[51], xit[51];
	
	*fret = (*func)(p);
	for(j = 1; j <= n; j++)
		pt[j] = p[j];
	for(*iter = 1; *iter <= ITMAX; ++(*iter)){
		fp =(*fret);
		ibig = 0;
		del = 0.0;
		for(i = 1; i <= n; i++){
			for(j = 1; j <= n; j++)
				xit[j] = xi[j][i];
			fptt = (*fret);
			linmin(p, xit, n, fret, func);

//			output<<*fret<<endl;

			if(fptt - (*fret) > del){
				del = fptt - (*fret);
				ibig = i;
			}
		}
		if(2.0 * (fp - (*fret)) <= ftol * (fabs(fp) + fabs(*fret)) + 1.0e-25)
		{
//			output.close();
			return;
		}

		for(j = 1; j <= n; j++){
			ptt[j] = 2.0 * p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		fptt = (*func)(ptt);
		if(fptt < fp){
			t = 2.0 * (fp - 2.0 * (*fret) + fptt) * SQR(fp - (*fret) - del) - del * SQR(fp - fptt);
			if(t < 0.0){
				linmin(p, xit, n, fret, func);
				for(j = 1; j <= n; j++){
					xi[j][ibig] = xi[j][n];
					xi[j][n] = xit[j];
				}
			}
		}
	}

//	output.close();
}

int ncom;
double pcom[51], xicom[51], (*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	double f1dim(double x);
	
	int j;
	double xx, xmin, fx, fb, fa, bx, ax;
	
	ncom = n;
	nrfunc = func;
	for(j = 1; j <= n; j++){
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0;
	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
	*fret = brent(ax, xx, bx, f1dim, 2.0e-4, &xmin);
	for(j = 1; j <= n; j++){
		xi[j] *= xmin;
		p[j] += xi[j];
	}
}

double f1dim(double x)
{
	int j;
	double f, xt[51];
	
	for(j = 1; j <= ncom; j++)
		xt[j] = pcom[j] + x * xicom[j];
	f = (*nrfunc)(xt);
	return f;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;
	
	m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
	if(!m)
		return NULL;
	
	m += NR_END;
	m -= nrl;
	
	m[nrl] = (double*) malloc ((size_t)((nrow * ncol + NR_END) * sizeof(double)));
	if(!m[nrl])
		return NULL;
	
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;
	
	return m;
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}