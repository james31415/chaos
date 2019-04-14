#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static float maxarg1, maxarg2;
#define FMAX(a, b) (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

#define NR_END 1
#define FREE_ARG char*

void
nrerror(char error_text[]) {
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}

float*
vector(long nl, long nh) {
	float *v;
	v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
	return v - nl + NR_END;
}

void
free_vector(float* v, long nl, long nh) {
	free((FREE_ARG)(v + nl - NR_END));
}

float**
matrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;
	m = (float**)malloc((size_t)((nrow + NR_END) * sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	m[nrl] = (float*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl+1; i <= nrh; ++i) m[i] = m[i - 1] + ncol;

	return m;
}

void
free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void
rk4(float y[], float dydx[], int n, float x, float h, float yout[], void (*derivs)(float, float [], float [])) {
	int i;
	float xh, hh, h6, *dym, *dyt, *yt;

	dym = vector(1, n);
	dyt = vector(1, n);
	yt = vector(1, n);
	hh = h * 0.5f;
	h6 = h / 6.0f;
	xh = x + hh;
	for (i = 1; i <= n; ++i) yt[i] = y[i] + hh * dydx[i];
	(*derivs)(xh, yt, dyt);
	for (i = 1; i <= n; ++i) yt[i] = y[i] + hh * dyt[i];
	(*derivs)(xh, yt, dym);
	for (i = 1; i <= n; ++i) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x + h, yt, dyt);
	for (i = 1; i <= n; ++i) yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0f * dym[i]);

	free_vector(yt, 1, n);
	free_vector(dyt, 1, n);
	free_vector(dym, 1, n);
}

void
rkdumb(float *xx, float **y, float vstart[], int nvar, float x1, float x2, int nstep, void (*derivs)(float, float[], float[])) {
	int i, k;
	float x, h;
	float *v, *vout, *dv;

	v = vector(1, nvar);
	vout = vector(1, nvar);
	dv = vector(1, nvar);
	for (i = 1; i <= nvar; ++i) {
		v[i] = vstart[i];
		y[i][1] = v[i];
	}
	xx[1] = x1;
	x = x1;
	h = (x2 - x1) / nstep;
	for (k = 1; k <= nstep; ++k) {
		(*derivs)(x, v, dv);
		rk4(v, dv, nvar, x, h, vout, derivs);
		if ((float)(x + h) == x) nrerror("Step size too small in routine rkdumb");
		x += h;
		xx[k+1] = x;
		for (i = 1; i <= nvar; ++i) {
			v[i] = vout[i];
			y[i][k+1] = v[i];
		}
	}
	free_vector(dv,1,nvar);
	free_vector(vout,1,nvar);
	free_vector(v,1,nvar);
}

void
rkck(float y[], float dydx[], int n, float x, float h, float yout[], float yerr[], void (*derivs)(float, float[], float[])) {
	int i;
	static float a2  = 0.2f, a3 = 0.3f, a4 = 0.6f, a5 = 1.0f, a6 = 0.875f, b21 = 0.2f,
		b31 = 3.0f / 40.0f, b32 = 9.0f / 40.0f, b41 = 0.3f, b42 = -0.9f, b43 = 1.2f,
		b51 = -11.0f / 54.0f, b52 = 2.5f, b53 = -70.0f / 27.0f, b54 = 35.0f / 27.0f,
		b61 = 1631.0f / 55296.0f, b62 = 175.0f / 512.0f, b63 = 575.0f / 13824.0f,
		b64 = 44275.0f / 110592.0f, b65 = 253.0f / 4096.0f, c1 = 37.0f / 378.0f,
		c3 = 250.0f / 621.0f, c4 = 125.0f / 594.0f, c6 = 512.0f / 1771.0f,
		dc5 = -277.0f / 14336.0f;
	float dc1 = c1 - 2825.0f / 27648.0f, dc3 = c3 - 18575.0f / 48384.0f,
		dc4 = c4 - 13525.0f / 55296.0f, dc6 = c6 - 0.25f;
	float *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;

	ak2 = vector(1, n);
	ak3 = vector(1, n);
	ak4 = vector(1, n);
	ak5 = vector(1, n);
	ak6 = vector(1, n);
	ytemp = vector(1, n);
	for (i = 1; i <= n; ++i) ytemp[i] = y[i] + b21 * h * dydx[i];
	(*derivs)(x + a2 * h, ytemp, ak2);
	for (i = 1; i <= n; ++i) ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
	(*derivs)(x + a3 * h, ytemp, ak3);
	for (i = 1; i <= n; ++i) ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
	(*derivs)(x + a4 * h, ytemp, ak4);
	for (i = 1; i <= n; ++i) ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
	(*derivs)(x + a5 * h, ytemp, ak5);
	for (i = 1; i <= n; ++i) ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
	(*derivs)(x + a6 * h, ytemp, ak6);
	for (i = 1; i <= n; ++i) yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
	for (i = 1; i <= n; ++i) yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);

	free_vector(ytemp, 1, n);
	free_vector(ak6, 1, n);
	free_vector(ak5, 1, n);
	free_vector(ak4, 1, n);
	free_vector(ak3, 1, n);
	free_vector(ak2, 1, n);
}

#define SAFETY 0.9f
#define PGROW -0.2f
#define PSHRNK -0.25f
#define ERRCON 1.89e-4
void rkqs(float y[], float dydx[], int n, float *x, float htry, float eps, float yscal[], float *hdid, float *hnext, void (*derivs)(float, float[], float[])) {
	int i;
	float errmax, h, xnew, *yerr, *ytemp;

	yerr = vector(1, n);
	ytemp = vector(1, n);
	h = htry;
	for (;;) {
		rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);
		errmax = 0.0f;
		for (i = 1; i <= n; ++i) errmax = FMAX(errmax, (float)fabs(yerr[i] / yscal[i]));
		errmax /= eps;
		if (errmax > 1.0f) {
			h = SAFETY * h * (float)pow(errmax, PSHRNK);
			if (h < 0.1f * h) h *= 0.1f;
			xnew = (*x) + h;
			if (xnew == *x) nrerror("stepsize underflow in rkqs");
			continue;
		} else {
			if (errmax > ERRCON) *hnext = SAFETY * h * (float)pow(errmax, PGROW);
			else *hnext = 5.0f * h;
			*x += (*hdid = h);
			for (i = 1; i <= n; ++i) y[i] = ytemp[i];
			break;
		}
	}

	free_vector(ytemp, 1, n);
	free_vector(yerr, 1, n);
}

#define MAXSTP 10000
#define TINY 1.0e-30f

int kmax, kount;
float *xp, **yp, dxsav;

void
odeint(float ystart[], int nvar, float x1, float x2, float eps, float h1,
	float hmin, int* nok, int* nbad,
	void (*derivs)(float, float[], float[]),
	void (*rkqs)(float[], float[], int, float*, float, float, float[], float*, float*, void (*)(float, float[], float[]))) {
	int nstp, i;
	float xsav, x, hnext, hdid, h;
	float *yscal, *y, *dydx;

	yscal = vector(1, nvar);
	y = vector(1, nvar);
	dydx = vector(1, nvar);
	x = x1;
	h = (x2 > x1) ? (float)fabs(h1) : -(float)fabs(h1);
	*nok = (*nbad) = kount = 0;
	for (i = 1; i <= nvar; ++i) y[i] = ystart[i];
	xsav = (kmax > 0) ? x - dxsav * 2.0f : 0;
	for (nstp = 1; nstp <= MAXSTP; ++nstp) {
		(*derivs)(x, y, dydx);
		for (i = 1; i <= nvar; ++i) yscal[i] = (float)fabs(y[i]) + (float)fabs(dydx[i] * h) + TINY;
		if (kmax > 0 && kount < kmax - 1 && fabs(x - xsav) > (float)fabs(dxsav)) {
			xp[++kount] = x;
			for (i = 1; i <= nvar; ++i) yp[i][kount] = y[i];
			xsav = x;
		}
		if ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;
		(*rkqs)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x - x2) * (x2 - x1) >= 0.0) {
			for (i = 1; i <= nvar; ++i) ystart[i] = y[i];
			if (kmax) {
				xp[++kount] = x;
				for (i = 1; i <= nvar; ++i) yp[i][kount] = y[i];
			}

			free_vector(dydx, 1, nvar);
			free_vector(y, 1, nvar);
			free_vector(yscal, 1, nvar);

			return;
		}
		if ((float)fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h = hnext;
	}
	nrerror("Too many steps in routine odeint");
}

void
derivative(float, float y[], float dydx[]) {
	dydx[1] = y[1];
}

int
main(void) {
	kmax = 10000;
	dxsav = 0.05f;
	xp = vector(1, kmax);
	yp = matrix(1, 1, 1, kmax);

	float* y0;
	y0 = vector(1, 1);
	y0[1] = 1;

	int nok, nbad;
	odeint(y0, 1, 0, 5, 1e-8f, 0.05f, 0, &nok, &nbad, derivative, rkqs);

	for (int i = 1; i < kount; ++i) {
		printf("%f: %f\n", xp[i], yp[1][i]);
	}
	printf("nOK: %d\n", nok);
	printf("nBad: %d\n", nbad);

	free_matrix(yp, 1, 1, 1, kmax);
	free_vector(xp, 1, kmax);

	return 0;
}

