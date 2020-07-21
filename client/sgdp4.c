/* > sgdp4.c
 *
 * 1.00 around 1980 - Felix R. Hoots & Ronald L. Roehrich, from original
 *                    SDP4.FOR and SGP4.FOR
 *
 ************************************************************************
 *
 *     Made famous by the spacetrack report No.3:
 *     "Models for Propogation of NORAD Element Sets"
 *     Edited and subsequently distributed by Dr. T. S. Kelso.
 *
 ************************************************************************
 *
 *     This conversion by:
 *     Paul S. Crawford and Andrew R. Brooks
 *     Dundee University
 *
 *                          NOTE !
 *  This code is supplied "as is" and without warranty of any sort.
 *
 * (c) 1994-2004, Paul Crawford, Andrew Brooks
 *
 ************************************************************************
 *
 * 1.07 arb     Oct    1994 - Transcribed by arb Oct 1994 into 'C', then
 *                            modified to fit Dundee systems by psc.
 *
 * 1.08 psc Mon Nov  7 1994 - replaced original satpos.c with SGP4 model.
 *
 * 1.09 psc Wed Nov  9 1994 - Corrected a few minor translation errors after
 *                            testing with example two-line elements.
 *
 * 1.10 psc Mon Nov 21 1994 - A few optimising tweeks.
 *
 * 1.11 psc Wed Nov 30 1994 - No longer uses eloset() and minor error in the
 *                            SGP4 code corrected.
 *
 * 2.00 psc Tue Dec 13 1994 - arb discovered the archive.afit.af.mil FTP site
 *                            with the original FORTRAN code in machine form.
 *                            Tidied up and added support for the SDP4 model.
 *
 * 2.01 psc Fri Dec 23 1994 - Tested out the combined SGP4/SDP4 code against
 *                            the original FORTRAN versions.
 *
 * 2.02 psc Mon Jan 02 1995 - Few more tweeks and tidied up the
 *                            doccumentation for more general use.
 *
 * 3.00 psc Mon May 29 1995 - Cleaned up for general use & distrabution (to
 *                            remove Dundee specific features).
 *
 * 3.01 psc Mon Jan 12 2004 - Minor bug fix for day calculation.
 *
 * 3.02 psc Mon Jul 10 2006 - Added if(rk < (real)1.0) test for sub-orbital decay.
 *
 * 3.03 psc Sat Aug 05 2006 - Added trap for divide-by-zero when calculating xlcof.
 *
 */

static const char SCCSid[] = "@(#)sgdp4.c       3.03 (C) 1995 psc SatLib: Orbital Model";

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/* ================ single / double precision fix-ups =============== */

#include "sgdp4h.h"

#define ECC_ZERO		((real)0.0)		/* Zero eccentricity case ? */
#define ECC_ALL			((real)1.0e-4)	/* For all drag terms in GSFC case. */
#define ECC_EPS			((real)1.0e-6)	/* Too low for computing further drops. */
#define ECC_LIMIT_LOW	((real)-1.0e-3)	/* Exit point for serious decaying of orbits. */
#define ECC_LIMIT_HIGH	((real)(1.0 - ECC_EPS))	/* Too close to 1 */

#define EPS_COSIO		(1.5e-12)	/* Minimum divisor allowed for (...)/(1+cos(IO)) */

#define TOTHRD  (2.0/3.0)

#if defined( SGDP4_SNGL ) || 0
#define NR_EPS  ((real)(1.0e-6))    /* Minimum ~1e-6 min for float. */
#else
#define NR_EPS  ((real)(1.0e-12))    /* Minimum ~1e-14 for double. */
//#define NR_EPS  ((real)(1.0e-14))    /* Minimum ~1e-14 for double. */
//#define NR_EPS  ((real)(1.0e-8))    /* Minimum ~1e-14 for double. */
#endif

#define Q0      ((real)120.0)
#define S0      ((real)78.0)
#define XJ2     ((real)1.082616e-3)
#define XJ3     ((real)-2.53881e-6)
#define XJ4     ((real)-1.65597e-6)
#define XKMPER  (6378.135)            /* Km per earth radii */
#define XMNPDA  (1440.0)              /* Minutes per day */
#define AE      (1.0)                 /* Earth radius in "chosen units". */

#if 0
/* Original code constants. */
#define XKE     (0.743669161e-1)
#define CK2     ((real)5.413080e-4)   /* (0.5 * XJ2 * AE * AE) */
#define CK4     ((real)0.62098875e-6) /* (-0.375 * XJ4 * AE * AE * AE * AE) */
#define QOMS2T  ((real)1.88027916e-9) /* (pow((Q0 - S0)*AE/XKMPER, 4.0)) */
#define KS      ((real)1.01222928)    /* (AE * (1.0 + S0/XKMPER)) */
#else
/* GSFC improved coeficient resolution. */
#define XKE     ((real)7.43669161331734132e-2)
#define CK2     ((real)(0.5 * XJ2 * AE * AE))
#define CK4     ((real)(-0.375 * XJ4 * AE * AE * AE * AE))
#define QOMS2T  ((real)1.880279159015270643865e-9) /* (pow((Q0 - S0)*AE/XKMPER, 4.0)) */
#define KS      ((real)(AE * (1.0 + S0/XKMPER)))
#endif
static const real a3ovk2 = (real)(-XJ3 / CK2 * (AE * AE * AE));
void fatal_error(const char *format, ...);
/* ================= Copy of the orbital elements ==================== */

static double xno;  /* Mean motion (rad/min) */
static real xmo;    /* Mean "mean anomaly" at epoch (rad). */
static real eo;     /* Eccentricity. */
static real xincl;  /* Equatorial inclination (rad). */
static real omegao; /* Mean argument of perigee at epoch (rad). */
static real xnodeo; /* Mean longitude of ascending node (rad, east). */
static real bstar;  /* Drag term. */

double SGDP4_jd0;  /* Julian Day for epoch (available to outside functions. */

/* ================== Local "global" variables for SGP4 ================= */

static int imode = SGDP4_NOT_INIT;
static real sinIO, cosIO, sinXMO, cosXMO;
static real c1, c2, c3, c4, c5, d2, d3, d4;
static real omgcof, xmcof, xlcof, aycof;
static real t2cof, t3cof, t4cof, t5cof;
static real xnodcf, delmo, x7thm1, x3thm1, x1mth2;
static real aodp, eta, omgdot, xnodot;
static double xnodp, xmdot;

static long Isat=0;	/* 16-bit compilers need 'long' integer for higher space catalogue numbers. */
double perigee, period, apogee;

long Icount = 0;
int MaxNR=0;
extern int Set_LS_zero;	/* From deep.c */


void fatal_error(const char *format, ...)
{
va_list arg_ptr;

    fflush(stdout);

    fprintf(stderr, "\nDundee Satellite Lab fatal run-time error:\n");

    va_start(arg_ptr, format);
    vfprintf(stderr, format, arg_ptr);
    va_end(arg_ptr);

    //fprintf(stderr, "\nNow terminating the program...\n");
    //    fflush(stderr);

    //    exit(5);
    return;
}


/* =======================================================================
   The init_sgdp4() function passes all of the required orbital elements to
   the sgdp4() function together with the pre-calculated constants. There is
   some basic error traps and the detemination of the orbital model is made.
   For near-earth satellites (xnodp < 225 minutes according to the NORAD
   classification) the SGP4 model is used, with truncated terms for low
   perigee heights when the drag terms are high. For deep-space satellites
   the SDP4 model is used and the deep-space terms initialised (a slow
   process). For orbits with an eccentricity of less than ECC_EPS the model
   reverts to a very basic circular model. This is not physically meaningfull
   but such a circluar orbit is not either! It is fast though.
   Callinr arguments:

   orb      : Input, structure with the orbital elements from NORAD 2-line
              element data in radian form.

   The return value indicates the orbital model used.
   ======================================================================= */

int init_sgdp4(orbit_t *orb)
{
LOCAL_REAL theta2, theta4, xhdot1, x1m5th;
LOCAL_REAL s4, del1, del0;
LOCAL_REAL betao, betao2, coef, coef1;
LOCAL_REAL etasq, eeta, qoms24;
LOCAL_REAL pinvsq, tsi, psisq, c1sq;
LOCAL_DOUBLE a0, a1, epoch;
real temp0, temp1, temp2, temp3;
long iday, iyear;

    /* Copy over elements. */
    /* Convert year to Gregorian with century as 1994 or 94 type ? */

    iyear  = (long)orb->ep_year;

    if (iyear < 1957)
        {
        /* Assume 0 and 100 both refer to 2000AD */
        iyear += (iyear < 57 ? 2000 : 1900);
        }

    if (iyear < 1901 || iyear > 2099)
        {
        fatal_error("init_sgdp4: Satellite ep_year error %ld", iyear);
        imode = SGDP4_ERROR;
        return imode;
        }

	Isat = orb->satno;

    /* Compute days from 1st Jan 1900 (works 1901 to 2099 only). */

    iday = ((iyear - 1901)*1461L)/4L + 364L + 1L;

    SGDP4_jd0 = JD1900 + iday + (orb->ep_day - 1.0);  /* Julian day number. */

    epoch  = (iyear - 1900) * 1.0e3 + orb->ep_day; /* YYDDD.DDDD as from 2-line. */

#ifdef DEBUG
    fprintf(stderr, "Epoch = %f SGDP4_jd0 = %f\n", epoch, SGDP4_jd0);
#endif

    eo     = (real)orb->ecc;
    xno    = (double)orb->rev * TWOPI/XMNPDA;   /* Radian / unit time. */
    xincl  = (real)orb->eqinc;
    xnodeo = (real)orb->ascn;
    omegao = (real)orb->argp;
    xmo    = (real)orb->mnan;
    bstar  = (real)orb->bstar;

    /* A few simple error checks here. */

    if (eo < (real)0.0 || eo > ECC_LIMIT_HIGH)
        {
        fatal_error("init_sgdp4: Eccentricity out of range for %ld (%le)", Isat, (double)eo);
        imode = SGDP4_ERROR;
        return imode;
        }

    if (xno < 0.035*TWOPI/XMNPDA || xno > 18.0*TWOPI/XMNPDA)
        {
        fatal_error("init_sgdp4: Mean motion out of range %ld (%le)", Isat, xno);
        imode = SGDP4_ERROR;
        return imode;
        }

    if (xincl < (real)0.0 || xincl > (real)PI)
        {
        fatal_error("init_sgdp4: Equatorial inclination out of range %ld (%le)", Isat, DEG(xincl));
        imode = SGDP4_ERROR;
        return imode;
        }

    /* Start the initialisation. */

    if (eo < ECC_ZERO)
        imode = SGDP4_ZERO_ECC; /* Special mode for "ideal" circular orbit. */
    else
        imode = SGDP4_NOT_INIT;

    /*
    Recover original mean motion (xnodp) and semimajor axis (aodp)
    from input elements.
    */

    SINCOS(xincl, &sinIO, &cosIO);

    theta2 = cosIO * cosIO;
    theta4 = theta2 * theta2;
    x3thm1 = (real)3.0 * theta2 - (real)1.0;
    x1mth2 = (real)1.0 - theta2;
    x7thm1 = (real)7.0 * theta2 - (real)1.0;

    a1 = pow(XKE / xno, TOTHRD);
    betao2 = (real)1.0 - eo * eo;
    betao = SQRT(betao2);
    temp0 = (real)(1.5 * CK2) * x3thm1 / (betao * betao2);
    del1 = temp0 / (a1 * a1);
    a0 = a1 * (1.0 - del1 * (1.0/3.0 + del1 * (1.0 + del1 * 134.0/81.0)));
    del0 = temp0 / (a0 * a0);
    xnodp = xno / (1.0 + del0);
    aodp = (real)(a0 / (1.0 - del0));
    perigee = (aodp * (1.0 - eo) - AE) * XKMPER;
    apogee = (aodp * (1.0 + eo) - AE) * XKMPER;
    period = (TWOPI * 1440.0 / XMNPDA) / xnodp;

    /*
    printf("Perigee = %lf km period = %lf min del0 = %e\n",
              perigee, period, del0);
    */
    if (perigee <= 0.0)
        {
		fprintf(stderr, "# Satellite %ld sub-orbital (apogee = %.1f km, perigee = %.1f km)\n", Isat, apogee, perigee);
        }

    if (imode == SGDP4_ZERO_ECC) return imode;

    if (period >= 225.0 && Set_LS_zero < 2)
        {
        imode = SGDP4_DEEP_NORM; /* Deep-Space model(s). */
        }
    else if (perigee < 220.0)
        {
        /*
        For perigee less than 220 km the imode flag is set so the
        equations are truncated to linear variation in sqrt A and
        quadratic variation in mean anomaly. Also the c3 term, the
        delta omega term and the delta m term are dropped.
        */
        imode = SGDP4_NEAR_SIMP;    /* Near-space, simplified equations. */
        }
    else
        {
        imode = SGDP4_NEAR_NORM;    /* Near-space, normal equations. */
        }

    /* For perigee below 156 km the values of S and QOMS2T are altered */

    if (perigee < 156.0)
        {
		s4 = (real)(perigee - 78.0);

        if(s4 < (real)20.0)
        	{
			fprintf(stderr, "# Very low s4 constant for sat %ld (perigee = %.2f)\n", Isat, perigee);
        	s4 = (real)20.0;
			}
		else
			{
			fprintf(stderr, "# Changing s4 constant for sat %ld (perigee = %.2f)\n", Isat, perigee);
			}

        qoms24 = POW4((real)((120.0 - s4) * (AE / XKMPER)));
        s4 = (real)(s4 / XKMPER + AE);
        }
    else
        {
        s4 = KS;
        qoms24 = QOMS2T;
        }

    pinvsq = (real)1.0 / (aodp * aodp * betao2 * betao2);
    tsi = (real)1.0 / (aodp - s4);
    eta = aodp * eo * tsi;
    etasq = eta * eta;
    eeta = eo * eta;
    psisq = FABS((real)1.0 - etasq);
    coef = qoms24 * POW4(tsi);
    coef1 = coef / POW(psisq, 3.5);

    c2 = coef1 * (real)xnodp * (aodp *
         ((real)1.0 + (real)1.5 * etasq + eeta * ((real)4.0 + etasq)) +
         (real)(0.75 * CK2) * tsi / psisq * x3thm1 *
         ((real)8.0 + (real)3.0 * etasq * ((real)8.0 + etasq)));

    c1 = bstar * c2;

    c4 = (real)2.0 * (real)xnodp * coef1 * aodp * betao2 * (eta *
         ((real)2.0 + (real)0.5 * etasq) + eo * ((real)0.5 + (real)2.0 *
         etasq) - (real)(2.0 * CK2) * tsi / (aodp * psisq) * ((real)-3.0 *
         x3thm1 * ((real)1.0 - (real)2.0 * eeta + etasq *
         ((real)1.5 - (real)0.5 * eeta)) + (real)0.75 * x1mth2 * ((real)2.0 *
         etasq - eeta * ((real)1.0 + etasq)) * COS((real)2.0 * omegao)));

	c5 = c3 = omgcof = (real)0.0;

    if (imode == SGDP4_NEAR_NORM)
        {
        /* BSTAR drag terms for normal near-space 'normal' model only. */
        c5 = (real)2.0 * coef1 * aodp * betao2 *
             ((real)1.0 + (real)2.75 * (etasq + eeta) + eeta * etasq);

        if(eo > ECC_ALL)
        	{
			c3 = coef * tsi * a3ovk2 * (real)xnodp * (real)AE * sinIO / eo;
			}

        omgcof = bstar * c3 * COS(omegao);
        }

    temp1 = (real)(3.0 * CK2) * pinvsq * (real)xnodp;
    temp2 = temp1 * CK2 * pinvsq;
    temp3 = (real)(1.25 * CK4) * pinvsq * pinvsq * (real)xnodp;

    xmdot = xnodp + ((real)0.5 * temp1 * betao * x3thm1 + (real)0.0625 *
            temp2 * betao * ((real)13.0 - (real)78.0 * theta2 +
            (real)137.0 * theta4));

    x1m5th = (real)1.0 - (real)5.0 * theta2;

    omgdot = (real)-0.5 * temp1 * x1m5th + (real)0.0625 * temp2 *
             ((real)7.0 - (real)114.0 * theta2 + (real)395.0 * theta4) +
             temp3 * ((real)3.0 - (real)36.0 * theta2 + (real)49.0 * theta4);

    xhdot1 = -temp1 * cosIO;
    xnodot = xhdot1 + ((real)0.5 * temp2 * ((real)4.0 - (real)19.0 * theta2) +
             (real)2.0 * temp3 * ((real)3.0 - (real)7.0 * theta2)) * cosIO;

	xmcof = (real)0.0;
    if(eo > ECC_ALL)
    	{
    	xmcof = (real)(-TOTHRD * AE) * coef * bstar / eeta;
		}

    xnodcf = (real)3.5 * betao2 * xhdot1 * c1;
    t2cof = (real)1.5 * c1;

	/* Check for possible divide-by-zero for X/(1+cosIO) when calculating xlcof */
	temp0 = (real)1.0 + cosIO;

	if(fabs(temp0) < EPS_COSIO) temp0 = (real)SIGN(EPS_COSIO, temp0);

    xlcof = (real)0.125 * a3ovk2 * sinIO *
            ((real)3.0 + (real)5.0 * cosIO) / temp0;

    aycof = (real)0.25 * a3ovk2 * sinIO;

    SINCOS(xmo, &sinXMO, &cosXMO);
    delmo = CUBE((real)1.0 + eta * cosXMO);

    if (imode == SGDP4_NEAR_NORM)
        {
        c1sq = c1 * c1;
        d2 = (real)4.0 * aodp * tsi * c1sq;
        temp0 = d2 * tsi * c1 / (real)3.0;
        d3 = ((real)17.0 * aodp + s4) * temp0;
        d4 = (real)0.5 * temp0 * aodp * tsi * ((real)221.0 * aodp +
             (real)31.0 * s4) * c1;
        t3cof = d2 + (real)2.0 * c1sq;
        t4cof = (real)0.25 * ((real)3.0 * d3 + c1 * ((real)12.0 * d2 +
                (real)10.0 * c1sq));
        t5cof = (real)0.2 * ((real)3.0 * d4 + (real)12.0 * c1 * d3 +
                (real)6.0 * d2 * d2 + (real)15.0 * c1sq * ((real)2.0 *
                d2 + c1sq));
        }
    else if (imode == SGDP4_DEEP_NORM)
        {
#ifdef NO_DEEP_SPACE
        fatal_error("init_sgdp4: Deep space equations not supported");
#else
        imode = SGDP4_dpinit(epoch, omegao, xnodeo, xmo, eo, xincl,
                              aodp, xmdot, omgdot, xnodot, xnodp);
#endif /* !NO_DEEP_SPACE */
        }

return imode;
}

/* =======================================================================
   The sgdp4() function computes the Keplarian elements that describe the
   position and velocity of the satellite. Depending on the initialisation
   (and the compile options) the deep-space perturbations are also included
   allowing sensible predictions for most satellites. These output elements
   can be transformed to Earth Centered Inertial coordinates (X-Y-Z) and/or
   to sub-satellite latitude and longitude as required. The terms for the
   velocity solution are often not required so the 'withvel' flag can be used
   to by-pass that step as required. This function is normally called through
   another since the input 'tsince' is the time from epoch.
   Calling arguments:

   tsince   : Input, time from epoch (minutes).

   withvel  : Input, non-zero if velocity terms required.

   kep      : Output, the Keplarian position / velocity of the satellite.

   The return value indicates the orbital mode used.

   ======================================================================= */

int sgdp4(double tsince, int withvel, kep_t *kep)
{
LOCAL_REAL rk, uk, xnodek, xinck, em, xinc;
LOCAL_REAL xnode, delm, axn, ayn, omega;
LOCAL_REAL capu, epw, elsq, invR, beta2, betal;
LOCAL_REAL sinu, sin2u, cosu, cos2u;
LOCAL_REAL a, e, r, u, pl;
LOCAL_REAL sinEPW, cosEPW, sinOMG, cosOMG;
LOCAL_DOUBLE xmp, xl, xlt;
const int MAXI = 10;

#ifndef NO_DEEP_SPACE
LOCAL_DOUBLE xn, xmam;
#endif  /* !NO_DEEP_SPACE */

real esinE, ecosE, maxnr;
real temp0, temp1, temp2, temp3;
real tempa, tempe, templ;
int ii;

#ifdef SGDP4_SNGL
real ts = (real)tsince;
#else
#define ts tsince
#endif /* ! SGDP4_SNGL */

	/* Update for secular gravity and atmospheric drag. */

	em = eo;
	xinc = xincl;

	xmp   = (double)xmo + xmdot * tsince;
	xnode = xnodeo + ts * (xnodot + ts * xnodcf);
	omega = omegao + omgdot * ts;

	switch(imode)
		{
		case SGDP4_ZERO_ECC:
			/* Not a "real" orbit but OK for fast computation searches. */
			kep->smjaxs = kep->radius = (double)aodp * XKMPER/AE;
			kep->theta = fmod(PI + xnodp * tsince, TWOPI) - PI;
			kep->eqinc = (double)xincl;
			kep->ascn = xnodeo;

			kep->argp = 0;
			kep->ecc = 0;

			kep->rfdotk = 0;
			if(withvel)
				 kep->rfdotk = aodp * xnodp * (XKMPER/AE*XMNPDA/86400.0); /* For km/sec */
			else
				 kep->rfdotk = 0;

		return imode;

		case SGDP4_NEAR_SIMP:
			tempa = (real)1.0 - ts * c1;
			tempe = bstar * ts * c4;
			templ = ts * ts * t2cof;
			a = aodp * tempa * tempa;
			e = em - tempe;
			xl = xmp + omega + xnode + xnodp * templ;
		break;

		case SGDP4_NEAR_NORM:
			delm  = xmcof * (CUBE((real)1.0 + eta * COS(xmp)) - delmo);
			temp0 = ts * omgcof + delm;
			xmp   += (double)temp0;
			omega -= temp0;
			tempa = (real)1.0 - (ts * (c1 + ts * (d2 + ts * (d3 + ts * d4))));
			tempe = bstar * (c4 * ts + c5 * (SIN(xmp) - sinXMO));
			templ = ts * ts * (t2cof + ts * (t3cof + ts * (t4cof + ts * t5cof)));
			//xmp   += (double)temp0;
			a = aodp * tempa * tempa;
			e = em - tempe;
			xl = xmp + omega + xnode + xnodp * templ;
		break;

#ifndef NO_DEEP_SPACE
		case SGDP4_DEEP_NORM:
		case SGDP4_DEEP_RESN:
		case SGDP4_DEEP_SYNC:
			tempa = (real)1.0 - ts * c1;
			tempe = bstar * ts * c4;
			templ = ts * ts * t2cof;
			xn = xnodp;

			SGDP4_dpsec(&xmp, &omega, &xnode, &em, &xinc, &xn, tsince);

			a = POW(XKE / xn, TOTHRD) * tempa * tempa;
			e = em - tempe;
			xmam = xmp + xnodp * templ;

			SGDP4_dpper(&e, &xinc, &omega, &xnode, &xmam, tsince);

			if (xinc < (real)0.0)
				{
				xinc = (-xinc);
				xnode += (real)PI;
				omega -= (real)PI;
				}

			xl = xmam + omega + xnode;

			/* Re-compute the perturbed values. */
			SINCOS(xinc, &sinIO, &cosIO);

			{
			real theta2 = cosIO * cosIO;

			x3thm1 = (real)3.0 * theta2 - (real)1.0;
			x1mth2 = (real)1.0 - theta2;
			x7thm1 = (real)7.0 * theta2 - (real)1.0;

			/* Check for possible divide-by-zero for X/(1+cosIO) when calculating xlcof */
			temp0 = (real)1.0 + cosIO;

			if(fabs(temp0) < EPS_COSIO) temp0 = (real)SIGN(EPS_COSIO, temp0);

    		xlcof = (real)0.125 * a3ovk2 * sinIO *
            			((real)3.0 + (real)5.0 * cosIO) / temp0;

			aycof = (real)0.25 * a3ovk2 * sinIO;
			}

		break;
#endif /* ! NO_DEEP_SPACE */

		default:
			fatal_error("sgdp4: Orbit not initialised");
			return SGDP4_ERROR;
		}

	if(a < (real)1.0)
		{
		fprintf(stderr, "sgdp4: Satellite %05ld crashed at %.3f (a = %.3f Earth radii)\n", Isat, ts, a);
        return SGDP4_ERROR;
		}

	if(e < ECC_LIMIT_LOW)
		{
		fprintf(stderr, "sgdp4: Satellite %05ld modified eccentricity too low (ts = %.3f, e = %e < %e)\n", Isat, ts, e, ECC_LIMIT_LOW);
        return SGDP4_ERROR;
		}

	if(e < ECC_EPS)
		{
		/*fprintf(stderr, "# ecc %f at %.3f for for %05ld\n", e, ts, Isat);*/
		e = ECC_EPS;
		}
	else if(e > ECC_LIMIT_HIGH)
		{
		/*fprintf(stderr, "# ecc %f at %.3f for for %05ld\n", e, ts, Isat);*/
		e = ECC_LIMIT_HIGH;
		}

	beta2 = (real)1.0 - e * e;

	/* Long period periodics */
	SINCOS(omega, &sinOMG, &cosOMG);

	temp0 = (real)1.0 / (a * beta2);
	axn = e * cosOMG;
	ayn = e * sinOMG + temp0 * aycof;
	xlt = xl + temp0 * xlcof * axn;

	elsq = axn * axn + ayn * ayn;
	if (elsq >= (real)1.0)
		{
		fprintf(stderr, "sgdp4: SQR(e) >= 1 (%.3f at tsince = %.3f for sat %05ld)\n", elsq, tsince, Isat);
		return SGDP4_ERROR;
		}

	/* Sensibility check for N-R correction. */
	kep->ecc = sqrt(elsq);

	/*
	 * Solve Kepler's equation using Newton-Raphson root solving. Here 'capu' is
	 * almost the "Mean anomaly", initialise the "Eccentric Anomaly" term 'epw'.
	 * The fmod() saves reduction of angle to +/-2pi in SINCOS() and prevents
	 * convergence problems.
	 *
	 * Later modified to support 2nd order NR method which saves roughly 1 iteration
	 * for only a couple of arithmetic operations.
	 */

	epw = capu = fmod(xlt - xnode, TWOPI);

	maxnr = kep->ecc;

	for(ii = 0; ii < MAXI; ii++)
		{
		double nr, f, df;
		SINCOS(epw, &sinEPW, &cosEPW);

		ecosE = axn * cosEPW + ayn * sinEPW;
		esinE = axn * sinEPW - ayn * cosEPW;

		f = capu - epw + esinE;
		if (fabs(f) < NR_EPS) break;

		df = 1.0 - ecosE;

		/* 1st order Newton-Raphson correction. */
		nr = f / df;

		if (ii == 0 && FABS(nr) > 1.25*maxnr)
			nr = SIGN(maxnr, nr);
#if 1
		/* 2nd order Newton-Raphson correction. */
		else
			nr = f / (df + 0.5*esinE*nr);	/* f/(df - 0.5*d2f*f/df) */
#endif

		epw += nr;  /* Newton-Raphson correction of -F/DF. */
		//if (fabs(nr) < NR_EPS) break;
		}

	/* Short period preliminary quantities */
	temp0 = (real)1.0 - elsq;
	betal = SQRT(temp0);
	pl = a * temp0;
	r = a * ((real)1.0 - ecosE);
	invR = (real)1.0 / r;
	temp2 = a * invR;
	temp3 = (real)1.0 / ((real)1.0 + betal);
	cosu = temp2 * (cosEPW - axn + ayn * esinE * temp3);
	sinu = temp2 * (sinEPW - ayn - axn * esinE * temp3);
	u = ATAN2(sinu, cosu);
	sin2u = (real)2.0 * sinu * cosu;
	cos2u = (real)2.0 * cosu * cosu - (real)1.0;
	temp0 = (real)1.0 / pl;
	temp1 = CK2 * temp0;
	temp2 = temp1 * temp0;

	/* Update for short term periodics to position terms. */

	rk = r * ((real)1.0 - (real)1.5 * temp2 * betal * x3thm1) + (real)0.5 * temp1 * x1mth2 * cos2u;
	uk = u - (real)0.25 * temp2 * x7thm1 * sin2u;
	xnodek = xnode + (real)1.5 * temp2 * cosIO * sin2u;
	xinck = xinc + (real)1.5 * temp2 * cosIO * sinIO * cos2u;

	if(rk < (real)1.0)
		{
#if 1
		fprintf(stderr, "sgdp4: Satellite %05ld crashed at %.3f (rk = %.3f Earth radii)\n", Isat, ts, rk);
#endif
        return SGDP4_ERROR;
		}

	kep->radius	= rk * XKMPER/AE;    /* Into km */
	kep->theta	= uk;
	kep->eqinc	= xinck;
	kep->ascn	= xnodek;
	kep->argp	= omega;
	kep->smjaxs	= a * XKMPER/AE;

	/* Short period velocity terms ?. */
	if (withvel)
		{
		/* xn = XKE / pow(a, 1.5); */
		temp0 = SQRT(a);
		temp2 = (real)XKE / (a * temp0);

		kep->rdotk = ((real)XKE * temp0 * esinE * invR -
					temp2 * temp1 * x1mth2 * sin2u) *
					(XKMPER/AE*XMNPDA/86400.0); /* Into km/sec */

		kep->rfdotk = ((real)XKE * SQRT(pl) * invR + temp2 * temp1 *
					(x1mth2 * cos2u + (real)1.5 * x3thm1)) *
					(XKMPER/AE*XMNPDA/86400.0);
		}
    else
		{
		kep->rdotk = kep->rfdotk = 0;
		}

#ifndef SGDP4_SNGL
#undef ts
#endif

return imode;
}

/* ====================================================================

   Transformation from "Kepler" type coordinates to cartesian XYZ form.
   Calling arguments:

   K    : Kepler structure as filled by sgdp4();

   pos  : XYZ structure for position.

   vel  : same for velocity.

   ==================================================================== */

void kep2xyz(kep_t *K, xyz_t *pos, xyz_t *vel)
{
real xmx, xmy;
real ux, uy, uz, vx, vy, vz;
real sinT, cosT, sinI, cosI, sinS, cosS;

	/* Orientation vectors for X-Y-Z format. */

	SINCOS((real)K->theta, &sinT, &cosT);
	SINCOS((real)K->eqinc, &sinI, &cosI);
	SINCOS((real)K->ascn,  &sinS, &cosS);

	xmx = -sinS * cosI;
	xmy =  cosS * cosI;

	ux =  xmx * sinT + cosS * cosT;
	uy =  xmy * sinT + sinS * cosT;
	uz = sinI * sinT;

	/* Position and velocity */

	if(pos != NULL)
		{
		pos->x = K->radius * ux;
		pos->y = K->radius * uy;
		pos->z = K->radius * uz;
		}

	if(vel != NULL)
		{
		vx =  xmx * cosT - cosS * sinT;
		vy =  xmy * cosT - sinS * sinT;
		vz = sinI * cosT;

		vel->x = K->rdotk * ux + K->rfdotk * vx;
		vel->y = K->rdotk * uy + K->rfdotk * vy;
		vel->z = K->rdotk * uz + K->rfdotk * vz;
		}

}

/* ======================================================================
   Compute the satellite position and/or velocity for a given time (in the
   form of Julian day number.)
   Calling arguments are:

   jd   : Time as Julian day number.

   pos  : Pointer to posiition vector, km (NULL if not required).

   vel  : Pointer to velocity vector, km/sec (NULL if not required).

   ====================================================================== */

int satpos_xyz(double jd, xyz_t *pos, xyz_t *vel)
{
kep_t K;
int withvel, rv;
double tsince;

	tsince = (jd - SGDP4_jd0) * XMNPDA;

#ifdef DEBUG
	fprintf(stderr, "Tsince = %f\n", tsince);
#endif

	if(vel != NULL)
		withvel = 1;
	else
		withvel = 0;

	rv = sgdp4(tsince, withvel, &K);

	kep2xyz(&K, pos, vel);

return rv;
}

/* ==================== End of file sgdp4.c ========================== */
