/* > deep.c
 *
 * 1.00 around 1980 - Felix R. Hoots & Ronald L. Roehrich, from original
 *                    DEEP.FOR used in the SGP deep-space models SDP4
 *                    and SDP8.
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
 * 2.00 psc Mon Dec 19 1994 - Translated from FORTRAN into 'C' (of sorts).
 *
 * 2.01 psc Wed Dec 21 1994 - Re-write of the secular integrator from a
 *                            messy FORTRAN block in to something which
 *                            (hopefully!) is understandable.
 *
 * 2.02 psc Thu Dec 22 1994 - Final mods and tested against the FORTRAN
 *                            version (using ~12 hour resonant and
 *                            geostationary (~24 hour) elements).
 *
 * 2.03 psc Mon Jan 02 1995 - Some additional refinements and error-traps.
 *
 * 3.00 psc Mon May 29 1995 - Cleaned up for general use & distrabution (to
 *                            remove Dundee specific features).
 *
 * 3.01 psc Mon Jan 12 2004 - Final fix agreed for "Lyddane bug".
 * 3.02 psc Mon Jul 03 2006 - Extended range "Lyddane bug" fix.
 * 3.03 psc Tue Jul 04 2006 - Bug fix for extended range "Lyddane bug" fix.
 */
#include <stdio.h>
#ifndef _MATH_H
#define _MATH_H
#include <math.h>
#endif
#ifndef _STDLIB_H
#define _STDLIB_H 
#include <stdlib.h>
#endif

static const char SCCSid[] = "@(#)deep.c        3.03 (C) 1995 psc SatLib: Deep Space effects";

#ifndef NO_DEEP_SPACE

#include "sgdp4h.h"

extern long Isat;
int Set_LS_zero = 0;	/* Set to 1 to zero Lunar-Solar terms at epoch. */

/* ======================= Function prototypes ====================== */

static void dot_terms_calculated(void);
static void compute_LunarSolar(double tsince);
static void thetag(double ep, real *thegr, double *days50);

/* ===================== Strange constants, etc ===================== */

#define ZNS         ((real)1.19459e-5)
#define C1SS        ((real)2.9864797e-6)
#define ZES         ((real)0.01675)

#define ZNL         ((real)1.5835218e-4)
#define C1L         ((real)4.7968065e-7)
#define ZEL         ((real)0.0549)

#define ZCOSIS      ((real)0.91744867)
#define ZSINIS      ((real)0.39785416)
#define ZCOSGS      ((real)0.1945905)
#define ZSINGS      ((real)-0.98088458)

#define Q22         ((real)1.7891679e-6)
#define Q31         ((real)2.1460748e-6)
#define Q33         ((real)2.2123015e-7)

#define G22         ((real)5.7686396)
#define G32         ((real)0.95240898)
#define G44         ((real)1.8014998)
#define G52         ((real)1.0508330)
#define G54         ((real)4.4108898)

#define ROOT22      ((real)1.7891679e-6)
#define ROOT32      ((real)3.7393792e-7)
#define ROOT44      ((real)7.3636953e-9)
#define ROOT52      ((real)1.1428639e-7)
#define ROOT54      ((real)2.1765803e-9)

#define THDT        ((real)4.37526908801129966e-3)
//#define THDT		((real)0.0043752691)

#define STEP            720.0
#define MAX_INTEGRATE	(STEP * 10000)
#define SIN_EPS			(real)(1.0e-12)

/* ======= Global variables used by dpsec(), from dpinit(). ======== */

static real eo;         /* copy of original eccentricity. */
static real xincl;      /* copy of original equatorial inclination. */

static int isynfl=0, iresfl=0;

static double atime, xli, xni, xnq, xfact;

static real ssl, ssg, ssh, sse, ssi;
static real xlamo, omegaq, omgdt, thgr;
static real del1, del2, del3, fasx2, fasx4, fasx6;
static real d2201, d2211, d3210, d3222, d4410, d4422;
static real d5220, d5232, d5421, d5433;

static real xnddt, xndot, xldot;	/* Integrator terms. */
static real xnddt0, xndot0, xldot0; /* Integrator at epoch. */

/* ======== Global Variables used by dpper(), from dpinit(). ======= */

static int ilsd=0, ilsz=0;

static real zmos, se2, se3, si2, si3, sl2, sl3, sl4;
static real sgh2, sgh3, sgh4, sh2, sh3;
static real zmol, ee2, e3 ,xi2, xi3, xl2, xl3, xl4;
static real xgh2, xgh3, xgh4, xh2, xh3;

static real pe, pinc, pgh, ph, pl;
static real pgh0, ph0, pe0, pinc0, pl0; /* Added terms to save the epoch values of perturbations. */


/* ==================================================================

   ----------------- DEEP SPACE INITIALIZATION ----------------------

  epoch     : Input, epoch time as YYDDD.DDDD as read from 2-line elements.
  omegao    : Input, argument of perigee from elements, radian.
  xnodeo    : Input, right asc. for ascn node from elements, radian.
  xmo       : Input, mean anomaly from elements, radian.
  orb_eo    : Input, eccentricity from elements, dimentionless.
  orb_xincl : Input, equatorial inclination from elements, radian.
  aodp      : Input, original semi-major axis, earth radii.
  xlldot    : Input, 1st derivative of "mean anomaly" (xmdot), radian/min.
  omgdot    : Input, 1st derivative of arg. per., radian/min.
  xnodot    : Input, 1st derivative of right asc., radian/min.
  xnodp     : Input, original mean motion, radian/min.

   ================================================================== */

int SGDP4_dpinit(double epoch, real omegao, real xnodeo, real xmo,
                 real orb_eo, real orb_xincl, real aodp, double xlldot,
                 real omgdot, real xnodot, double xnodp)
{
LOCAL_DOUBLE ds50, day, xnodce, bfact=0, gam, c;
LOCAL_REAL ctem, sinq, cosq, aqnv, xmao, stem, eqsq, xnoi, ainv2;
LOCAL_REAL zcosg, zsing, zcosi, zsini, zcosh, zsinh;
LOCAL_REAL cosomo, zcosgl, zcoshl, zcosil, sinomo;
LOCAL_REAL xpidot, zsinil, siniq2, cosiq2;
LOCAL_REAL rteqsq, zsinhl, zsingl;
LOCAL_REAL eoc, sgh, g200, bsq, zmo, xno2;
LOCAL_REAL a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
LOCAL_REAL x1, x2, x3, x4, x5, x6, x7, x8;
LOCAL_REAL z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33;
LOCAL_REAL s1, s2, s3, s4, s5, s6, s7, cc, ao, eq, se, shdq, si, sl;
LOCAL_REAL zx, zy, ze, zn;
LOCAL_REAL g201, g211, g310, g300, g322, g410, g422, g520, g533, g521, g532;
LOCAL_REAL f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543;
real siniq, cosiq;
real temp0, temp1;
int ls, imode=0;
int ishq;

    /*
    Copy the supplied orbital elements to "local" (static to this file)
    variables and compute common trig values.
    */
    eq = eo = orb_eo;
    xincl = orb_xincl;

    /* Decide on direct or Lyddane Lunar-Solar perturbations. */
    ilsd = 0;
    if(xincl >= (real)0.2) ilsd = 1;

	/* Drop some terms below 3 deg inclination. */
	ishq = 0;
#define SHQT 0.052359877
	if (xincl >= (real)SHQT) ishq = 1; /* As per reoprt #3. */

    SINCOS(omegao, &sinomo, &cosomo);
    SINCOS(xnodeo, &sinq, &cosq);
    SINCOS(xincl, &siniq, &cosiq);

    if (fabs(siniq) <= SIN_EPS)
        {
        siniq = SIGN(SIN_EPS, siniq);
        }

    cosiq2 = cosiq * cosiq;
    siniq2 = siniq * siniq;

    ao = aodp;
    omgdt = omgdot;
    eqsq = eo * eo;
    bsq = (real)1.0 - eqsq;
    rteqsq = SQRT(bsq);
    thetag(epoch, &thgr, &ds50);

    /*printf("# epoch = %.8f ds50 = %.8f thgr = %f\n", epoch, ds50, DEG(thgr));*/

    xnq = xnodp;
    aqnv = (real)1.0 / ao;
    xmao = xmo;
    xpidot = omgdt + xnodot;
    omegaq = omegao;

    /* INITIALIZE LUNAR SOLAR TERMS */

    day = ds50 + 18261.5;

    xnodce = 4.523602 - day * 9.2422029e-4;
    temp0 = (real)fmod(xnodce, TWOPI);
    SINCOS(temp0, &stem, &ctem);

    zcosil = (real)0.91375164 - ctem * (real)0.03568096;
    zsinil = SQRT((real)1.0 - zcosil * zcosil);
    zsinhl = stem * (real)0.089683511 / zsinil;
    zcoshl = SQRT((real)1.0 - zsinhl * zsinhl);
    c = day * 0.2299715 + 4.7199672;
    gam = day * 0.001944368 + 5.8351514;
    zmol = (real)MOD2PI(c - gam);
    zx = stem * (real)0.39785416 / zsinil;
    zy = zcoshl * ctem + zsinhl * (real)0.91744867 * stem;
    zx = ATAN2(zx, zy);
    zx = (real)fmod(gam + zx - xnodce, TWOPI);
    SINCOS(zx, &zsingl, &zcosgl);
    zmos = (real)MOD2PI(day * 0.017201977 + 6.2565837);

    /* DO SOLAR TERMS */

    zcosg = ZCOSGS;
    zsing = ZSINGS;
    zcosi = ZCOSIS;
    zsini = ZSINIS;
    zcosh = cosq;
    zsinh = sinq;
    cc = C1SS;
    zn = ZNS;
    ze = ZES;
    zmo = zmos;
    xnoi = (real)(1.0 / xnq);

    for(ls = 0; ls < 2; ls++)
        {
        a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        a8 = zsing * zsini;
        a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        a10 = zcosg * zsini;
        a2 = cosiq * a7 + siniq * a8;
        a4 = cosiq * a9 + siniq * a10;
        a5 = -siniq * a7 + cosiq * a8;
        a6 = -siniq * a9 + cosiq * a10;

        x1 = a1 * cosomo + a2 * sinomo;
        x2 = a3 * cosomo + a4 * sinomo;
        x3 = -a1 * sinomo + a2 * cosomo;
        x4 = -a3 * sinomo + a4 * cosomo;
        x5 = a5 * sinomo;
        x6 = a6 * sinomo;
        x7 = a5 * cosomo;
        x8 = a6 * cosomo;

        z31 = x1 * (real)12.0 * x1 - x3 * (real)3.0 * x3;
        z32 = x1 * (real)24.0 * x2 - x3 * (real)6.0 * x4;
        z33 = x2 * (real)12.0 * x2 - x4 * (real)3.0 * x4;
        z1 = (a1 * a1 + a2 * a2) * (real)3.0 + z31 * eqsq;
        z2 = (a1 * a3 + a2 * a4) * (real)6.0 + z32 * eqsq;
        z3 = (a3 * a3 + a4 * a4) * (real)3.0 + z33 * eqsq;
        z11 = a1 * (real)-6.0 * a5 + eqsq * (x1 * (real)-24.0 * x7 - x3 *
                (real)6.0 * x5);
        z12 = (a1 * a6 + a3 * a5) * (real)-6.0 + eqsq * ((x2 * x7 +
              x1 * x8) * (real)-24.0 - (x3 * x6 + x4 * x5) * (real)6.0);
        z13 = a3 * (real)-6.0 * a6 + eqsq * (x2 * (real)-24.0 * x8 - x4 *
                (real)6.0 * x6);
        z21 = a2 * (real)6.0 * a5 + eqsq * (x1 * (real)24.0 * x5 -
              x3 * (real)6.0 * x7);
        z22 = (a4 * a5 + a2 * a6) * (real)6.0 + eqsq * ((x2 * x5 + x1 * x6) *
                (real)24.0 - (x4 * x7 + x3 * x8) * (real)6.0);
        z23 = a4 * (real)6.0 * a6 + eqsq * (x2 * (real)24.0 * x6 - x4 *
                (real)6.0 * x8);
        z1 = z1 + z1 + bsq * z31;
        z2 = z2 + z2 + bsq * z32;
        z3 = z3 + z3 + bsq * z33;
        s3 = cc * xnoi;
        s2 = s3 * (real)-0.5 / rteqsq;
        s4 = s3 * rteqsq;
        s1 = eq * (real)-15.0 * s4;
        s5 = x1 * x3 + x2 * x4;
        s6 = x2 * x3 + x1 * x4;
        s7 = x2 * x4 - x1 * x3;
        se = s1 * zn * s5;
        si = s2 * zn * (z11 + z13);
        sl = -zn * s3 * (z1 + z3 - (real)14.0 - eqsq * (real)6.0);
        sgh = s4 * zn * (z31 + z33 - (real)6.0);

        shdq = 0;
        if(ishq)
        	{
			real sh = -zn * s2 * (z21 + z23);
			shdq = sh / siniq;
			}

        ee2 = s1 * (real)2.0 * s6;
        e3 = s1 * (real)2.0 * s7;
        xi2 = s2 * (real)2.0 * z12;
        xi3 = s2 * (real)2.0 * (z13 - z11);
        xl2 = s3 * (real)-2.0 * z2;
        xl3 = s3 * (real)-2.0 * (z3 - z1);
        xl4 = s3 * (real)-2.0 * ((real)-21.0 - eqsq * (real)9.0) * ze;
        xgh2 = s4 * (real)2.0 * z32;
        xgh3 = s4 * (real)2.0 * (z33 - z31);
        xgh4 = s4 * (real)-18.0 * ze;
        xh2 = s2 * (real)-2.0 * z22;
        xh3 = s2 * (real)-2.0 * (z23 - z21);

        if (ls == 1) break;

        /* DO LUNAR TERMS */

        sse = se;
        ssi = si;
        ssl = sl;
        ssh = shdq;
        ssg = sgh - cosiq * ssh;
        se2 = ee2;
        si2 = xi2;
        sl2 = xl2;
        sgh2 = xgh2;
        sh2 = xh2;
        se3 = e3;
        si3 = xi3;
        sl3 = xl3;
        sgh3 = xgh3;
        sh3 = xh3;
        sl4 = xl4;
        sgh4 = xgh4;
        zcosg = zcosgl;
        zsing = zsingl;
        zcosi = zcosil;
        zsini = zsinil;
        zcosh = zcoshl * cosq + zsinhl * sinq;
        zsinh = sinq * zcoshl - cosq * zsinhl;
        zn = ZNL;
        cc = C1L;
        ze = ZEL;
        zmo = zmol;
        }

    sse += se;
    ssi += si;
    ssl += sl;
    ssg += sgh - cosiq * shdq;
    ssh += shdq;

    if (xnq < 0.0052359877 && xnq > 0.0034906585)
        {
        /* 24h SYNCHRONOUS RESONANCE TERMS INITIALIZATION */
        iresfl = 1;
        isynfl = 1;
        g200 = eqsq * (eqsq * (real)0.8125 - (real)2.5) + (real)1.0;
        g310 = eqsq * (real)2.0 + (real)1.0;
        g300 = eqsq * (eqsq * (real)6.60937 - (real)6.0) + (real)1.0;
        f220 = (cosiq + (real)1.0) * (real)0.75 * (cosiq + (real)1.0);
        f311 = siniq * (real)0.9375 * siniq * (cosiq * (real)3.0 +
                (real)1.0) - (cosiq + (real)1.0) * (real)0.75;
        f330 = cosiq + (real)1.0;
        f330 = f330 * (real)1.875 * f330 * f330;
        del1 = (real)3.0 * (real)(xnq * xnq * aqnv * aqnv);
        del2 = del1 * (real)2.0 * f220 * g200 * Q22;
        del3 = del1 * (real)3.0 * f330 * g300 * Q33 * aqnv;
        del1 = del1 * f311 * g310 * Q31 * aqnv;
        fasx2 = (real)0.13130908;
        fasx4 = (real)2.8843198;
        fasx6 = (real)0.37448087;
        xlamo = xmao + xnodeo + omegao - thgr;
        bfact = xlldot + xpidot - THDT;
        bfact += (double)(ssl + ssg + ssh);
        }
    else if (xnq >= 0.00826 && xnq <= 0.00924 && eq >= (real)0.5)
        {
        /* GEOPOTENTIAL RESONANCE INITIALIZATION FOR 12 HOUR ORBITS */
        iresfl = 1;
        isynfl = 0;
        eoc = eq * eqsq;
        g201 = (real)-0.306 - (eq - (real)0.64) * (real)0.44;

        if (eq <= (real)0.65)
            {
            g211 = (real)3.616 - eq * (real)13.247 + eqsq * (real)16.29;
            g310 = eq * (real)117.39 - (real)19.302 - eqsq * (real)228.419 + eoc * (real)156.591;
            g322 = eq * (real)109.7927 - (real)18.9068 - eqsq * (real)214.6334 + eoc * (real)146.5816;
            g410 = eq * (real)242.694 - (real)41.122 - eqsq * (real)471.094 + eoc * (real)313.953;
            g422 = eq * (real)841.88 - (real)146.407 - eqsq * (real)1629.014 + eoc * (real)1083.435;
            g520 = eq * (real)3017.977 - (real)532.114 - eqsq * 5740.032 + eoc * (real)3708.276;
            }
        else
            {
            g211 = eq * (real)331.819 - (real)72.099 - eqsq * (real)508.738 + eoc * (real)266.724;
            g310 = eq * (real)1582.851 - (real)346.844 - eqsq * (real)2415.925 + eoc * (real)1246.113;
            g322 = eq * (real)1554.908 - (real)342.585 - eqsq * (real)2366.899 + eoc * (real)1215.972;
            g410 = eq * (real)4758.686 - (real)1052.797 - eqsq * (real)7193.992 + eoc * (real)3651.957;
            g422 = eq * (real)16178.11 - (real)3581.69 - eqsq * (real)24462.77 + eoc * (real)12422.52;

            if (eq <= (real)0.715)
                {
                g520 = (real)1464.74 - eq * (real)4664.75 + eqsq * (real)3763.64;
                }
            else
                {
                g520 = eq * (real)29936.92 - (real)5149.66 - eqsq * (real)54087.36 + eoc * (real)31324.56;
                }
            }

        if (eq < (real)0.7)
            {
            g533 = eq * (real)4988.61 - (real)919.2277 - eqsq * (real)9064.77 + eoc * (real)5542.21;
            g521 = eq * (real)4568.6173 - (real)822.71072 - eqsq * (real)8491.4146 + eoc * (real)5337.524;
            g532 = eq * (real)4690.25 - (real)853.666 - eqsq * (real)8624.77 + eoc * (real)5341.4;
            }
        else
            {
            g533 = eq * (real)161616.52 - (real)37995.78 - eqsq * (real)229838.2 + eoc * (real)109377.94;
            g521 = eq * (real)218913.95 - (real)51752.104 - eqsq * (real)309468.16 + eoc * (real)146349.42;
            g532 = eq * (real)170470.89 - (real)40023.88 - eqsq * (real)242699.48 + eoc * (real)115605.82;
            }

        f220 = (cosiq * (real)2.0 + (real)1.0 + cosiq2) * (real)0.75;
        f221 = siniq2 * (real)1.5;
        f321 = siniq * (real)1.875 * ((real)1.0 - cosiq * (real)2.0 - cosiq2 * (real)3.0);
        f322 = siniq * (real)-1.875 * (cosiq * (real)2.0 + (real)1.0 - cosiq2 * (real)3.0);
        f441 = siniq2 * (real)35.0 * f220;
        f442 = siniq2 * (real)39.375 * siniq2;
        f522 = siniq * (real)9.84375 * (siniq2 * ((real)1.0 - cosiq *
                (real)2.0 - cosiq2 * (real)5.0) + (cosiq * (real)4.0 -
                (real)2.0 + cosiq2 * (real)6.0) * (real)0.33333333);
        f523 = siniq * (siniq2 * (real)4.92187512 * ((real)-2.0 - cosiq *
                (real)4.0 + cosiq2 * (real)10.0) + (cosiq * (real)2.0 +
                (real)1.0 - cosiq2 * (real)3.0) * (real)6.56250012);
        f542 = siniq * (real)29.53125 * ((real)2.0 - cosiq * (real)8.0 +
                cosiq2 * (cosiq * (real)8.0 - (real)12.0 + cosiq2 *
                (real)10.0));
        f543 = siniq * (real)29.53125 * ((real)-2.0 - cosiq * (real)8.0 +
                cosiq2 * (cosiq * (real)8.0 + (real)12.0 - cosiq2 *
                (real)10.0));
        xno2 = (real)(xnq * xnq);
        ainv2 = aqnv * aqnv;
        temp1 = xno2 * (real)3.0 * ainv2;
        temp0 = temp1 * ROOT22;
        d2201 = temp0 * f220 * g201;
        d2211 = temp0 * f221 * g211;
        temp1 *= aqnv;
        temp0 = temp1 * ROOT32;
        d3210 = temp0 * f321 * g310;
        d3222 = temp0 * f322 * g322;
        temp1 *= aqnv;
        temp0 = temp1 * (real)2.0 * ROOT44;
        d4410 = temp0 * f441 * g410;
        d4422 = temp0 * f442 * g422;
        temp1 *= aqnv;
        temp0 = temp1 * ROOT52;
        d5220 = temp0 * f522 * g520;
        d5232 = temp0 * f523 * g532;
        temp0 = temp1 * (real)2.0 * ROOT54;
        d5421 = temp0 * f542 * g521;
        d5433 = temp0 * f543 * g533;
        xlamo = xmao + xnodeo + xnodeo - thgr - thgr;
        bfact = xlldot + xnodot + xnodot - THDT - THDT;
        bfact += (double)(ssl + ssh + ssh);
        }
    else
        {
        /* NON RESONANT ORBITS */
        iresfl = 0;
        isynfl = 0;
        }

	if(iresfl == 0)
		{
		/* Non-resonant orbits. */
        imode = SGDP4_DEEP_NORM;
		}
	else
		{
		/* INITIALIZE INTEGRATOR */
		xfact = bfact - xnq;
		xli = (double)xlamo;
		xni = xnq;
		atime = 0.0;

		dot_terms_calculated();

		/* Save the "dot" terms for integrator re-start. */
		xnddt0 = xnddt;
		xndot0 = xndot;
		xldot0 = xldot;

		if (isynfl)
			imode = SGDP4_DEEP_SYNC;
		else
			imode = SGDP4_DEEP_RESN;
		}

	/* Set up for original mode (LS terms at epoch non-zero). */
	ilsz = 0;
	pgh0 = ph0 = pe0 = pinc0 = pl0 = (real)0.0;

	if(Set_LS_zero)
		{
		/* Save the epoch case Lunar-Solar terms to remove this bias for
		 * actual computations later on.
		 * Not sure if this is a good idea.
		 */
		compute_LunarSolar(0.0);

		pgh0	= pgh;
		ph0		= ph;
		pe0		= pe;
		pinc0	= pinc;
		pl0		= pl;
		ilsz	= 1;
		}


return imode;
} /* SGDP4_dpinit */

/* =====================================================================

   ------------- ENTRANCE FOR DEEP SPACE SECULAR EFFECTS ---------------

   xll      : Input/Output, modified "mean anomaly" or "mean longitude".
   omgasm   : Input/Output, modified argument of perigee.
   xnodes   : Input/Output, modified right asc of ascn node.
   em       : Input/Output, modified eccentricity.
   xinc     : Input/Output, modified inclination.

   xn       : Output, modified period from 'xnodp'.

   tsince   : Input, time from epoch (minutes).

   ===================================================================== */

int SGDP4_dpsec(double *xll, real *omgasm, real *xnodes, real *em,
                real *xinc, double *xn, double tsince)
{
LOCAL_DOUBLE delt, ft, xl;
real temp0;

    *xll 	+= ssl * tsince;
    *omgasm += ssg * tsince;
    *xnodes += ssh * tsince;
    *em 	+= sse * tsince;
    *xinc 	+= ssi * tsince;

    if (iresfl == 0) return 0;

    /*
	 * A minor increase in some efficiency can be had by restarting if
	 * the new time is closer to epoch than to the old integrated
	 * time. This also forces a re-start on a change in sign (i.e. going
	 * through zero time) as then we have |tsince - atime| > |tsince|
	 * as well. Second test is for stepping back towards zero, forcing a restart
	 * if close enough rather than integrating to zero.
	 */
#define AHYST 1.0
	/* Most accurate (OK, most _consistant_) method. Restart if need to
	 * integrate 'backwards' significantly from current point.
	 */
	if(fabs(tsince) < STEP ||
	   (atime > 0.0 && tsince < atime - AHYST) ||
	   (atime < 0.0 && tsince > atime + AHYST))
       {
       /* Epoch restart if we are at, or have crossed, tsince==0 */
       atime = 0.0;
       xni = xnq;
       xli = (double)xlamo;
       /* Restore the old "dot" terms. */
       xnddt = xnddt0;
       xndot = xndot0;
       xldot = xldot0;
       }

    ft = tsince - atime;

	if (fabs(ft) > MAX_INTEGRATE)
		{
		printf("SGDP4_dpsec: Integration limit reached\n");
        return -1;
		}

    if (fabs(ft) >= STEP)
        {
        /*
        Do integration if required. Find the step direction to
        make 'atime' catch up with 'tsince'.
        */
        delt = (tsince >= atime ? STEP : -STEP);

        do {
            /* INTEGRATOR (using the last "dot" terms). */
            xli += delt * (xldot + delt * (real)0.5 * xndot);
            xni += delt * (xndot + delt * (real)0.5 * xnddt);
            atime += delt;

            dot_terms_calculated();

            /* Are we close enough now ? */
            ft = tsince - atime;
            } while (fabs(ft) >= STEP);
        }

    xl  = xli + ft * (xldot + ft * (real)0.5 * xndot);
    *xn = xni + ft * (xndot + ft * (real)0.5 * xnddt);

    temp0 = -(*xnodes) + thgr + tsince * THDT;

    if (isynfl == 0)
        *xll = xl + temp0 + temp0;
    else
        *xll = xl - *omgasm + temp0;

return 0;
}  /* SGDP4_dpsec */

/* =====================================================================

   Here we do the "dot" terms for the integrator. Separate function so we
   can call when initialising and save the atime==0.0 values for later
   epoch re-start of the integrator.

   ===================================================================== */

static void dot_terms_calculated(void)
{
LOCAL_DOUBLE x2li, x2omi, xomi;

    /* DOT TERMS CALCULATED */
    if (isynfl)
        {
        xndot = del1 * SIN(xli - fasx2)
              + del2 * SIN((xli - fasx4) * (real)2.0)
              + del3 * SIN((xli - fasx6) * (real)3.0);

        xnddt = del1 * COS(xli - fasx2)
              + del2 * COS((xli - fasx4) * (real)2.0) * (real)2.0
              + del3 * COS((xli - fasx6) * (real)3.0) * (real)3.0;
        }
    else
        {
        xomi = omegaq + omgdt * atime;
        x2omi = xomi + xomi;
        x2li = xli + xli;

        xndot = d2201 * SIN(x2omi + xli - G22)
              + d2211 * SIN(xli - G22)
              + d3210 * SIN(xomi + xli - G32)
              + d3222 * SIN(-xomi + xli - G32)
              + d5220 * SIN(xomi + xli - G52)
              + d5232 * SIN(-xomi + xli - G52)
              + d4410 * SIN(x2omi + x2li - G44)
              + d4422 * SIN(x2li - G44)
              + d5421 * SIN(xomi + x2li - G54)
              + d5433 * SIN(-xomi + x2li - G54);

        xnddt = d2201 * COS(x2omi + xli - G22)
              + d2211 * COS(xli - G22)
              + d3210 * COS(xomi + xli - G32)
              + d3222 * COS(-xomi + xli - G32)
              + d5220 * COS(xomi + xli - G52)
              + d5232 * COS(-xomi + xli - G52)
              + (d4410 * COS(x2omi + x2li - G44)
              +  d4422 * COS(x2li - G44)
              +  d5421 * COS(xomi + x2li - G54)
              +  d5433 * COS(-xomi + x2li - G54)) * (real)2.0;
        }

    xldot = (real)(xni + xfact);
    xnddt *= xldot;

} /* dot_terms_calculated */

/* =====================================================================

   ---------------- ENTRANCES FOR LUNAR-SOLAR PERIODICS ----------------

   em       : Input/Output, modified eccentricity.
   xinc     : Input/Output, modified inclination.
   omgasm   : Input/Output, modified argument of perigee.
   xnodes   : Input/Output, modified right asc of ascn node.
   xll      : Input/Output, modified "mean anomaly" or "mean longitude".
   tsince   : Input, time from epoch (minutes).

   ===================================================================== */

int SGDP4_dpper(real *em, real *xinc, real *omgasm, real *xnodes,
                double *xll, double tsince)
{
real sinis, cosis;

	compute_LunarSolar(tsince);

    *xinc += pinc;
    *em += pe;

    /* Spacetrack report #3 has sin/cos from before perturbations
     * added to xinc (oldxinc), but apparently report # 6 has then
     * from after they are added.
     */
	SINCOS(*xinc, &sinis, &cosis);

    if (ilsd)
		{
		/* APPLY PERIODICS DIRECTLY */
		real tmp_ph;
		tmp_ph = ph / sinis;

		*omgasm += pgh - cosis * tmp_ph;
		*xnodes += tmp_ph;
		*xll	+= pl;
		}
    else
		{
		/* APPLY PERIODICS WITH LYDDANE MODIFICATION */
		LOCAL_REAL alfdp, betdp, dalf, dbet, xls, dls;
		LOCAL_REAL sinok, cosok;
		int ishift;
		real tmp, oldxnode = (*xnodes);

		SINCOS(*xnodes, &sinok, &cosok);
		alfdp = sinis * sinok;
		betdp = sinis * cosok;
		dalf = ph * cosok + pinc * cosis * sinok;
		dbet = -ph * sinok + pinc * cosis * cosok;
		alfdp += dalf;
		betdp += dbet;
		xls = (real)*xll + *omgasm + cosis * *xnodes;
		dls = pl + pgh - pinc * *xnodes * sinis;
		xls += dls;
		*xnodes = ATAN2(alfdp, betdp);

		/* Get perturbed xnodes in to same quadrant as original. */
		ishift = NINT((oldxnode - (*xnodes))/TWOPI);
		*xnodes += (real)(TWOPI * ishift);

		*xll += (double)pl;
		*omgasm = xls - (real)*xll - cosis * (*xnodes);
		}

return 0;
} /* SGDP4_dpper */

/* =====================================================================
   Do the Lunar-Solar terms for the SGDP4_dpper() function (normally only
   every 1/2 hour needed. Seperate function so initialisng could save the
   epoch terms to zero them. Not sure if this is a good thing (some believe
   it the way the equations were intended) as the two-line elements may
   be computed to give the right answer with out this (which I would hope
   as it would make predictions consistant with the 'official' model
   code).
   ===================================================================== */

static void compute_LunarSolar(double tsince)
{
LOCAL_REAL sinzf, coszf;
LOCAL_REAL f2, f3, zf, zm;
LOCAL_REAL sel, sil, ses, sll, sis, sls;
LOCAL_REAL sghs, shs, sghl, shl;

	/* Update Solar terms. */
	zm = zmos + ZNS * tsince;
	zf = zm + ZES * (real)2.0 * SIN(zm);
	SINCOS(zf, &sinzf, &coszf);
	f2 = sinzf * (real)0.5 * sinzf - (real)0.25;
	f3 = sinzf * (real)-0.5 * coszf;
	ses  = se2 * f2 + se3 * f3;
	sis  = si2 * f2 + si3 * f3;
	sls  = sl2 * f2 + sl3 * f3 + sl4 * sinzf;

	sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
	shs  = sh2  * f2 + sh3  * f3;

	/* Update Lunar terms. */
	zm = zmol + ZNL * tsince;
	zf = zm + ZEL * (real)2.0 * SIN(zm);
	SINCOS(zf, &sinzf, &coszf);
	f2 = sinzf * (real)0.5 * sinzf - (real)0.25;
	f3 = sinzf * (real)-0.5 * coszf;
	sel = ee2 * f2 + e3 * f3;
	sil = xi2 * f2 + xi3 * f3;
	sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;

	sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
	shl  = xh2  * f2 + xh3  * f3;

	/* Save computed values to calling structure. */
	pgh  = sghs + sghl;
	ph   = shs + shl;
	pe   = ses + sel;
	pinc = sis + sil;
	pl   = sls + sll;

	if (ilsz)
		{
		/* Correct for previously saved epoch terms. */
		pgh  -= pgh0;
		ph   -= ph0;
		pe   -= pe0;
		pinc -= pinc0;
		pl   -= pl0;
		}

}

/* =====================================================================
   This function converts the epoch time (in the form of YYDDD.DDDDDDDD,
   exactly as it appears in the two-line elements) into days from 00:00:00
   hours Jan 1st 1950 UTC. Also it computes the right ascencion of Greenwich
   at the epoch time, but not in a very accurate manner. However, the same
   method is used here to allow exact comparason with the original FORTRAN
   versions of the programs. The calling arguments are:

   ep       : Input, epoch time of elements (as read from 2-line data).
   thegr    : Output, right ascensionm of Greenwich at epoch, radian.
   days50   : Output, days from Jan 1st 1950 00:00:00 UTC.

   ===================================================================== */

#define THETAG 2

/* Version like sat_code. */
#define J1900	(2451545.5 - 36525. - 1.)
#define SECDAY	(86400.0)

#define C1		(1.72027916940703639E-2)
#define C1P2P	(C1 + TWOPI)
#define THGR70	(1.7321343856509374)
#define FK5R	(5.07551419432269442E-15)


static void thetag(double ep, real *thegr, double *days50)
{
double d;
long n, jy;
double jd, theta;

    jy = (long)((ep + 2.0e-7) * 0.001); /* Extract the year. */
    d = ep - jy * 1.0e3;                /* And then the day of year. */

    /* Assume " 8" is 1980, or more sensibly 2008 ? */
    /*
    if (jy < 10) jy += 80;
    */
    if (jy < 50) jy += 100;

    if (jy < 70)                        /* Fix for leap years ? */
        n = (jy - 72) / 4;
    else
        n = (jy - 69) / 4;

    *days50 = (jy - 70) * 365.0 + 7305.0 + n + d;

	jd = d + J1900 + jy * 365. + ((jy - 1) / 4);

#if THETAG == 0
	/* Original report #3 code. */
    theta = *days50 * 6.3003880987 + 1.72944494;
#elif THETAG == 1
   	{
	/* Method from project pluto code. */
	/* Reference:  The 1992 Astronomical Almanac, page B6. */
 	const double omega_E = 1.00273790934; /* Earth rotations per sidereal day (non-constant) */
  	const double UT = fmod(jd + 0.5, 1.0);
  	double t_cen, GMST;

  	t_cen = (jd - UT - 2451545.0) / 36525.0;
  	GMST = 24110.54841 + t_cen * (8640184.812866 + t_cen * (0.093104 - t_cen * 6.2E-6));

  	GMST = fmod( GMST + SECDAY * omega_E * UT, SECDAY);

  	if(GMST < 0.0) GMST += SECDAY;

  	theta = TWOPI * GMST / SECDAY;
	}
#elif THETAG == 2
	{
	/* Method from SGP4SUB.F code. */
	double ts70, ds70, trfac;
	long ids70;

    ts70 = (*days50) - 7305.0;
    ids70 = (long)(ts70 + 1.0e-8);
    ds70 = ids70;

    trfac = ts70 - ds70;

	/* CALCULATE GREENWICH LOCATION AT EPOCH */
	theta = THGR70 + C1*ds70 + C1P2P*trfac + ts70*ts70*FK5R;
	}
#else
#error 'Unknown method for theta-G calculation'
#endif

	theta = fmod(theta, TWOPI);
    if (theta < 0.0) theta += TWOPI;

    *thegr = (real)theta;

} /* thetag */


#endif /* !NO_DEEP_SPACE */
