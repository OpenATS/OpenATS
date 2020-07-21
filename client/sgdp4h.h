/* > sgdp4h.h
 *
 *
 *     Paul S. Crawford and Andrew R. Brooks
 *     Dundee University
 *
 *                          NOTE !
 *  This code is supplied "as is" and without warranty of any sort.
 *
 * (c) 1994-2004, Paul Crawford, Andrew Brooks
 *
 *
 * 2.00 psc Sun May 28 1995 - Modifed for non-Dundee use.
 *
 */

#ifndef _SGDP4H_H
#define _SGDP4H_H

/*
 * Set up standard system-dependent names UNIX, LINUX, RISCOS, MSDOS, WIN32
 */

#if defined( unix )
# define UNIX
# if defined( linux ) && !defined( LINUX )
# define LINUX
# endif
#elif defined( __riscos ) && !defined( RISCOS )
# define RISCOS
#elif !defined( MSDOS ) && !defined( WIN32 ) && !defined( __CYGWIN__ )
# define MSDOS
#endif

/*
 * Include files
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <memory.h>
#ifndef _TIME_H
#define _TIME_H
#include <time.h>
#endif
#include <sys/types.h>
#ifdef UNIX
#include <unistd.h>
#endif

#ifdef SUN4
#include <memory.h>
#endif

#ifdef sun
#include <sys/time.h>    /* solaris 7 has struct timeval in here */
#include <sunmath.h>     /* for sincos() which is in libsunmath */
#endif

#ifdef linux
#include <stdint.h>
void sincos(double x, double *s, double *c); /* declared where? */
#endif

/*
 * ================= SYSTEM SPECIFIC DEFINITIONS =====================
 */

/* Use INLINE keyword when declaring inline functions */
#ifdef WIN32
#define INLINE __inline
#elif defined( MSDOS )
#define INLINE
#else
/*UNIX?*/
#define INLINE inline
#endif

/* Sun C compiler has automatic inline and doesn't understand inline keyword */
#ifdef __SUNPRO_C
#undef INLINE
#define INLINE
#define MACROS_ARE_SAFE
#endif

/* Some very common constants. */

#ifndef M_PI
#define M_PI  3.141592653589793
#endif /* MSDOS */

#ifndef PI
#define PI M_PI
#endif

#define TWOPI   (2.0*PI)    /* Optimising compiler will deal with this! */
#define PB2     (0.5*PI)
#define PI180   (PI/180.0)

#define SOLAR_DAY       (1440.0)             /* Minutes per 24 hours */
#define SIDERIAL_DAY    (23.0*60.0 + 56.0 + 4.09054/60.0)   /* Against stars */

#define EQRAD   (6378.137)                   /* Earth radius at equator, km */
#define LATCON  (1.0/298.257)                /* Latitude radius constant */
#define ECON    ((1.0-LATCON)*(1.0-LATCON))

#define JD1900 2415020.5    /* Julian day number for Jan 1st, 00:00 hours 1900 */


/*
 * =============================== MACROS ============================
 *
 *
 *  Define macro for sign transfer, double to nearest (long) integer,
 *  to square an expression (not nested), and A "safe" square, uses test
 *  to force correct sequence of evaluation when the macro is nested.
 */

/*
 * These macros are safe since they make no assignments.
 */
#define SIGN(a, b)  ((b) >= 0 ? fabs(a) : -fabs(a))
/* Coordinate conversion macros */
#define DEG(x) ((x)/PI180)
#define RAD(x) ((x)*PI180)
#define GEOC(x) (atan(ECON*tan(x))) /* Geographic to geocentric. */
#define GEOG(x) (atan(tan(x)/ECON))

/*
 * All other compilers can have static inline functions.
 * (SQR is used badly here: do_cal.c, glat2lat.c, satpos.c, vmath.h).
 */
static INLINE int       NINT(double  a) { return (int)(a > 0 ? a+0.5 : a-0.5); }
static INLINE long      NLONG(double a) { return (long)(a > 0 ? a+0.5 : a-0.5); }

static INLINE double    DSQR(double a) { return(a*a); }
static INLINE float     FSQR(float a)  { return(a*a); }
static INLINE int       ISQR(int   a)  { return(a*a); }

static INLINE double    DCUBE(double a) { return(a*a*a); }
static INLINE float     FCUBE(float a)  { return(a*a*a); }
static INLINE int       ICUBE(int   a)  { return(a*a*a); }

static INLINE double    DPOW4(double a) { a*=a; return(a*a); }
static INLINE float     FPOW4(float a)  { a*=a; return(a*a); }
static INLINE int       IPOW4(int   a)  { a*=a; return(a*a); }

static INLINE double    DMAX(double a,double b) { if (a>b) return  a; else return b; }
static INLINE float     FMAX(float a, float b)  { if (a>b) return  a; else return b; }
static INLINE int       IMAX(int   a, int   b)  { if (a>b) return  a; else return b; }

static INLINE double    DMIN(double a,double b) { if (a<b) return  a; else return b; }
static INLINE float     FMIN(float a, float b)  { if (a<b) return  a; else return b; }
static INLINE int       IMIN(int   a, int   b)  { if (a<b) return  a; else return b; }

static INLINE double    MOD2PI(double a) { a=fmod(a, TWOPI); return a < 0.0 ? a+TWOPI : a; }
static INLINE double    MOD360(double a) { a=fmod(a, 360.0); return a < 0.0 ? a+360.0 : a; }

/*
 * Unless you have higher than default optimisation the Sun compiler
 * would prefer to be told explicitly about inline functions after their
 * declaration.
 */
#if defined(__SUNPRO_C) && !defined(MACROS_ARE_SAFE)
#pragma inline_routines(NINT, NLONG, DSQR, FSQR, ISQR, DCUBE, FCUBE, ICUBE, DPOW4, FPOW4, IPOW4)
#pragma inline_routines(DMAX, FMAX, IMAX, DMIN, FMIN, IMIN, MOD2PI, MOD360, S_GEOC, S_GEOG)
#endif

/* ==================================================================== */

typedef struct orbit_s
{
    /* Add the epoch time if required. */

  int ep_year;/* Year of epoch, e.g. 94 for 1994, 100 for 2000AD */
  double ep_day;	/* Day of epoch from 00:00 Jan 1st ( = 1.0 ) */
  double rev;	/* Mean motion, revolutions per day */
  double bstar;	/* Drag term .*/
  double eqinc;	/* Equatorial inclination, radians */
  double ecc;	/* Eccentricity */
  double mnan;	/* Mean anomaly at epoch from elements, radians */
  double argp;	/* Argument of perigee, radians */
  double ascn;	/* Right ascension (ascending node), radians */
  double smjaxs;	/* Semi-major axis, km */
  double ndot2,nddot6;  /* Mean motion derivatives */
  char desig[10]; /* International designation */
  long norb;	/* Orbit number, for elements */
  int satno;	/* Satellite number. */

} orbit_t;

typedef struct xyz_s
{
    double x;
    double y;
    double z;
} xyz_t;

typedef struct kep_s
{
    double theta;     /* Angle "theta" from equatorial plane (rad) = U. */
    double ascn;      /* Right ascension (rad). */
    double eqinc;     /* Equatorial inclination (rad). */
    double radius;    /* Radius (km). */
    double rdotk;
    double rfdotk;

	/*
	 * Following are without short-term perturbations but used to
	 * speed searchs.
	 */

	double argp;	/* Argument of perigee at 'tsince' (rad). */
	double smjaxs;	/* Semi-major axis at 'tsince' (km). */
	double ecc;		/* Eccentricity at 'tsince'. */

} kep_t;

/* ================ Single or Double precision options. ================= */

#define DEFAULT_TO_SNGL 0

#if defined( SGDP4_SNGL ) || (DEFAULT_TO_SNGL && !defined( SGDP4_DBLE ))
/* Single precision option. */
typedef float real;
#ifndef SGDP4_SNGL
#define SGDP4_SNGL
#endif

#else
/* Double precision option. */
typedef double real;
#ifndef SGDP4_DBLE
#define SGDP4_DBLE
#endif

#endif  /* Single or double choice. */

/* Something silly ? */
#if defined( SGDP4_SNGL ) && defined( SGDP4_DBLE )
#error sgdp4h.h - Cannot have both single and double precision defined
#endif

/* =========== Do we have sincos() functions available or not ? ======= */
/*
We can use the normal ANSI 'C' library functions in sincos() macros, but if
we have sincos() functions they are much faster (25% under some tests). For
DOS programs we use our assembly language functions using the 80387 (and
higher) coprocessor FSINCOS instruction:

void sincos(double x, double *s, double *c);
void sincosf(float x, float *s, float *c);

For the Sun 'C' compiler there is only the system supplied double precision
version of these functions.
*/

#ifdef MACRO_SINCOS
#define sincos(x,s,c) {double sc__tmp=(x);\
                       *(s)=sin(sc__tmp);\
                       *(c)=cos(sc__tmp);}

#define SINCOS(x,s,c) {double sc__tmp=(double)(x);\
                       *(s)=(real)sin(sc__tmp);\
                       *(c)=(real)cos(sc__tmp);}

#elif !defined( sun )

/* For Microsoft C6.0 compiler, etc. */
#ifdef SGDP4_SNGL
#define SINCOS sincosf
#else
#define SINCOS sincos
#endif /* ! SGDP4_SNGL */

void sincos(double, double *, double *);
void sincosf(float, float *, float *);

#else
/* Sun 'C' compiler. */
#ifdef SGDP4_SNGL
/* Use double function and cast results to single precision. */
#define SINCOS(x,s,c) {double s__tmp, c__tmp;\
                       sincos((double)(x), &s__tmp, &c__tmp);\
                       *(s)=(real)s__tmp;\
                       *(c)=(real)c__tmp);}
#else
#define SINCOS sincos
#endif /* ! SGDP4_SNGL */
#endif /* ! MACRO_SINCOS */

/* ================= Stack space problems ? ======================== */

#if !defined( MSDOS )
/* Automatic variables, faster (?) but needs more stack space. */
#define LOCAL_REAL   real
#define LOCAL_DOUBLE double
#else
/* Static variables, slower (?) but little stack space. */
#define LOCAL_REAL   static real
#define LOCAL_DOUBLE static double
#endif

/* ======== Macro fixes for float/double in math.h type functions. ===== */

#define SIN(x)      (real)sin((double)(x))
#define COS(x)      (real)cos((double)(x))
#define SQRT(x)     (real)sqrt((double)(x))
#define FABS(x)     (real)fabs((double)(x))
#define POW(x,y)    (real)pow((double)(x), (double)(y))
#define FMOD(x,y)   (real)fmod((double)(x), (double)(y))
#define ATAN2(x,y)  (real)atan2((double)(x), (double)(y))

#ifdef SGDP4_SNGL
#define CUBE FCUBE
#define POW4 FPOW4
#else
#define CUBE DCUBE
#define POW4 DPOW4
#endif

/* SGDP4 function return values. */

#define SGDP4_ERROR     (-1)
#define SGDP4_NOT_INIT  0
#define SGDP4_ZERO_ECC  1
#define SGDP4_NEAR_SIMP 2
#define SGDP4_NEAR_NORM 3
#define SGDP4_DEEP_NORM 4
#define SGDP4_DEEP_RESN 5
#define SGDP4_DEEP_SYNC 6

//#include "satutl.h"

/* ======================= Function prototypes ====================== */

#ifdef __cplusplus
extern "C" {
#endif

/** deep.c **/

int SGDP4_dpinit(double epoch, real omegao, real xnodeo, real xmo,
                 real orb_eo, real orb_xincl, real aodp, double xmdot,
                 real omgdot, real xnodot, double xnodp);

int SGDP4_dpsec(double *xll, real *omgasm, real *xnodes, real *em,
                real *xinc, double *xn, double tsince);

int SGDP4_dpper(real *em, real *xinc, real *omgasm, real *xnodes,
                double *xll, double tsince);

/** sgdp4.c **/

int init_sgdp4(orbit_t *orb);
int sgdp4(double tsince, int withvel, kep_t *kep);
void kep2xyz(kep_t *K, xyz_t *pos, xyz_t *vel);
int satpos_xyz(double jd, xyz_t *pos, xyz_t *vel);

#ifdef __cplusplus
}
#endif

#endif /* !_SGDP4H_H */
