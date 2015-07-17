#include <stdio.h>
#include <math.h>
#define RADIAN 0.0174532925200
#define PI 3.14159265359


/*   precess  Calculates the correction in RA and DEC to be added when 
 *            precessing coordinates.  All angles given in radians. 
 *            PRECESS calculates the corrections to be added to the
 *            mean coordinates for epoch to give the apparent 
 *            coordinates for epoch1.  PRECESS also calculates the 
 *            equation of the equinoxes (DC in minutes of time)
 *            which may be added to the mean sideraeal time to
 *            give the apparent sidereal time.  deld and delra 
 *            (the corrections) contain corrections for precession,
 *            annual abberation, and some terms of nutation.  If
 *            RA and DEC are for the mean epoch (i.e. halfway between
 *            epoch and epoch teh the precission of delra and deld is
 *            about 2 arcseconds.  If RA and DEC are either of the
 *            endpoints, the precision is somewhat worse.
 *            j is the day of the year, epoch is the epoch of the
 *            given ra and dec.  Epoch1 is the year that the coord.
 *            should be precessed to.
 *
 *
 *         Taken from:
 *         Methods of Experimental Physics
 *         Volume 12, Part C: Radio Observations
 *         M.L. Meeks, Editor
 */

void precess(j, epoch1, epoch, ra, dec, pra, pdec)
double j, epoch1, epoch, ra, dec, *pra, *pdec;
{
double snd, csd, tnd, csr, snr, al, jul, jep, j00, t0, t,
       zeta0, z, theta, am, an, alam, snl, csl, omega, arg, dlong, doblq;
double deld, delra, delep;

snd = sin(dec);
csd = cos(dec);
tnd = snd / csd;

csr = cos(ra);
snr = sin(ra);

/*  al is an approximate day number (i.e. the number of days since
       January 0 of the year epoch1.  */
al=j;

/*  t0 is the time from 1900 to epoch (centuries)  */ 
t0 = (epoch-1900.0)/100.0;

/*  t is the time from epoch to epoch1 (centuries)  */
t = (epoch1-epoch)/100.0 + al/(365.2421988*100.0);

/*  zeta0, z, and theta are precessional angles  */
zeta0 = (2304.250 + 1.396 * t0) * t + 0.302 * t * t + 0.018 * t * t  * t;
z = zeta0 + 0.791 * t * t;
theta = (2004.682 - 0.853 * t0) * t - 0.426 * t * t - 0.042 * t * t * t;

/*  am and an are the precessional numbers  */
am = (zeta0+z) * 4.848136811e-6;
an = theta * 4.848136811E-6;

/*  alam is an approximate mean longitude for the sun  */
alam = (0.985647 * al + 278.5) * 0.0174532925;
snl = sin(alam);
csl = cos(alam);

/*  delra and deld are the annual abberation terms in radians  */
/*  0.91745051 is the cos(obliquity of ecliptic)  */

delra = -9.92413605e-5 * (snl * snr + 0.91745051 * csl *csr) / csd
        + am + an * snr * tnd;
deld = -9.92413605e-5 * (snl * csr * snd - 0.91745051 * csl * snr * snd
        + 0.39784993 * csl * csd) + an * csr;

/*  omega is the angle of the first term of nutation in degrees  */
omega = 259.183275 - 1934.142 * (t0 + t);
arg = omega * RADIAN;

/*  dlong is the nutation in longitude in radians  */
dlong = -8.3597e-5 * sin(arg);
/*  doblq is the nutation in obliquity  */
doblq = 4.4678e-5 * cos(arg);

/*  Add in nutation terms  */
delra += dlong * (0.91745051 + 0.39784993 * snr * tnd) - csr * tnd * doblq;
deld += 0.39784993 * csr * dlong + snr * doblq;

/*  Correct ra and dec  */
*pra = ra + delra;
*pdec = dec + deld;

}

