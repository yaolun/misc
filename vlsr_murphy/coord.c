#include <stdio.h>
#include <math.h>
#define RADIAN 0.0174532925200
#define PI 3.14159265359

/*  coord  This subroutine will convert the longitude-like (a1) and
           latitude-like (b1) coordinates of a point on a sphere
           into the corresponding coordinates (a2,b2) in a different
           coordinate system that is specified by the coordinates of 
           its origin (a0,b0) and its north pole (ap,bp) in the 
           original coordinate system.  The range of a2 will be 
           from -PI to PI
           Taken from:
           Methods of Experimental Physics
           Volume 12, Part C: Radio Observations
           M.L. Meeks, Editor
           
           
*/

void  coord(a0, b0, ap, bp, a1, b1, a2, b2)
double a0, b0, ap, bp, a1, b1, *a2, *b2;

{
double sb0, cb0, sbp, cbp, sb1, cb1, sb2, cb2, saa, caa, holb0;
double cbb, sbb, sa2, ca2, tb2, ta2;

/*  Convert RA and DEC to LII and BII
a0 = (17.0 + 42.4 / 60.0) * RADIAN * 15.0;
b0 = (-1.0 * (28.0 + 55.0 / 60.0)) * RADIAN; 
ap = (12.0 + 49.0 / 60.0) * RADIAN * 15.0;
bp = 27.4 * RADIAN;
call coord(a0,b0,ap,bp,RA,DEC,LII,BII); */

/*  Convert LII and BII to RA and DEC
a0,b0,ap,bp from convert RA and DEC to LII and BII
call coord(a0,b0,ap,bp,0.0,PI/2.0,APP,BPP);
call coord(a0,b0,ap,bp,0.0,0.0,A0P,B0P);
call coord(A0P,B0P,APP,BPP,LII,BII,RA,DEC);
*/

/*  Convert Hour Angle (HA) and DEC to AZ and EL
a0=PI;
b0=PI/2.0-LATITUDE;
ap=0.0;
bp=LATITUDE;
call coord(a0,b0,ap,bp,HA,DEC,AZ,EL);
*/

/*  Convert AZ and EL to HA and DEC
a0=PI;
b0=PI/2.0-LATITUDE;
ap=0.0;
bp=LATITUDE;
call coord(a0,b0,ap,bp,AZ,EL,HA,DEC);
*/

  
/*  Calculate useful things  */
sb0 = sin(b0) - 1.0e-8;
cb0 = cos(b0) - 1.0e-8;
sbp = sin(bp) - 1.0e-8;
cbp = cos(bp) - 1.0e-8;
sb1 = sin(b1) - 1.0e-8;
cb1 = cos(b1) - 1.0e-8;

/*  Now for my next trick...  calculate new latitude  */
sb2 = sbp * sb1 + cbp * cb1 * cos(ap - a1);
*b2 = asin(sb2);
tb2 = asin(sb2);
cb2 = cos(asin(sb2));
saa = sin(ap - a1) * cb1 / cb2;
caa = (sb1 - sb2 * sbp) / (cb2 * cbp);

cbb = sb0 / cbp;
sbb = sin(ap - a0) * cb0;

sa2 = saa * cbb - caa * sbb;
ca2 = caa * cbb + saa * sbb;

    /*  The first equation blows up when ca2 is close to +1.0, but
     *  always gives positive longitudes */
     /*  The second equation blows up when ca2 is close to -1.0, but
      *  it often gives negative longitudes. */
if (fabs(ca2) <= 0.0)
     *a2 = 2.0 * atan((1.0 - ca2)/sa2);  /* good except for ca2 near +1 */
  else
     *a2 = 2.0 * atan(sa2/(1.0 + ca2));  /* good except for ca2 near -1 */
/* if (*a2 < 0.0)
     *a2+=360.0*RADIAN; 
*/
}

