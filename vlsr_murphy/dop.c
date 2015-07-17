#include <stdio.h>
#include <math.h>
#define RADIAN 0.0174532925200
#define PI 3.14159265359

/*  dop  This subroutine will calculate the motion of the earth 
 *       with respect to the Local Standard of Rest.  This program 
 *       omits the planetary perturbations on the earth's orbit.  
 *       These amount to about 0.003 km/s and are thought to be the
 *       largest contribution to the error in the velocity.
 *       The ra and dec must be specified in degress.
 *       The location of the observer is specified by the latitude
 *       (alat), the geodetic longitude in degrees (olong), and
 *       the elevation in meters above sea level (elev).
 *       The program gives as output the local mean sidereal time
 *       in days (xlst), the velocity component in km/s of the
 *       sun's motion with respect to the LSR as projected onto
 *       the line of sight to the source (vsun), and the total
 *       velocity component (v1).  Positive velocity corresponds
 *       to increasing distance between source and observer.
 *
 *       ra_sys and dec_sys are the direction in which the Sun is 
 *       moving at velocity vel_sys.  ra_sys is in hours and dec_sys
 *       is in degrees.  For LSR, ra_sys=18.0, dec_sys=30.0 and 
 *       vel_sys=20.0;
 *       
 *         Taken from:
 *         Methods of Experimental Physics
 *         Volume 12, Part C: Radio Observations
 *         M.L. Meeks, Editor
 *
 *   Author of c version:  Edward M. Murphy
 *   Modified:  1995 March 14
 */

void dop_ed(ra1,dec,epoch,yr,day,ut,alat,olong,elev,xlst,vsun,vmon,vobs,v1,ra_sys,dec_sys,vel_sys)
double ra1,dec,epoch,yr,day,ut,alat,olong,elev,*xlst,*vsun,*vmon,*vobs,*v1,ra_sys,dec_sys,vel_sys;

{
double aaa,dd,x0,y0,z0,ra,cc,cs,s,cat,wlong,du,tu,utda,smd;
double t,start,c1,gst,r5ho,vrho,dlat,am,mpi,e,a1,vs,xlsta;
double xlam,alam,aa,an,hop,v,omga,omgar,amon,gamp,pim,em,olamm;
double aim,amm,vsm,alamm,anm,aam,hopm,junk1,junk2;
double dela,deldd,junk,rho,along,beta,algm,betam;
int julda();

/*  The following deal with Solar motion towards RA=AAA, DEC=DD */
aaa=18.0*PI/12.0;
dd=30.0*RADIAN;

/*  Precess this direction to date */
precess(day,yr,1900.0,aaa,dd,&dela,&deldd);
aaa=dela;
dd=deldd;

/*  This velocity is converted to cartesian coordinates. */
x0=vel_sys*cos(aaa)*cos(dd);
y0=vel_sys*sin(aaa)*cos(dd);
z0=vel_sys*sin(dd);

/*  Convert ra and dec to radians */
ra=ra1*RADIAN;
ra1/=360.0;
dec*=RADIAN;

precess(0.0,yr,epoch,ra,dec,&junk1,&junk2);
ra=junk1;
ra1=ra/RADIAN/360.0;
dec=junk2;

/*  cc,cs and s are the direction cosines corresponding to ra and dec */
cc=cos(dec)*cos(ra);
cs=cos(dec)*sin(ra);
s=sin(dec);

/*  vsun is the projection onto the line of sight to the star of the
    Sun's motion with respect to the Local Standard of Rest (km/sec) */
*vsun=-1.0*x0*cc-y0*cs-z0*s;

/*  Coord. of observer, lat (radians), and long (revs=days) */
cat=alat*RADIAN;
wlong=olong/360.0;

/*  The following calculations deal with time.
 *  The epoch is 1900 January 0.5 UT = Julian Day 2415020.0
 *  du is the time from the epoch to Jan 0.0 of the current year (days)
 */
du=((double) julda((int) yr)-2415020.0)-0.500;

/*  tu is du converted to Julian centuries */
tu=du/36525.0;

/*  utda is the GMT from Jan 0.0 to the present (days) */
utda=day+ut/24.0;

/*  smd is the time from the epoch to the present (days) */
smd=du+utda;

/*  t is smd converted to Julian centuries */
t=smd/36525.0;

/*  start is the Greenwich mean sidereal time on Jan 0.0 (days)
 *  (the extra 129.174 secs corresponds to the 0.7 century subracted
 *  from tu.  The precision is thereby improved. */
start=(6.0+38.0/60.0+(45.836+129.1794+8640184.542*(tu-0.7)+0.0929*tu*tu)/3600.0)/24.0;

/*  c1 is the conversion factor from solar time to sidereal time */
c1=0.997269566414;

/*  GST is the Greenwich Mean Sidereal Time (days) */
gst=start+utda/c1;

/*  xlst is the Local Mean Sidereal Time (from jan 0) (days) */
xlsta=gst-wlong;
xlsta=modf(xlsta, &junk);

/*  The following deal with the observer's motion with respect
 *  to the earth's center.  rho is the radius vector from the
 *  earth's center to the observer (meters) */
rho=6378160.00*(0.998327073+0.001676438*cos(2.0*cat)-0.00000351*cos(4.0*cat)+0.000000008*cos(6.0*cat))+elev;

/*  vrho is the corresponding circular velocity (meters/sidereal day) */
vrho=2.0*PI*rho;

/*  converted to km/s */
vrho/=(24000.0*3600.0*c1);

/*  reduction of geodetic latitude to geocentric latitude (arcsec) */
dlat=-1.0*(11.0*60.0+32.7430)*sin(2.0*cat)+1.1633*sin(4.0*cat)-0.0026*sin(6.0*cat);

/*  convert cat to geocentric latitude (radians) */
cat+=dlat*RADIAN/3600.0;

/*  vobs is the projection onto the line of sight to the star of the
 *  velocity of the observer with respect to the earth's center (km/sec)
 */
*vobs=vrho*cos(cat)*cos(dec)*sin(2.0*PI*(xlsta-ra1));

/*  The following calculations deal with the earth's orbit about the
 *  Sun.  AM is the mean anomaly of the earth's orbit (radians) */
am=(358.47583+0.9856002670*smd-0.000150*t*t-0.000003*t*t*t)*RADIAN;

/* mpi is the mean longitude of perihelion (radians) */
mpi=(101.22083+0.0000470684*smd+0.000453*t*t+0.000003*t*t*t)*RADIAN;

/* e is the eccentricity of the orbit (dimensionless) */
e=0.01675104-0.00004180*t-0.000000126*t*t;

/* a1 is the mean obliquity of the ecliptic (radians) */
a1=(23.452294-0.0130125*t-0.00000164*t*t+0.000000503*t*t*t)*RADIAN;

/*  VS is the true anomaly (approx formula) (radians)
 *  (equation of the center) */
vs=am+(2.0*e-0.250*e*e*e)*sin(am)+1.25*e*e*sin(2.0*am)+13.0/12.0*e*e*e*sin(3.0*am);

/*  xlam is the true longitude of the earth as seen by the sun (radians)*/
xlam=mpi+vs;

/*  alam is the true longitude of the sun as seen from the earth (radians) */
alam=xlam+PI;

/*  beta is the latitude of the star (radians)
 *  along is the longitude of the star (radians) */
coord(0.0,0.0,-1.0*PI/2.0,PI/2.0-a1,ra,dec,&along,&beta);

/*  aa is the semi-major axis of the earth's orbit */
aa=149598500.0;

/*  an is the mean angular rate of the earth about the Sun (radians/day) */
an=2.0*PI/365.2564;

/*  hop is h/p from start= the component of the earth's velocity 
 *  perpendicular to the radius vector (km/day) */
hop=an*aa/sqrt(1.0-e*e);

/*  converted to km/s */
hop/=86400.0;

/* V is the projection onto the line of sight to the star of the velocity
 * of the earth-moon barycenter with respect to the Sun (km/s) */
v=-1.0*hop*cos(beta)*(sin(alam-along)-e*sin(mpi-along));

/*  The following calculations deal with the moon's orbit around the
 *  earth moon barycenter
 *  Omga is the longitude of the mean ascending node of the lunar orbit
 *  degrees */
omga=259.183275-0.0529539222*smd+0.002078*t*t+0.000002*t*t*t;

/*  omgar is omga in radians */
omgar=omga*RADIAN;

/*  amon is omga plus the mean lunar longitude of the moon (degrees) */
amon=270.434164+13.1763965268*smd-0.001133*t*t+0.0000019*t*t*t;

/*  gamp is omga plus the lunar long of the lunar perigee (degrees) */
gamp=334.329556+0.1114040803*smd-0.010325*t*t-0.000012*t*t*t;

/*  pim is the mean lunar long of lunar perigee (to radians) */
pim=(gamp-omga)*RADIAN;

/*  em is the eccentricity of the moon's orbit */
em=0.054900489;

/*  olamm is the mean lunar longitude of the moon (to radians) */
olamm=(amon-omga)*RADIAN;

/*  AIM is the inclination of the lunar orbit to the ecliptic (radians) */
aim=5.1453964*RADIAN;

/*  amm is the approx mean anomaly (radians), it is approx because
    pim should be the true rather than mean lunar long. of perigee */
amm=olamm-pim;

/*  vsm is the true anomaly (approx form) (radians)
 *  (equation of the center) */
vsm=amm+(2.0*em-0.25*em*em*em)*sin(amm)+1.25*em*em*sin(2.0*amm)+13.0/12.0*em*em*em*sin(3.0*amm);

/* alamm is the true lunar longitude of the moon (radians) */
alamm=pim+vsm;

/* anm is the mean angular rate of the the lunar rotation (rad/day) */
anm=2.0*PI/27.321661;

/* aam is the semi major axis of the lunar orbit (km) */
aam=60.2665*6378.388;

/* betam is the lunar lat of the star (radians)
 * algm is the lunar long of the star (radians) */
coord(omgar,0.0,omgar-PI/2.0,PI/2.0-aim,along,beta,&algm,&betam);

/* hopm is h/p from start = the component of the lunar velocity 
   pependicular to the radius vector. (km/day). */
hopm=anm*aam/sqrt(1.0-em*em);

/* converted to (km/sec) */
hopm=hopm/86400.0;

/* vmon is the projection onto the line of sight to the star of the 
 * velocity of the earth's center with respect to the earth moon
 * barycenter (km/sec).  The 81.30 is the ratio of the earth's
 * mass to the moon's mass */
*vmon=-1.0*hopm/81.30*cos(betam)*(sin(alamm-algm)-em*sin(pim-algm));

*v1=v+*vsun+*vmon+*vobs;
*xlst=xlsta;

} 


/*  julda   computes the julian day number at 12 hrs ut on January 0
 *          of the year nyr (Gregorian calendar).
 *          For example, julda= 2439856 for nyr=1968
 *
 *         Taken from:
 *         Methods of Experimental Physics
 *         Volume 12, Part C: Radio Observations
 *         M.L. Meeks, Editor
 *
 *   Author of c version:  Edward M. Murphy
 *   Modified:  1994 Spetember 29
 */

int julda(nyr)
int nyr;

{
int nyrm1,ic,jul;

nyrm1=nyr-1;
ic=nyrm1/100;
jul=1721425+365*nyrm1+nyrm1/4-ic+ic/4;
return jul;
}


						       


