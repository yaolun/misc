#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG 0
#define MAX_ENTRIES 10
#define RADIAN 0.0174532925200
#define PI 3.14159265359
#define MAXCHARS 80
#define LESCHARS 20

typedef struct {
     char *name;
     char *val;
   } entry;

/*  vlsr.c is a program to calculate the velocity component of the
 *        observer with respect to the Local Standard of Rest.
 *
 *  Author:  Edward M. Murphy
 *  Written:  1998 February 17
 *  
 */

main(argc, argv)
int argc;
char *argv[];
{

/*  Define variable and arrays */
register int x,m=0;
entry entries[MAX_ENTRIES];
char *year_in, *month_in, *day_in, *uttime_in, *right_ascension,*declination;
char *equinox, *latitude, *longitude, *elevation, *earthframe;
char *makeword();
char *fmakeword();
char epochv[2];
int i,j,k,cl;
double ra, dec, year, month, day, ut, lat, lon, elev, doy, rain, decin, utin;
double xlst, vsun, vmon, vobs,v1,v2,v3,mxv,mxh,mxdv,mxdh, epoch, uttime;
double ra_sys, dec_sys, vel_sys;
float epch1;

static char daytab[2][13] = {
  {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
  {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};
 
int day_of_year(year, month, day)
int year, month, day;
{
int i, leap;
 
     leap = year%4 == 0 && year%100 != 0 || year%400 == 0;
     for (i=1; i < month; i++)
        day += daytab[leap][i];
     return day;
};

// /*  Set up output HTML form */
// printf("Content-type: text/html%c%c",10,10);

// if(strcmp(getenv("REQUEST_METHOD"),"POST")) {
//         printf("This script should be referenced with a METHOD of POST.\n");
//         printf("If you don't understand this, see this ");
//         printf("<A HREF=\"http://www.ncsa.uiuc.edu/SDG/Software/Mosaic/Docs/fill-out-forms/overview.html\">forms overview</A>.%c",10);
//         exit(1);
//     }
// if(strcmp(getenv("CONTENT_TYPE"),"application/x-www-form-urlencoded")) {
//         printf("This script can only be used to decode form results. \n");
//         exit(1);
//     }
// printf("<HTML><HEAD></HEAD><BODY>\n");
// printf("<TITLE>Ed Murphy's VLSR Calculator</TITLE>\n");
// printf("<h1>V<sub>LSR</sub> Calculator Version 1.0</h1>\n");

// /*  Determine the length of the input from Mosaic */    
// cl = atoi(getenv("CONTENT_LENGTH"));

// /*  Read in the input one character at a time.  The fields are delimeted
//  *  by &'s. */
// m=0;
// for(x=0;cl && (!feof(stdin));x++) {
//         m=x;
//         entries[x].val = fmakeword(stdin,'&',&cl);
//         plustospace(entries[x].val);
//         nospace(entries[x].val);
//         entries[x].name = makeword(entries[x].val,'=');
//     }

// if (DEBUG) printf("The form contained the following name/value pairs:<p>%c",10);
// if (DEBUG) printf("<ul>%c",10);

// /*  Now, make sure that entries[0] points to RA, [1] to dec, etc. */
// year_in = month_in = day_in = uttime_in = right_ascension = NULL;
// declination = equinox = latitude = longitude = elevation = earthframe = NULL;
// for(x=0; x <= m; x++) {
//      if (DEBUG) printf("<li> <code>\"%s\" = \"%s\"</code>%c",
//           entries[x].name, entries[x].val, 10);
//      if (strcmp(entries[x].name, "year") == 0) 
// 	  year_in = entries[x].val;
//      if (strcmp(entries[x].name, "month") == 0) 
// 	  month_in = entries[x].val;
//      if (strcmp(entries[x].name, "day") == 0) 
// 	  day_in = entries[x].val;
//      if (strcmp(entries[x].name, "uttime") == 0) 
// 	  uttime_in = entries[x].val;
//      if (strcmp(entries[x].name, "right_ascension") == 0) 
// 	  right_ascension = entries[x].val;
//      if (strcmp(entries[x].name, "declination") == 0) 
// 	  declination = entries[x].val;
//      if (strcmp(entries[x].name, "equinox") == 0) 
// 	  equinox = entries[x].val;
//      if (strcmp(entries[x].name, "longitude") == 0) 
// 	  longitude = entries[x].val;
//      if (strcmp(entries[x].name, "latitude") == 0) 
// 	  latitude = entries[x].val;
//      if (strcmp(entries[x].name, "elevation") == 0) 
// 	  elevation = entries[x].val;
//      if (strcmp(entries[x].name, "earthframe") == 0) 
// 	  earthframe = entries[x].val;
//     }
//     if (DEBUG) printf("</ul>%c",10);

// year=atof(year_in);
// if (strncmp(month_in,"Jan",3)==0) month=1.0;
//    else if (strncmp(month_in,"Feb",3)==0) month=2.0;
//    else if (strncmp(month_in,"Mar",3)==0) month=3.0;
//    else if (strncmp(month_in,"Apr",3)==0) month=4.0;
//    else if (strncmp(month_in,"May",3)==0) month=5.0;
//    else if (strncmp(month_in,"Jun",3)==0) month=6.0;
//    else if (strncmp(month_in,"Jul",3)==0) month=7.0;
//    else if (strncmp(month_in,"Aug",3)==0) month=8.0;
//    else if (strncmp(month_in,"Sep",3)==0) month=9.0;
//    else if (strncmp(month_in,"Oct",3)==0) month=10.0;
//    else if (strncmp(month_in,"Nov",3)==0) month=11.0;
//    else if (strncmp(month_in,"Dec",3)==0) month=12.0;
//    else if (strncmp(month_in,"Mon",3)==0) {
//         printf("\nYou did not enter a Month, using January as default\n");
//         printf("<BR>");
//         month=1.0;
//         }
//    else month=1.0;
// day=atof(day_in);
// utin=atof(uttime_in); 
// rain=atof(right_ascension);
// decin=atof(declination);
// lat=atof(latitude);
// lon=atof(longitude);
// elev=atof(elevation);
// sscanf(equinox,"%1s%f",epochv,&epch1);
// epoch=(double) epch1;

// /*  Convert ra and dec to degrees. */
// ra=0.0;
// dec=0.0;
// parnum(rain,&ra);
// ra*=15.0;
// parnum(decin,&dec);
// parnum(utin,&uttime);
 
// // doy=(double) day_of_year((int) year, (int) month, (int) day);

// /*  Run dop for IAU defined standard of rest */
ra_sys=18.0;
dec_sys=30.0;
vel_sys=20.0;
// if (DEBUG) printf("DOP INPUT: %10.6f %10.6f %5.1f %5.2f %10.6f\n",ra,dec,
//                    year,doy,uttime);
// if (DEBUG) printf("DOP INPUT: %10.6f %10.6f %10.6f\n",lon,lat,elev);

/* ra and dec are in degree */ 
ra = 67.89196135;
dec = 18.13468391;
epoch = 2000;
year = 2016;
month = 11;
day = 11;
doy=(double) day_of_year((int) year, (int) month, (int) day);
uttime = 6;
/* lat and lon in degree */
lat = 34.58;
lon = 118.1;
elev = 12500;

dop_ed(ra,dec,epoch,year,doy,uttime,lat,lon,elev,&xlst,&vsun,&vmon,&vobs,&v1,
       ra_sys,dec_sys,vel_sys);

/* Geocentric */
printf("Geocentric\n");
printf("Standard IAU Solar Motion V_LSR =%+10.2f\n",v1-vobs); 
/* Topocentric */
printf("Topocentric\n");
printf("Standard IAU Solar Motion V_LSR =%+10.2f\n",v1); 

/*  Run dop for Mihalas and Binney standard */
ra_sys=17.765618;  /* l=53.0, b=25.0 in epoch 1900 coords */
dec_sys=28.028649;
vel_sys=16.5;
dop_ed(ra,dec,epoch,year,doy,uttime,lat,lon,elev,&xlst,&vsun,&vmon,&vobs,&v1,
       ra_sys,dec_sys,vel_sys);

/* Geocentric */
printf("Geocentric\n");
printf("Mihalas & Binney Solar Motion V_LSR =%+10.2f\n",v1-vobs); 
/* Topocentric */
printf("Topocentric\n");
printf("Mihalas & Binney Solar Motion V_LSR =%+10.2f\n",v1);

// if (earthframe[0]=='g') 
//      printf("<H3>Geocentric</H3>\n");
//    else
//      printf("<H3>Topocentric</H3>\n");


// if (earthframe[0]=='g') 
//      printf("<H3>Standard IAU Solar Motion V<sub>LSR</sub>=%+10.2f</H3>\n",v1-vobs); 
//  else
//      printf("<H3>Standard IAU Solar Motion V<sub>LSR</sub>=%+10.2f</H3>\n",v1); 

//      printf("<BR>\n");

/*  Run dop for Mihalas and Binney standard */
// ra_sys=17.765618;  /* l=53.0, b=25.0 in epoch 1900 coords */
// dec_sys=28.028649;
// vel_sys=16.5;
// dop_ed(ra,dec,epoch,year,doy,uttime,lat,lon,elev,&xlst,&vsun,&vmon,&vobs,&v1,
//        ra_sys,dec_sys,vel_sys);

// if (earthframe[0]=='g') 
//      printf("<H3>Mihalas & Binney Solar Motion V<sub>LSR</sub>=%+10.2f</H3>\n",v1-vobs); 
//  else
//      printf("<H3>Mihalas & Binney Solar Motion V<sub>LSR</sub>=%+10.2f</H3>\n",v1); 


// /*  End HTML stuff */   
//     printf("<HR>\n");
//     printf("The algorithm for calculating V<sub>LSR</sub> was written by\n");
//     printf("M. A. Gordon and can be found in <EM> Methods of Experimental Physics,\n");
//     printf("Volume 12, Part C: Radio Observations </EM>, Ed. M. L. Meeks.\n");
//     printf("<BR> The IAU Standard Solar Motion assumes:\n");
//     printf("<UL><LI> V = 20 km s<sup>-1</sup> toward \n");
//     printf("    <LI> RA = 18<sup>h</sup> (Equinox B1900.0)\n");
//     printf("    <LI> DEC = 30&deg; (Equinox B1900.0)\n");
//     printf("</UL>\n");
//     printf("The Mihalas and Binney Solar Motion assumes:\n");
//     printf("<UL><LI> V = 16.5 km s<sup>-1</sup> toward \n");
//     printf("    <LI> <em>l</em> = 53&deg;\n");
//     printf("    <LI> <em>b</em> = 25&deg;\n");
//     printf("</UL>\n");

//     printf("The B in the epoch field refers to the FK4 system (Besselian), J\n");
//     printf("refers to the FK5 system (Julian).\n");
//     printf("The quoted accuracy is 0.02 km s<sup>-1</sup>.\n");


//     printf("<ADDRESS><A HREF=\"mailto:emurphy@pha.jhu.edu\">\n");
//     printf("emurphy@pha.jhu.edu</A></ADDRESS>\n");
//     printf("</BODY></HTML>\n");
}/*  End main */

/*  Procedure day_of_year calculates the day of the year. */
 
// static char daytab[2][13] = {
//   {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
//   {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
// };
 
// int day_of_year(year, month, day)
// int year, month, day;
// {
// int i, leap;
 
//      leap = year%4 == 0 && year%100 != 0 || year%400 == 0;
//      for (i=1; i < month; i++)
//         day += daytab[leap][i];
//      return day;
// }












