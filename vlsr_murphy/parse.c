#include <stdio.h>
#include <math.h>

#define RADIAN 0.0174532925200
#define PI 3.14159265359

/*  parse.c  is a subroutine to divide the number ddmmss.ss into
             its component parts dd mm ss.s.  All numbers are double.
 *           Each of the returned numbers will have the same sign as
 *           ddmmss.ss.
 *
 *  Author:  Edward M. Murphy
 *  Written:  1994 April 10
 *  
 */

void parse(dec,dd,mm,ss)
double dec,*dd,*mm,*ss;

{

/*  Define variable and arrays */
int i,j,k;
double hold1;

hold1=modf(dec/10000.0, dd);
hold1=modf((dec-*dd*10000.0)/100.0, mm);
*ss=dec-*dd*10000.0-*mm*100.0;
if ((fabs(*mm) >= 60.0) || (fabs(*ss) >= 60.0))
     printf("Error in parse: deg=%10.5f  min=%10.5f  sec=%10.5f\n", 
           *dd,*mm,*ss); 
}



/*  parnum.c  is a subroutine to convert the number ddmmss.ss into
 *            dd.ddddd.  It calls parse to split up ddmmss.s.
 *    
 *
 *  Author:  Edward M. Murphy
 *  Written:  1994 April 10
 *  
 */

void parnum(ddmmss,decd)
double ddmmss,*decd;

{

/*  Define variable and arrays */
int i,j,k;
double dd,mm,ss;

parse(ddmmss,&dd, &mm, &ss);
*decd=dd+mm/60.0+ss/3600.0;

}

/*  inparnum.c  is a subroutine to convert the number dd.ddddd into
 *            ddmmss.
 *    
 *
 *  Author:  Edward M. Murphy
 *  Written:  1994 April 10
 *  
 */

void inparnum(decd,ddmmss)
double decd,*ddmmss;

{

/*  Define variable and arrays */
int i,j,k;
double dd,mm,ss,hold1,hold2;

hold1=modf(decd, &dd);
hold1=modf((decd-dd)*60.0, &mm);
ss=(((decd-dd)*60.0)-mm)*60.0;
if ((fabs(mm) >= 60.0) || (fabs(ss) >= 60.0))
     printf("Error in inparnum: deg=%10.5f  min=%10.5f  sec=%10.5f\n", 
           dd,mm,ss); 
*ddmmss=dd*10000.0+mm*100.0+ss;
}






