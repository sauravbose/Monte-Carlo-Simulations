#include<stdio.h>
#include<math.h>
#include "common.h"

double cor_U(double R, double Rho)
{ 

    double ri3,coru;
    ri3 = sig3/(R*R*R);
    coru = 2*M_PI*eps4*(Rho*sig3)*((ri3*ri3*ri3)/9-ri3/3);
    return coru;
  
}

double cor_P(double R, double Rho)
{
    double ri3, corp;
    ri3 = sig3/(R*R*R);
    corp = 4*M_PI*eps4*(Rho*Rho)*sig3*((2*ri3*ri3*ri3)/9-ri3/3);
    return corp;
}
