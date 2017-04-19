#include<stdio.h>
#include<math.h>
#include "common.h"
#include<assert.h>
#include<time.h>

void chem_potential()
{ double bias_main;
  nhatav_main = nhatav;
  k_main = k;

  bias_freenerg();
  bias_main = bias_av;
  
  ub_nhatav = ener_bias(nhatav);
  printf("\nUb_nhatav: %lf\n",ub_nhatav);

  therm_integral = bias_onoff();
  printf("\nThermal integral: %lf\n",therm_integral);

  A_nhatavg = Ab_nhatavg + therm_integral - ub_nhatav;
  printf("\nFree energy at nhatav = %lf   is %lf\n",nhatav_main,A_nhatavg);
  printf("\nbias: %lf\n",bias_main);
  FILE *writeA_nav;
  writeA_nav = fopen("A_nav_vs_nav.txt","a");
  fprintf(writeA_nav,"%lf %lf \n",nhatav_main,A_nhatavg);
  fclose(writeA_nav);
    
  
}

double bias_freenerg()                 //Computes Ab_nhatav
{
  /*FILE *writeAb;
  writeAb = fopen("Ab_nt_vs_nt.txt","a");

  for(i=0;i<5;i++)
  {
       p_nhat[i] = count[i]/nattempt;
       //    printf("count: %lf\nnatt:%lf\n",count[i], nattempt);
       //printf("Probability of nhat particles in the box: %lf\n",p_nhat[i]);

      Ab_nhat[i] = -1*(1.0/beta)*log(p_nhat[i]);
      //printf("Ab_nhat: %lf\n",Ab_nhat[i]);

      fprintf(writeAb,"%lf %lf \n",n[i],Ab_nhat[i]);

  }

  fclose(writeAb);
  */

  
  //nt=80
  Ab_nhatavg =  0.009983429*pow(nhatav,2.0) -1.613203171*nhatav + 65.1782744;



  printf("Abnhatav: %lf\n", Ab_nhatavg);
  // printf("nhatav: %lf\n", nhatav);


}

double bias_onoff()
{
  double ll, ul, w1, w2, w3, w4, w5, w6,w7, w8, w9, w10, x1, x2, x3, x4, x5, x6,x7, x8, x9, x10, arg1, arg2, arg3, arg4, arg5, arg6,arg7, arg8, arg9, arg10, integral, integral1;
  double val1;
  double l,u,ww1,ww2,ww3,ww4,ww5,ww6,xx1,xx2,xx3,xx4,xx5,xx6,a1,a2,a3,a4,a5,a6,integral2;
  ll = 0.0;
  ul = k;
  //10 point quadrature
  w1 = 0.295524224715;
  w2 = 0.295524224715;
  w3 = 0.269266719309;
  w4 = 0.269266719309;
  w5 = 0.219086362516;
  w6 = 0.219086362516;
  w7 = 0.149451349151;
  w8 = 0.149451349151;
  w9 = 0.066671344309;
  w10 = 0.066671344309; 

  x1 = -0.148874338982;
  x2 = 0.148874338982;
  x3 = -0.433395394129;
  x4 = 0.433395394129;
  x5 = -0.679409568299;
  x6 = 0.679409568299;
  x7 = -0.865063366689;
  x8 = 0.865063366689;
  x9 = -0.973906528517;
  x10 = 0.973906528517;

  arg1 = ((ul-ll)/2.0)*x1 + ((ul+ll)/2.0);
  arg2 = ((ul-ll)/2.0)*x2 + ((ul+ll)/2.0);
  arg3 = ((ul-ll)/2.0)*x3 + ((ul+ll)/2.0);
  arg4 = ((ul-ll)/2.0)*x4 + ((ul+ll)/2.0);
  arg5 = ((ul-ll)/2.0)*x5 + ((ul+ll)/2.0);
  arg6 = ((ul-ll)/2.0)*x6 + ((ul+ll)/2.0);
  arg7 = ((ul-ll)/2.0)*x7 + ((ul+ll)/2.0);
  arg8 = ((ul-ll)/2.0)*x8 + ((ul+ll)/2.0);
  arg9 = ((ul-ll)/2.0)*x9 + ((ul+ll)/2.0);
  arg10 = ((ul-ll)/2.0)*x10 + ((ul+ll)/2.0);

  //printf("arg1: %lf\targ2: %lf\targ3:%lf\targ4:%lf\targ5: %lf\targ6: %lf\n",arg1,arg2,arg3,arg4,arg5,arg6);

  integral1 =((ul-ll)/2.0)*(w1*integrand(arg1) + w2*integrand(arg2) + w3*integrand(arg3) + w4*integrand(arg4) + w5*integrand(arg5) + w6*integrand(arg6) + w7*integrand(arg7) + w8*integrand(arg8) + w9*integrand(arg9) + w10*integrand(arg10));

 
  /*
 
  //6 point quadrature
  l = 0.005;
  u = 0.02;
  ww1 = 0.360761573048;
  ww2 = 0.360761573048;
  ww3 = 0.467913934573;
  ww4 = 0.467913934573;
  ww5 = 0.171324492379;
  ww6 = 0.171324492379;

  xx1 = 0.661209386466;
  xx2 = -0.661209386466;
  xx3 = -0.238619186083;
  xx4 = 0.238619186083;
  xx5 = -0.932469514203;
  xx6 = 0.932469514203;
  
  a1 = ((u-l)/2.0)*xx1 + ((u+l)/2.0);
  a2 = ((u-l)/2.0)*xx2 + ((u+l)/2.0);
  a3 = ((u-l)/2.0)*xx3 + ((u+l)/2.0);
  a4 = ((u-l)/2.0)*xx4 + ((u+l)/2.0);
  a5 = ((u-l)/2.0)*xx5 + ((u+l)/2.0);
  a6 = ((u-l)/2.0)*xx6 + ((u+l)/2.0);

  //  printf("arg1: %lf\targ2: %lf\targ3:%lf\targ4:%lf\targ5: %lf\targ6: %lf\n",arg1,arg2,arg3,arg4,arg5,arg6);

integral2=((u-l)/2.0)*(ww1*integrand(a1) + ww2*integrand(a2) + ww3*integrand(a3) + ww4*integrand(a4) + ww5*integrand(a5) + ww6*integrand(a6));
  */
  //integral = integral1 + integral2;
  integral = integral1;
  /*//4 point quadrature  
  w1 = 0.652145;
  w2 = 0.652145;
  w3 = 0.347855;
  w4 = 0.347855;

  x1 = -0.339981;
  x2 = 0.339981;
  x3 = -0.861136;
  x4 = 0.861136;

  arg1 = ((ul-ll)/2.0)*x1 + ((ul+ll)/2.0);
  arg2 = ((ul-ll)/2.0)*x2 + ((ul+ll)/2.0);
  arg3 = ((ul-ll)/2.0)*x3 + ((ul+ll)/2.0);
  arg4 = ((ul-ll)/2.0)*x4 + ((ul+ll)/2.0);

  printf("arg1: %lf\targ2: %lf\targ3:%lf\targ4:%lf\n",arg1,arg2,arg3,arg4);

  integral =((ul-ll)/2.0)*(w1*integrand(arg1) + w2*integrand(arg2) + w3*integrand(arg3) + w4*integrand(arg4));
  */
  
    

  //2 point quadrature
  /*w1 = 1.0;
  w2 = 1.0;

  x1 = -0.5773502692;
  x2 = 0.5773502692;

  arg1 = ((ul-ll)/2.0)*x1 + ((ul+ll)/2.0);
  arg2 = ((ul-ll)/2.0)*x2 + ((ul+ll)/2.0);

  printf("arg1: %lf\targ2: %lf\n",arg1,arg2);
  integral =((ul-ll)/2.0)*(w1*integrand(arg1) + w2*integrand(arg2));
  */
  return integral;
}

double integrand(double k_arg)
{
  double f_k, nhatav_karg;

  k = k_arg;
  
  sys_set();
  set_box();
  mmc();
  
  f_k = bias_av;
    
  FILE *writen_k;
  writen_k = fopen("nhatav_vs_k.txt","a");
  fprintf(writen_k,"%lf %lf %lf %lf\n",nt,nhatav, k, f_k);
  fclose(writen_k);
    
  return f_k;
}
