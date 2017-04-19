#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "common.h"

main()
{ 
  prog_flag =1;
  
  //Read inputs
  read_input();
  
  //Set up the system's initial configuration
  sweep = 0;
  sys_set();
  
  //Write initial configuration
  coord_write(sweep);
    
  //Set up the biased box 
  set_box();
    
  
  //Calculate Tail Corrections
  V_corr = cor_U(rc,rho2);            //In units of eV
  printf("\nV_corr (eV): %lf\n",V_corr);
    
  P_corr = cor_P(rc,rho2);            //In units of ev/Angstrom^3
  printf("\nP_corr: %lf\n",P_corr);
    
  //Compute initial energy of the system
  V = tot_energy();
  printf("\nInitial sysytem potential: %lf\n",V);
  
  //THE MMC
  mmc();

  prog_flag = 0;
 
 //Compute the chemical potential    
  chem_potential();
  
 
      
}
