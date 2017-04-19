#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "common.h"

double tot_energy()
{
  double xi, yi, zi, tot_pot, V, eni;
  int jb;
  tot_pot = 0.0;

  for(i=1;i<=(np-1);i++)
    {
      xi = x_coord[i];
      yi = y_coord[i];
      zi = z_coord[i];

      jb = i+1;

      eni = ener_part(i, xi, yi, zi, jb);
     
      tot_pot+=eni;

    }
    
  V = tot_pot + np*V_corr;

  return V;
}

double ener_part(int I, double Xi, double Yi, double Zi, int Jb)
{
  int j, jj;
  double dx, dy, dz, pot, part_energy,r2;

  part_energy = 0.0;

  if(vlist_flg == 1)
    {
      for(jj=1; jj<=nlist[I]; jj++)
	{
	  j = list[I][jj];
	  if((j!=I) &&(j>=Jb))
	    {
	      r2 = part_dist(Xi, Yi, Zi, j);
	      pot = calc_pot(r2);
	      part_energy+=pot;
	    }
	}
      return part_energy;
    }

  else 
    {
      for(j=Jb; j<=np; j++)
	{
	  if(j!=I)
	    {
	      r2 = part_dist(Xi, Yi, Zi, j);
	      pot = calc_pot(r2);
      
	      part_energy+=pot;
	      
	    }
	}
      return part_energy;
    }
}

double calc_pot(double R2)
{
  double r2i, r6i, pair_pot;

  if(R2<rc2)
    {
      r2i = sig2/R2;
      r6i = r2i*r2i*r2i;
      
      pair_pot = eps4*(r6i*r6i - r6i);
    }
  else
    {
      pair_pot = 0.0;

    }
  return pair_pot; 

}
