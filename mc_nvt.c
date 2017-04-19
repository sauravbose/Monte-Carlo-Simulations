#include<stdio.h>
#include<math.h>
#include "common.h"
#include<assert.h>


inline long int *vector(long int n, char *s )
{
  long int *v ;
  v = (long int *)calloc( (n), sizeof(long int) ) ;
  if (!v)
    {
      printf(" Error in %s , Exiting ... \n",s) ;
      exit(1) ;
    }
  assert(v) ;
  return(v) ;
}
 

void mmc()
{
  nhatsum = 0.0;
  bias_sum = 0.0;
  Vsum = 0.0;
//*********************************************START OF MMC****************************************************//

    for(process=1;process<=2;process++)                //process 1 - equilibration; process 2 - production
      {
        if(process == 1)
          {
            nsweep = min_sweep;

            if(nsweep!=0)
              printf("\nEquilibration Starts..\nNumber of sweeps in this process: %d\n",nsweep);

          }
        else if(process == 2)
          {
            nsweep = totsweep - min_sweep;
	    
	    for(i=0;i<5;i++)
	      {
		count[i] = 0.0;
		n[i] = 0.0;
		p_nhat[i] = 0.0;
		Ab_nhat[i] = 0.0;
	      }
            
	    if(nsweep!=0)
              printf("\nProduction Starts..\nNumber of sweeps in this process: %d\n",nsweep);
            else
              printf("\nSimulation over at the equilibration stage. Increase total number of sweeps\n");
          }
	
        nattempt = 0.0;
        nacc = 0.0;
	
	n[0] = nt-2.0;
	n[1] = nt-1.0;
	n[2] = nt;
	n[3] = nt+1.0;
	n[4] = nt+2.0;
	
	move_max = adjust(nattempt,nacc,move_max);                    //To adjust move size

	printf("\nmove: %lf\n",move_max);
	
	for(sweep=1; sweep<=nsweep; sweep++)
	  {
	    if(sweep%1000 == 0)
	      printf("\nSweep: %d...\n",sweep);
	   
	    for(move=1;move<=np;move++)
	      {
		mcmove();
				
	      }//End of all moves or one sweep

	    if(process == 2)
	      { 
		nhatsum += nhattemp;
		nhatav = nhatsum/sweep;
		
		
		bias_sum += dener_bias(nhattemp);
		bias_av = bias_sum/sweep;
		
		Vsum += V;
		Vav = Vsum/sweep;
	      }
	    
	    /*	    if((process == 2)&&(prog_flag == 1))
	      {
		if(sweep%20 == 0)
		  {
		    FILE *nhat_av;
		    nhat_av = fopen("avnhat_vs_sweep.ods","a");
		    fprintf(nhat_av,"%lf %lf\n",(sweep+min_sweep),nhatav);
		    fclose(nhat_av);

		    FILE *nhat_inst;
		    nhat_inst = fopen("inst_nhat_vs_sweep.ods","a");
		    fprintf(nhat_inst,"%lf %lf \n",(sweep+min_sweep),nhat);
		    fclose(nhat_inst);
		    
		    FILE *V_inst;
                    V_inst = fopen("inst_V_vs_sweep.ods","a");
                    fprintf(V_inst,"%lf %lf \n",(sweep+min_sweep),V);
                    fclose(V_inst);

		    FILE *V_av;
                    V_av = fopen("avV_vs_sweep.ods","a");
                    fprintf(V_av,"%lf %lf\n",(sweep+min_sweep),Vav);
                    fclose(V_av);


		    
		    FILE *b_av;
                    b_av = fopen("avbias_vs_sweep.ods","a");
                    fprintf(b_av,"%lf %lf\n",(sweep+min_sweep),bias_av);
                    fclose(b_av);
                    
		  }
	      }
	    */
	    if(sweep%100==0)
	      {
		//coord_write(sweep);
		move_max = adjust(nattempt,nacc,move_max); 
	      }
	    
	  }//end of all sweeps
      }//end of process loop

//*********************************************END OF MMC****************************************************//

}

double adjust(int nattempt, int nacc, double move_max)
{ 
  if(nattempt == 0)
    {
      Att = nattempt;
      Acc = nacc;
    }
  
  else
    { 
     perc = (double)(nacc-Acc)/(nattempt-Att);
     printf("\nperc: %lf\n",perc);
      
      if(perc>0.5)
	move_max = move_max*1.05;
      if(perc<0.5)
	move_max = move_max*0.95;
      
    }
  Att = nattempt;
  Acc = nacc;
  
  return move_max;
  
}

void mcmove()
{ 
  int jmove;
  double V_old, V_new, xn,yn,zn, rmove, EXP, del;
  nattempt++;
  jmove = 1;

  r = range_random(0,1);
  o = (int)((np-1)*r)+1;
  
  flag = 5;
  if(vlist_flg == 1)
    {
      del = pow((x_coord[o]-xv_coord[o]),2.0) +  pow((y_coord[o]-yv_coord[o]),2.0) +  pow((z_coord[o]-zv_coord[o]),2.0);
      if(del > skin2)
	verlet();
    }

 V_old = ener_part(o,x_coord[o],y_coord[o],z_coord[o],jmove);

  bias_old = ener_bias(nhattemp);
   
  xn = x_coord[o] + (range_random(0,1)-0.5)*move_max;
  yn = y_coord[o] + (range_random(0,1)-0.5)*move_max;
  zn = z_coord[o] + (range_random(0,1)-0.5)*move_max;

  if(vlist_flg == 1)
    {
      del = pow((xn-xv_coord[o]),2.0) +  pow((yn-yv_coord[o]),2.0) +  pow((zn-zv_coord[o]),2.0);
      if(del > skin2)
        verlet();
      del = pow((xn-xv_coord[o]),2.0) +  pow((yn-yv_coord[o]),2.0) +  pow((zn-zv_coord[o]),2.0);
      if(del > skin2)
	{
	  printf("\nError!! Displacement too large for Verlet\n");
	  exit(101);
	}
    }

  if(xn<0)
    xn = xn + boxL;
  if(xn>boxL)
    xn = xn - boxL;
  if(yn<0)
    yn = yn + boxL;
  if(yn>boxL)
    yn = yn - boxL;
  if(zn<0)
    zn = zn + boxL;
  if(zn>boxL)
    zn = zn - boxL;
  
  nhattemp =  update_nhat(x_coord[o],y_coord[o], z_coord[o], xn, yn, zn);
 
V_new = ener_part(o,xn,yn,zn,jmove);
  bias_new = ener_bias(nhattemp);  
 
  if(process == 2)
    {
      for(i=0;i<5;i++)
        {
          if (nhattemp == n[i])  //nhattemp is the no. of particles in the bias box once the move has been made irrespective of accepted or rejected
            count[i]++;
        }

    }

  rmove = RANDOM();
  
   EXP = exp(((V_old - V_new)+(bias_old - bias_new))/(Tstar*epsilon));

   
  if(EXP>rmove)
    {
      nacc++;
     V = V + (V_new - V_old);

      x_coord[o] = xn;
      y_coord[o] = yn;
      z_coord[o] = zn;
      
      if (flag == 1)
	nhat++;
      else if (flag == 0)
	nhat--;
     
    }
  
  nhattemp = nhat;
 
}

void verlet()
{
  double dist2;
  int j;

  nlist = vector(np+1,"mem_alloc of nlist");
  
  for(i = 1; i<=np; i++)
    {
      nlist[i] = 0;
      xv_coord[i] = x_coord[i];
      yv_coord[i] = y_coord[i];
      zv_coord[i] = z_coord[i];

    }

  for(i = 1; i<=(np-1); i++)
    {
      for(j = i+1; j<=np; j++)
	{
	  dist2 = part_dist(x_coord[i], y_coord[i], z_coord[i], j);
	  if(dist2<rv2)
	    {
	      nlist[i]++;
	      nlist[j]++;
	     
	      list[i][nlist[i]] = j;
	      list[j][nlist[j]] = i;
	      
	    }
	}


    }

}

double part_dist(double Xi, double Yi, double Zi, int j)
{
  double dx, dy, dz, r2;

  dx = Xi - x_coord[j];
  dy = Yi - y_coord[j];
  dz = Zi - z_coord[j];

  if(dx>hboxL)
    dx = boxL - dx;
  else if(dx<(-1.0*hboxL))
    dx = dx + boxL;

  if(dy>hboxL)
    dy = boxL - dy;
  else if(dy<(-1.0*hboxL))
    dy = dy + boxL;

  if(dz>hboxL)
    dz = boxL - dz;
  else if(dz<(-1.0*hboxL))
    dz = dz + boxL;

  r2 = dx*dx + dy*dy + dz*dz;

  return r2;
  
}

double ener_bias(double nh)
{
  
  ub = 0.5*k*(nh-nt)*(nh-nt);
  
  return ub;

}

double dener_bias(double nh)
{
  dub = 0.5*(nh-nt)*(nh-nt);

  return dub;
}

double update_nhat(double xo, double yo, double zo, double xc, double yc, double zc)     //o -> original (before the move), c->current(after the move)
{
  
  if((xo<bbox_disp)||(xo>(bboxL+bbox_disp)) || (yo<bbox_disp)||(yo>(bboxL+bbox_disp)) || (zo<bbox_disp)||(zo>(bboxL+bbox_disp)))
    {
      if((xc>bbox_disp)&&(xc<bbox_disp+bboxL))
	if((yc>bbox_disp)&&(yc<bbox_disp+bboxL))
	    if((zc>bbox_disp)&&(zc<bbox_disp+bboxL))
		{
		 
		  nhattemp++;
		  flag = 1;
		}	
	    
	    
    }

    
  else if(((xo>bbox_disp)&&(xo<(bbox_disp+bboxL))) && ((yo>bbox_disp)&&(yo<(bbox_disp+bboxL))) && ((zo>bbox_disp)&&(zo<(bbox_disp+bboxL))))
    {
      if((xc<bbox_disp)||(xc>(bboxL+bbox_disp)) || (yc<bbox_disp)||(yc>(bboxL+bbox_disp)) || (zc<bbox_disp)||(zc>(bboxL+bbox_disp)))
	{
	  
	  nhattemp--;
	  flag = 0;
	}
    }
  
  
  return nhattemp;
  
}
