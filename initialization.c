#include<stdio.h>
#include<math.h>
#include<assert.h>		//for assert() in dvector()
#include "common.h"

void read_input()
{
//*****************System Data*******************//
    FILE *sys_in;
    sys_in = fopen("sys_input.dat","r");

    fscanf(sys_in,"%lf",&sigma);
    printf("\nSigma (Angstrom): %lf\n",sigma);

    fscanf(sys_in,"%lf",&epsilon);
    printf("\nEpsilon (eV): %lf\n",epsilon);

    fscanf(sys_in,"%lf",&boxL_sigma);
    printf("\nBox length in sigma units: %lf\n",boxL_sigma);

    fscanf(sys_in,"%lf",&rho2star);
    printf("\nRho_env_star: %lf\n",rho2star);

    fscanf(sys_in,"%lf",&Tstar);
    printf("\nTstar: %lf\n",Tstar);

    fscanf(sys_in,"%lf",&rc_sigma);
    printf("\nRc: %lfsigma\n",rc_sigma);

    fscanf(sys_in,"%lf",&rv_sigma);
    printf("\nRv: %lfsigma\n",rv_sigma);

    fscanf(sys_in,"%lf",&k);
    printf("\nBias Strength (k): %lf\n",k);

    fscanf(sys_in,"%lf",&nt);
    printf("\nTarget number of particles: %lf\n",nt);


    fclose(sys_in);
//**********************************************//

//*****************Simulation Data**************************************************//
    FILE *sim_in;
    sim_in = fopen("sim_input.dat","r");

    fscanf(sim_in,"%lf",&totsweep);
    printf("\nTotal number of sweeps: %lf\n",totsweep);

    fscanf(sim_in,"%lf",&min_sweep);
    printf("\nSweeps after which data collection starts: %lf\n",min_sweep);

    fscanf(sim_in,"%lf",&move_max_sigma);
    printf("\nMax move size in units of sigma: %lfsigma\n",move_max_sigma);

    fscanf(sim_in,"%d",&vlist_flg);
    printf("\nVerlet List flag: %d\n",vlist_flg);


    fclose(sim_in);
//**************************************************************************************//

//*******************Defining some constants*****************************//

    move_max = move_max_sigma*sigma;
    printf("\nMax move size (in Angstrom): %lf\n",move_max);

    rc = rc_sigma*sigma;
    printf("\nCut off distance (in Angstrom): %lf\n",rc);

    rv = rv_sigma*sigma;
    printf("\nVerlet Radius (in Angstrom): %lf\n",rv);

    sig2 = sigma*sigma;
    sig3 = sig2*sigma;

    rc2 = rc*rc;

    rv2 = rv*rv;

    skin = rv - rc;
    skin2 = skin*skin;
    
    eps4 = 4*epsilon;
    eps48 = 48*epsilon;

    rho2 = rho2star/sig3;
    
    Vsum = 0.0;
    Vav = 0.0;

    beta = 1.0/(Tstar*epsilon);
}

inline double *dvector(long int n, char *s )
{
  double *dv ;
  dv = (double *)calloc( (n), sizeof(double) ) ;
   
  if (!dv)
    {
      printf(" Error in %s , Exiting ... \n",s) ;
      exit(1) ;
    }
  assert(dv) ;
  return(dv) ;
}


void sys_set()
{   
  //   srand((unsigned int)time(NULL));
  srand(238746213);  

  boxL = boxL_sigma*sigma;
    printf("\nBox length = %lfsigma = %lf\n",boxL_sigma,boxL);
   
    hboxL = boxL/2.0;
    
    bboxL_sigma = boxL_sigma/2.0;

    npenv = rho2star*(pow(boxL_sigma,3.0) - pow(bboxL_sigma,3.0));
    npenv = (int)npenv;
    npenv = (double)npenv;

    np = npenv+nt;
  
    printf("Number of particles in the box: %lf",np);
    init_config_lattice();

}

void init_config_lattice()         //Sets initial configuration of atoms in a simple cubic lattice
{
    int n, index, i, j, k;
    double del,x ,y, z;
    double test;
    
    x_coord = dvector(np+1,"mem_alloc of x_coord");       
    y_coord = dvector(np+1,"mem_alloc of y_coord");
    z_coord = dvector(np+1,"mem_alloc of z_coord");

    /*np+1 and not np is used above because x_coord[1] is the first element in the definition below, thus memory
      for the 0th element is not accounted for if np is used resulting in a glibc detected error*/
   
    //For the Verlet List
    if(vlist_flg == 1)
      {
	xv_coord = dvector(np+1,"mem_alloc of xv_coord");
	yv_coord = dvector(np+1,"mem_alloc of yv_coord");
	zv_coord = dvector(np+1,"mem_alloc of zv_coord");

	list  = (int **)malloc(sizeof(int *)*(np+1));
	for(i=1; i <=np; i++)
	  {
	    list[i] = (int *)malloc(sizeof(int)*(np+1));
	  } 

      }
    
    n = (int)(pow(np,(1.0/3.0))+1);

    if(n==0)
        n=1;
    del = boxL/(double)n;
    index = 0;
    x = -1.0*del;
    for(i=1;i<=n;i++)
    {
        x = x+del;
        y = -1*del;

        for(j=1;j<=n;j++)
        {
            y = y+del;
            z = -1*del;

            for(k=1;k<=n;k++)
            {
                z = z+del;
                if(index<np)
                {
                    index++;
                    x_coord[index] = x;
                    y_coord[index] = y;
                    z_coord[index] = z;
                }
            }
        }
    }

    if(vlist_flg == 1)
      verlet();
    
}

void set_box()
{
  nhat = 0.0;
  bboxL = hboxL;
  printf("\nBiased Box length = %lfsigma = %lf\n",bboxL_sigma,bboxL);

  bbox_disp = hboxL/2.0;
  
  for(i=1;i<=np;i++)
    {
      if((x_coord[i]>bbox_disp)&&(x_coord[i]<bbox_disp+bboxL))
	if((y_coord[i]>bbox_disp)&&(y_coord[i]<bbox_disp+bboxL))
	  if((z_coord[i]>bbox_disp)&&(z_coord[i]<bbox_disp+bboxL))
	    nhat++;
    }
  nhattemp = nhat;
  printf("\nNumber of particles in the inner box: %lf\n",nhat);

}


