#include<assert.h>


//Functions
void read_input();                  //Reads input from sys_input.dat(system data) and sim_input.dat(simulation data)
void sys_set();                      //Sets up the system configuration
void init_config_lattice();
void coord_write(int);
double cor_U(double, double);
double cor_P(double, double);
double tot_energy();
double ener_part(int, double, double, double, int);
double calc_pot(double);
void mmc();
double adjust(int, int, double);
void mcmove();
double range_random(double, double);
double RANDOM();
void verlet();
double part_dist(double, double, double, int);
void set_box();
double ener_bias(double);
double dener_bias(double);
double update_nhat(double, double, double, double, double, double);
void chem_potential();
double bias_freenerg();
double bias_onoff();
double integrand(double);

//void file_init();

//*************************************Variables****************************************************************//

//Defined in main()
int i, process, nsweep,run_num,prog_flag;  
double V, V_corr, P_corr;                                  

//Defined in read_input()
double sigma, epsilon, np, npenv, rho2star, Tstar, rc_sigma, rc, rv_sigma, rv;    //Section 1
int vlist_flg;
double totsweep, min_sweep, write_sweep, move_max_sigma, boxL_sigma, bboxL_sigma, boxL, hboxL, *x_coord, *y_coord, *z_coord, *xv_coord, *yv_coord, *zv_coord; //Section 2
double move_max, sig2, sig3, rho2, eps4, eps48, rc2, rv2, skin, skin2, Vsum, Vav, beta;     //Section 3
int **list;

//defined in set_box()
double bboxL, bbox_disp,nhat, nhattemp;

//defined in ener_bias()
double nt, k, ub;

// defined in dener_bias()
double dub;

//defined in mmc()
int sweep, move;
double nhatsum, nhatav, nacc, nattempt,bias_sum, bias_av,Vsum,Vav;

//defined in adjust()
double perc;
int Att, Acc;
//defined in mcmove()
double r, bias_old, bias_new;
int o, flag;

//defined in verlet()
long int *nlist;

//defined in chem_potential()
double Ab_nhatavg, ub_nhatav, A_nhatavg, nhatav_main, k_main;

//defined in bias_freenerg()
double count[5], n[5], p_nhat[5], Ab_nhat[5];

//defined in bias_onoff()
double therm_integral;

//extern FILE *nhat_av;

//**************************************************************************************************************//
