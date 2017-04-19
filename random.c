#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include "common.h" 

double range_random(double min, double max)
{
  double rn;
  rn = min + rand()/(RAND_MAX/(max-min));
  return rn;

}

double RANDOM()
{
  double rn;

  rn = (double) rand()/(RAND_MAX+1.0);

  return rn;
}
