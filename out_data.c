#include<stdio.h>
#include "common.h"
#include<math.h>

void coord_write(int sweep)
{   int np_write;
    np_write = (int)np;
    char filename[] = "coord_";
    char extension[] = ".xyz";
    char fsweep[1000000];
    sprintf(fsweep,"%d",sweep);

    strcat(filename, fsweep);
    strcat(filename, extension);

    FILE *writecoord;
    writecoord = fopen(filename,"w");
    fprintf(writecoord,"%i\n\n",np_write);
    for(i=1;i<=np;i++)
    fprintf(writecoord,"%lf %lf %lf\n",x_coord[i],y_coord[i],z_coord[i]);
    fclose(writecoord);
}
