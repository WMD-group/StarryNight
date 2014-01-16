/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#include <math.h>
#include <stdio.h>
#include "mt19937ar-cok.c"

#define LATTICE_SIZE 500  // Malloc is for losers.

struct dipole
{
    float angle;
} lattice[LATTICE_SIZE][LATTICE_SIZE];

int main(void)
{
    int i,j,k; //for loop iterators

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

    //Fire up the twister!
    init_genrand(0);

    //Random initial lattice
     for (i=0;i<LATTICE_SIZE;i++)
        for (k=0;k<LATTICE_SIZE;k++)
            lattice[i][k].angle=2*M_PI*genrand_real2();

    //Print lattice
    for (i=0;i<LATTICE_SIZE;i++)
        for (k=0;k<LATTICE_SIZE;k++)
            printf(" %f",lattice[i][k].angle);

    return 0;
}
