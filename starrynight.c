/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mt19937ar-cok.c"

#define X 500  // Malloc is for losers.
#define Y 500

struct dipole
{
    float angle;
} lattice[X][Y];

float beta=1.0;  // beta=1/T  T=temperature of the lattice, in units of k_B
float Efield=1.0;

unsigned long ACCEPT=0; //counters for MC moves
unsigned long REJECT=0;

// Prototypes...
int rand_int(int SPAN);
void MC_move();
void outputlattice_pnm(char * filename);
void outputlattice_ppm_hsv(char * filename);
void outputlattice_svg(char * filename);

int main(void)
{
    int i,j,k; //for loop iterators

    char name[50]; //for output filenames

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

    //Fire up the twister!
    //init_genrand(0);  // reproducible
    init_genrand(time(NULL)); // seeded with current time

    //Random initial lattice
     for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            lattice[i][k].angle=2*M_PI*genrand_real2(); // randomised initial orientation of dipoles
//           lattice[i][k].angle=2*M_PI*(i*X+k)/((float)X*Y); // continous set
//           of dipole orientations to test colour output (should appear as
//           spectrum)

    //Print lattice
/*    for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            printf(" %f",lattice[i][k].angle);
*/
    outputlattice_ppm_hsv("initial.pnm");

    for (i=0;i<1000;i++)
    {
         sprintf(name,"MC_step_%.4d.pnm",i);
        outputlattice_ppm_hsv(name);
        
        fprintf(stderr,".");

        if (i%200==0) Efield=0.0-Efield;

        for (k=0;k<1e5;k++)
            MC_move();
    }

    fprintf(stderr,"\n");

    outputlattice_ppm_hsv("final.pnm");
    outputlattice_svg("final.svg");

    fprintf(stderr,"ACCEPT: %lu REJECT: %lu ratio: %f",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));

    return 0;
}

int rand_int(int SPAN)
{
    return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

void MC_move()
{
    int x, y;
    int dx, dy;
    float newangle,oldangle,testangle,d;
    float dE=0.0;

    // Choose random dipole / lattice location

    x=rand_int(X);
    y=rand_int(Y);

    // random new orientation
    newangle=2*M_PI*genrand_real2();
    oldangle=lattice[x][y].angle;


    // Sum over near neighbours for dipole-dipole interaction
    for (dx=-2;dx<=2;dx++)
        for (dy=-2;dy<=2;dy++)
        {
            if (dx==0 && dy==0)
                break; //no infinities / self interactions please!

            d=sqrt((float) dx*dx + dy*dy); //that old chestnut

            testangle=lattice[(x+dx)%X][(y+dy)%Y].angle;

            //it goes without saying that the following line is the single
            //most important in the program... Energy calculation!
            dE+=  + cos(newangle-testangle)/(d*d*d)
                  - cos(oldangle-testangle)/(d*d*d);
        }

    // Interaction of dipole with (unshielded) E-field
    dE+= + Efield*cos(newangle)
         - Efield*cos(oldangle);

    if (dE > 0.0 || exp(dE * beta) > genrand_real2() )
    {
        lattice[x][y].angle=newangle;
        ACCEPT++;
    }
    else
        REJECT++;

    // DEBUGGING / OUTPUT PRINTS
/*
    if (lattice[x][y].angle==newangle)
        fprintf(stderr,"A "); //i.e. accepted move
    else
        fprintf(stderr,"R ");

    fprintf(stderr,"MC: %d X %d Y oldangle %f newangle %f dE: %f\n",x,y,oldangle,newangle,dE);
*/
}

void outputlattice_png(char * filename)
{
    int i,k;
    FILE *fo;
    fo=fopen(filename,"w");

    fprintf (fo,"P2\n%d %d\n%d\n", X, Y, SHRT_MAX);

    for (i=0;i<X;i++)
    {
        for (k=0;k<Y;k++)
            fprintf(fo,"%d ",(int)(SHRT_MAX*lattice[i][k].angle/(2*M_PI)));
        fprintf(fo,"\n");
    }

}

void outputlattice_ppm_hsv(char * filename)
{
    int i,k;
    float angle;

    float r,g,b; // RGB
    float h,s,v; // HSV
    float p,t,q,f; // intemediates for HSV->RGB conversion
    int hp;
    
    FILE *fo;
    fo=fopen(filename,"w");

//Set Saturation + Value, vary hue
    s=0.8; v=0.8;

    fprintf (fo,"P6\n%d %d\n255\n", X, Y);

    for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
        {
            h=lattice[i][k].angle;

            // http://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
            hp=(int)floor(h/(M_PI/3.0))%6; //radians, woo
            f=h/(M_PI/3.0)-floor(h/(M_PI/3.0));
            p=v*(1.0-s);
            q=v*(1.0-f*s);
            t=v*(1.0-(1.0-f)*s);
            
            switch (hp){
                case 0: r=v; g=t; b=p; break;
                case 1: r=q; g=v; b=p; break;
                case 2: r=p; g=v; b=t; break;
                case 3: r=p; g=q; b=v; break;
                case 4: r=t; g=p; b=v; break;
                case 5: r=v; g=p; b=q; break;
            }

            fprintf(fo,"%c%c%c",(char)(254.0*r),(char)(254.0*g),(char)(254.0*b));
        }
    fclose(fo); //don't forget :^)
}

void outputlattice_svg(char * filename)
{
    int i,k;

    FILE *fo;
    fo=fopen(filename,"w");

    fprintf(fo,"<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");

    //our arrow marker...
    fprintf(fo," <marker id=\"triangle\" viewBox=\"0 0 10 10\" refX=\"0\" refY=\"5\" markerUnits=\"strokeWidth\" markerWidth=\"4\" markerHeight=\"3\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" /></marker>\n");

     for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            fprintf(fo," <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:rgb(0,0,0);stroke-width:2.0\" marker-end=\"url(#triangle)\"/>\n",
                    (double) 10.0*i, (double) 10.0*k,
                    10.0*(i+0.9*sin(lattice[i][k].angle)),
                    10.0*(k+0.9*cos(lattice[i][k].angle))
                   );
    
    fprintf(fo,"</svg>\n");

    fclose(fo);
}


