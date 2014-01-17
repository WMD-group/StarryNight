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

#define X 100  // Malloc is for losers.
#define Y 100

struct dipole
{
    float angle;
} lattice[X][Y];

float beta=1.0;  // beta=1/T  T=temperature of the lattice, in units of k_B

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
//            lattice[i][k].angle=M_PI/2;
//            lattice[i][k].angle=2*M_PI*(i*X+k)/((float)X*Y); // continous set
//           of dipole orientations to test colour output (should appear as
//           spectrum)

    //Print lattice
/*    for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            printf(" %f",lattice[i][k].angle);
*/
    outputlattice_ppm_hsv("initial.png");

    for (i=0;i<200;i++)
    {
        sprintf(name,"MC-PNG_step_%.4d.png",i);
        outputlattice_ppm_hsv(name);

        sprintf(name,"MC-SVG_step_%.4d.svg",i);
        outputlattice_svg(name);

        fprintf(stderr,".");

        for (k=0;k<10*X*Y;k++)
            MC_move();
    }

    fprintf(stderr,"\n");

    outputlattice_ppm_hsv("final.png");
    outputlattice_svg("final.svg");

    fprintf(stderr,"ACCEPT: %lu REJECT: %lu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\b");

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

    for (dx=-2;dx<=2;dx++)
        for (dy=-2;dy<=2;dy++)
        {
            if (dx==0 && dy==0)
                continue; //no infinities / self interactions please!

            d=sqrt((float) dx*dx + dy*dy); //that old chestnut

            if (d>2.0) continue; // Cutoff in d

            testangle=lattice[(x+dx)%X][(y+dy)%Y].angle;

            dE+=  + cos(newangle-testangle)/(d*d*d)
                  - cos(oldangle-testangle)/(d*d*d);
        }

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

    for (i=0;i<X;i++) //force same ordering as SVG...
        for (k=0;k<Y;k++)
        {
            h=fmod(lattice[i][k].angle,M_PI*2.0);

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

//            fprintf(stderr,"h: %f r: %f g: %f b: %f\n",h,r,g,b);

            fprintf(fo,"%c%c%c",(char)(254.0*r),(char)(254.0*g),(char)(254.0*b));
        }
    fclose(fo); //don't forget :^)
}

void outputlattice_svg(char * filename)
{
    int i,k;

    FILE *fo;
    fo=fopen(filename,"w");

    fprintf(fo,"<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" height=\"%d\" width=\"%d\">\n",X,Y);

    //our arrow marker...
    fprintf(fo," <marker id=\"triangle\" viewBox=\"0 0 10 10\" refX=\"0\" refY=\"5\" markerUnits=\"strokeWidth\" markerWidth=\"4\" markerHeight=\"3\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" /></marker>\n");

    //No markers...  marker-end=\"url(#triangle)\"

     for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            fprintf(fo," <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:rgb(0,0,0);stroke-width:0.1\" marker-end=\"url(#triangle)\"  />\n",
                    i+0.5 - 0.4*sin(lattice[k][i].angle), 
                    k+0.5 - 0.4*cos(lattice[k][i].angle),
                    i+0.5 + 0.4*sin(lattice[k][i].angle),
                    k+0.5 + 0.9*cos(lattice[k][i].angle)
                   );
    
    fprintf(fo,"</svg>\n");

    fclose(fo);
}


