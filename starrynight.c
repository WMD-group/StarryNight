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

#define X 50  // Malloc is for losers.
#define Y 50

struct dipole
{
    float angle;
} lattice[X][Y];

float beta=2.0;  // beta=1/T  T=temperature of the lattice, in units of k_B

float Efield=0.0;
float Eangle=0.0;

float K=50.0; //elastic coupling constant for dipole moving within cage

unsigned long ACCEPT=0; //counters for MC moves
unsigned long REJECT=0;

// Prototypes...
int rand_int(int SPAN);
double site_energy(int x, int y, double newangle, double oldangle);
void MC_move();
double lattice_energy();
void outputlattice_pnm(char * filename);
void outputlattice_ppm_hsv(char * filename);
void outputlattice_svg(char * filename);

int main(void)
{
    int i,j,k; //for loop iterators

    char name[50]; //for output filenames

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

// If we're going to do some actual science, we better have a logfile...
    FILE *log;
    log=fopen("starry.log","w");

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

    for (i=0;i<400;i++)
    {
        // Log some data...
        fprintf(log,"%d %f\n",ACCEPT+REJECT,lattice_energy());

        sprintf(name,"MC-PNG_step_%.4d.png",i);
        outputlattice_ppm_hsv(name);

        sprintf(name,"MC-SVG_step_%.4d.svg",i);
        outputlattice_svg(name);

        fprintf(stderr,".");

        if (i==50)  { Efield=1.0; Eangle=M_PI/2;}
        if (i%100) { Efield=-Efield;}

        for (k=0;k<X*Y;k++)
            MC_move();
    }

    fprintf(stderr,"\n");

    outputlattice_ppm_hsv("final.png");
    outputlattice_svg("final.svg");

    fprintf(stderr,"ACCEPT: %lu REJECT: %lu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

int rand_int(int SPAN)
{
    return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

double site_energy(int x, int y, double newangle, double oldangle)
{
    int dx,dy;
    float d;
    double dE=0.0;
    double testangle;

    // Sum over near neighbours for dipole-dipole interaction
    for (dx=-2;dx<=2;dx++)
        for (dy=-2;dy<=2;dy++)
        {
            if (dx==0 && dy==0)
                continue; //no infinities / self interactions please!

            d=sqrt((float) dx*dx + dy*dy); //that old chestnut

            if (d>2.0) continue; // Cutoff in d

            testangle=lattice[(X+x+dx)%X][(Y+y+dy)%Y].angle;

            //it goes without saying that the following line is the single
            //most important in the program... Energy calculation!
            dE+=  + cos(newangle-testangle)/(d*d*d)
                  - cos(oldangle-testangle)/(d*d*d);
        }

    // Interaction of dipole with (unshielded) E-field
    dE+= + Efield*cos(newangle-Eangle)
         - Efield*cos(oldangle-Eangle);

    // interaction with strain of cage
    dE += + K*sin(2*newangle)*sin(2*newangle)
          - K*sin(2*oldangle)*sin(2*oldangle);


    return(dE); 
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

    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well
    newangle=2*M_PI*genrand_real2();

    // comparison point for the dE - the present configuration
    oldangle=lattice[x][y].angle;

    //calc site energy
    //double site_energy(int x, int y, double newangle, double oldangle);
    dE=site_energy(x,y,newangle,oldangle);

    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
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

double lattice_energy()
{
    int x,y,dx,dy;
    double E=0.0,d,oldangle,testangle;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
        {

// NB: just copied + pasted this code :| - should probably generalise to
// a function, otherwise variations in cutoff / Hamiltonian will have to be in
// two places, ugh.
        oldangle=lattice[x][y].angle;
        
            // Sum over near neighbours for dipole-dipole interaction
            for (dx=-2;dx<=2;dx++)
                for (dy=-2;dy<=2;dy++)
                {
                    if (dx==0 && dy==0)
                        continue; //no infinities / self interactions please!

                    d=sqrt((float) dx*dx + dy*dy); //that old chestnut

                    if (d>2.0) continue; // Cutoff in d

                    testangle=lattice[(X+x+dx)%X][(Y+y+dy)%Y].angle;

                    //it goes without saying that the following line is the single
                    //most important in the program... Energy calculation!
                    E+=   cos(oldangle-testangle)/(d*d*d);
                }

            // Interaction of dipole with (unshielded) E-field
            E+=   Efield*cos(oldangle-Eangle);
        }

//    fprintf(stderr,"Energy of lattice: %f\n",E);

    return(E);
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
    s=0.4; v=0.8;

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
    fprintf(fo," <marker id=\"triangle\" viewBox=\"0 0 10 10\" refX=\"10\" refY=\"5\" markerUnits=\"strokeWidth\" markerWidth=\"4\" markerHeight=\"3\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" /></marker>\n");

    //No markers...  marker-end=\"url(#triangle)\"

     for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            fprintf(fo," <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:rgb(0,0,0);stroke-width:0.1\" marker-end=\"url(#triangle)\" />\n",
                    i+0.5 - 0.5*sin(lattice[k][i].angle), 
                    k+0.5 - 0.5*cos(lattice[k][i].angle),
                    i+0.5 + 0.5*sin(lattice[k][i].angle),
                    k+0.5 + 0.5*cos(lattice[k][i].angle)
                   );
    
    fprintf(fo,"</svg>\n");

    fclose(fo);
}


