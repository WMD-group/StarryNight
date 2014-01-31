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
#include <libconfig.h>

#include "mt19937ar-cok.c"

#define X 50  // Malloc is for losers.
#define Y 50

struct dipole
{
    float x,y,z;
    float angle;
} lattice[X][Y];

// SIMULATION PARAMETERS

double beta=1.0;  // beta=1/T  T=temperature of the lattice, in units of k_B

struct dipole Efield; //now a vector, still k_B.T units per lattice unit
//double Efield=0.01; // units k_B.T per lattice unit
double Eangle=0.0;

double K=1.0; //elastic coupling constant for dipole moving within cage

double Dipole=1.0; //units of k_B.T for spacing = 1 lattice unit

int DIM=2; 

//END OF SIMULATION PARAMETERS
// Except for the ones hardcoded into the algorithm :^)

unsigned long ACCEPT=0; //counters for MC moves
unsigned long REJECT=0;


// Prototypes...
static int rand_int(int SPAN);
static double site_energy(int x, int y, struct dipole *newdipole, struct dipole *olddipole);
static void MC_move();
void initialise_lattice();
static double lattice_energy_log(FILE *log);
void outputlattice_pnm(char * filename);
void outputlattice_ppm_hsv(char * filename);
void outputlattice_svg(char * filename);

int main(void)
{
    int i,j,k; //for loop iterators
    int MCMegaSteps=400;
    double MCMegaMultiplier=1.0;
    config_t cfg, *cf; //libconfig config structure
    double tmp;

    char name[50]; //for output filenames

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

//Load and parse config file
    cf = &cfg;
    config_init(cf);

    if (!config_read_file(cf,"starrynight.cfg")) 
    {
        fprintf(stderr, "%s:%d - %s\n",
                config_error_file(cf),
                config_error_line(cf),
                config_error_text(cf));
        config_destroy(cf);
        return(EXIT_FAILURE);
    }

    config_lookup_float(cf,"beta",&beta);

    config_lookup_float(cf,"Efield.x",&tmp);  Efield.x=(float)tmp;
    config_lookup_float(cf,"Efield.y",&tmp);  Efield.y=(float)tmp;
    config_lookup_float(cf,"Efield.z",&tmp);  Efield.z=(float)tmp;

//    config_lookup_float(cf,"Eangle",&Eangle);

    config_lookup_float(cf,"K",&K);
    config_lookup_float(cf,"Dipole",&Dipole);

    config_lookup_int(cf,"MCMegaSteps",&MCMegaSteps);
    config_lookup_float(cf,"MCMegaMultiplier",&MCMegaMultiplier);

    fprintf(stderr,"Config loaded. ");

// If we're going to do some actual science, we better have a logfile...
    FILE *log;
    log=fopen("starrynight.log","w");
    fprintf(stderr,"Log file opened. ");

    //Fire up the twister!
    //init_genrand(0);  // reproducible
    init_genrand(time(NULL)); // seeded with current time
    fprintf(stderr,"Twister initialised. ");

    initialise_lattice();
    fprintf(stderr,"Lattice initialised.");

    outputlattice_ppm_hsv("initial.png");

    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");

    fprintf(stderr,"'.' is %d MC moves attempted.\n",(int)(X*Y*MCMegaMultiplier));

    fprintf(log,"# ACCEPT+REJECT, Efield, Eangle, E_dipole, E_strain, E_field, (E_dipole+E_strain+E_field)\n");

    for (i=0;i<MCMegaSteps;i++)
    {
        // Log some data... Nb: Slow as does a NxN summation of lattice energy
        // contributions!
        lattice_energy_log(log);
        //fprintf(log,"%lu %f %f %f\n",ACCEPT+REJECT,lattice_energy(),Efield,Eangle); //FIXME: lattice_energy all broken, data worthless presently.
        // TODO: some kind of dipole distribution? Would I have to bin it
        // myself? (boring.)
        // TODO: Split Energy into different contributions... would be nice to
        // see polarisation delta.E spike when the field flips

        // Log some pretty pictures...
        sprintf(name,"MC-PNG_step_%.4d.png",i);
        outputlattice_ppm_hsv(name);

        sprintf(name,"MC-SVG_step_%.4d.svg",i);
        outputlattice_svg(name);

        // Update the (interactive) user what we're up to
        fprintf(stderr,".");

        // Manipulate the run conditions depending on simulation time
        if (i==100) { DIM=3;}  // ESCAPE FROM FLATLAND
        if (i==200) { Efield.z=1.0;}      // relax back to nothing
        if (i==300) {Efield.z=0.0; Efield.x=1.0;}

        // Do some MC moves!
        for (k=0;k<X*Y*MCMegaMultiplier;k++) //let's hope the compiler inlines this to avoid stack abuse. Alternatively move core loop to MC_move fn?
            MC_move();
    }

    // OK; we're finished...

    fprintf(stderr,"\n");

    // Final data output / summaries.
    outputlattice_ppm_hsv("final.png");
    outputlattice_svg("final.svg");

    fprintf(stderr,"ACCEPT: %lu REJECT: %lu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

static void random_sphere_point(struct dipole *p)
{
    int i;
    // Marsaglia 1972 
    float x1,x2;
    do {
        x1=2.0*genrand_real1() - 1.0;
        x2=2.0*genrand_real1() - 1.0;
    } while (x1*x1 + x2*x2 > 1.0);

    if (DIM<3){
        // Circle picking, after Cook 1957
        // http://mathworld.wolfram.com/CirclePointPicking.html
        p->x = (x1*x1 - x2*x2)  / (x1*x1 + x2*x2);
        p->y =      2*x1*x2     / (x1*x1 + x2*x2);
        p->z = 0.0;
    }
    else
    {
        // Sphere picking
        p->x = 2*x1*sqrt(1-x1*x1-x2*x2);
        p->y = 2*x2*sqrt(1-x1*x1-x2*x2);
        p->z = 1.0 - 2.0* (x1*x1+x2*x2);
    }
}

static float dot(struct dipole *a, struct dipole *b)
{
    int D;
    float sum=0.0;

    sum+=a->x*b->x;
    sum+=a->y*b->y;
    sum+=a->z*b->z;

    return(sum);
}

void initialise_lattice()
{
    int i,k;
    //Random initial lattice
     for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
//            random_sphere_point(& lattice[i][k]);
//            lattice[i][k].angle=2*M_PI*genrand_real2(); // randomised initial orientation of dipoles
//            lattice[i][k].angle=M_PI/2;
     {
         lattice[i][k].angle=2*M_PI*(i*X+k)/((float)X*Y); // continous set
//           of dipole orientations to test colour output (should appear as
//           spectrum)
         lattice[i][k].x = sin(lattice[i][k].angle);
         lattice[i][k].y=cos(lattice[i][k].angle);
         lattice[i][k].z=0.0;
     }

    //Print lattice
    for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            printf("\n %f %f %f %f",lattice[i][k].x,lattice[i][k].y,lattice[i][k].z,
                    dot(&lattice[i][k],&lattice[i][k]));
}

static int rand_int(int SPAN) // TODO: profile this to make sure it runs at an OK speed.
{
    return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

static double site_energy(int x, int y, struct dipole *newdipole, struct dipole *olddipole)
{
    int dx,dy;
    float d;
    double dE=0.0;
    struct dipole *testdipole, n;

    // Sum over near neighbours for dipole-dipole interaction
    for (dx=-2;dx<=2;dx++)
        for (dy=-2;dy<=2;dy++)
        {
            if (dx==0 && dy==0)
                continue; //no infinities / self interactions please!

            d=sqrt((float) dx*dx + dy*dy); //that old chestnut

            if (d>2.0) continue; // Cutoff in d

            testdipole=& lattice[(X+x+dx)%X][(Y+y+dy)%Y];
//            testangle=lattice[(X+x+dx)%X][(Y+y+dy)%Y].angle;

            //it goes without saying that the following line is the single
            //most important in the program... Energy calculation!

            n.x=(float)dx/d; n.y=(float)dy/d; //normalised diff. vector

//            n=atan2((float)dy,(float)dx); //angle of normal vector between test points
            // Anti-ferroelectric (dipole like)
            //  - this now contains a lot of trig to do the dot products. Maybe
            //  faster to generate the vectors and do it component wise?
//            dE+=  + Dipole * ( cos(newangle-testangle) - 3.* cos(n-newangle) * cos(n-testangle) ) /(d*d*d) 
//                  - Dipole * ( cos(oldangle-testangle) - 3.* cos(n-oldangle) * cos(n-testangle) ) /(d*d*d) ;
            // Ferroelectric / Potts model
//            dE+=  - Dipole * cos(newangle-testangle)/(d*d*d)
//                  + Dipole * cos(oldangle-testangle)/(d*d*d);

            //True dipole like
            dE+= - Dipole * ( dot(newdipole,testdipole) - 3*dot(&n,newdipole)*dot(&n,testdipole) ) / (d*d*d)
                 + Dipole * ( dot(olddipole,testdipole) - 3*dot(&n,olddipole)*dot(&n,testdipole) ) / (d*d*d); 

            // Ferroelectric / Potts model - vector form
//            dE+= - Dipole * dot(newdipole,testdipole) / (d*d*d)
//                + Dipole * dot(olddipole,testdipole) / (d*d*d);
        }

    // Interaction of dipole with (unshielded) E-field
//    dE+= + Efield*cos(newangle-Eangle)
//         - Efield*cos(oldangle-Eangle);
    dE+= + dot(newdipole, & Efield)
         - dot(olddipole, & Efield);
    //fprintf(stderr,"%f\n",dot(newdipole, & Efield));

    // interaction with strain of cage modelled as cos^2 function (low energy
    // is diagonal with MA ion along hypotenuse)
//    dE += + K*cos(2*newangle)*cos(2*newangle)
//          - K*cos(2*oldangle)*cos(2*oldangle);

    return(dE); 
}

static void MC_move()
{
    int x, y;
    int dx, dy;
    float d;
    float dE=0.0;
    struct dipole newdipole, *olddipole;

    // Choose random dipole / lattice location

    x=rand_int(X);
    y=rand_int(Y);

    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well
//    newangle=2*M_PI*genrand_real2();
    random_sphere_point(& newdipole);    

    // comparison point for the dE - the present configuration
//    oldangle=lattice[x][y].angle;

    olddipole=& lattice[x][y];

    //calc site energy
    //double site_energy(int x, int y, double newangle, double oldangle);
    dE=site_energy(x,y,& newdipole,olddipole);

    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x][y].x=newdipole.x;
        lattice[x][y].y=newdipole.y;
        lattice[x][y].z=newdipole.z;

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

static double lattice_energy_log(FILE *log)
{
    int x,y,dx,dy;
    double E_dipole=0.0,E_strain=0.0,E_field=0.0;
    double d,oldangle,testangle,n;

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
            n=atan2((float)dy,(float)dx); //angle of normal vector between test points
            // Anti-ferroelectric (dipole like)
            //  - this now contains a lot of trig to do the dot products. Maybe
            //  faster to generate the vectors and do it component wise?
            E_dipole+=   Dipole * ( cos(oldangle-testangle) - 3.* cos(n-oldangle) * cos(n-testangle) ) /(d*d*d) ;
     
            // Ferroelectric / Potts model
//            dE+=  - Dipole * cos(newangle-testangle)/(d*d*d);
 
                }

            // Interaction of dipole with (unshielded) E-field
//            E_field+=  Efield*cos(oldangle-Eangle);

            //Interaction with cage
            E_strain+=  K*sin(2*oldangle)*sin(2*oldangle);
        }


//    fprintf(stderr,"Energy of lattice: %f\n",E);

    fprintf(log,"%lu %f %f %f %f %f %f\n",ACCEPT+REJECT,Efield.x,Eangle,E_dipole,E_strain,E_field,E_dipole+E_strain+E_field);

    return(E_dipole+E_strain+E_field); //FIXME: is this still useful?
}

// TODO: move these output routines to a separate file...

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
    s=0.6; v=0.8;

    fprintf (fo,"P6\n%d %d\n255\n", X, Y);

    for (i=0;i<X;i++) //force same ordering as SVG...
        for (k=0;k<Y;k++)
        {
            h=M_PI+atan2(lattice[i][k].y,lattice[i][k].x); //Nb: assumes 0->2PI interval!
            //h=fmod(lattice[i][k].angle,M_PI*2.0); //old angle code
            v=0.5+0.4*lattice[i][k].z; //darken towards the south (-z) pole
            s=0.6-0.6*fabs(lattice[i][k].z); //desaturate towards the poles

            // http://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
            hp=(int)floor(h/(M_PI/3.0)); //radians, woo
            f=h/(M_PI/3.0)-(float)hp;
            
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
    fprintf(fo," <marker id=\"triangle\" viewBox=\"0 0 10 10\" refX=\"7\" refY=\"5\" markerUnits=\"strokeWidth\" markerWidth=\"4\" markerHeight=\"3\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" /></marker>\n");

    //No markers...  marker-end=\"url(#triangle)\"

     for (i=0;i<X;i++)
        for (k=0;k<Y;k++)
            fprintf(fo," <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:rgb(%d,%d,%d);stroke-width:0.17\" marker-end=\"url(#triangle)\" />\n",
                    i+0.5 - 0.4*lattice[k][i].x, 
                    k+0.5 - 0.4*lattice[k][i].y,
                    i+0.5 + 0.4*lattice[k][i].x,
                    k+0.5 + 0.4*lattice[k][i].y,
                    (int)((-lattice[k][i].z+1.0)*127.0),
                    (int)((-lattice[k][i].z+1.0)*127.0),
                    (int)((-lattice[k][i].z+1.0)*127.0)
                   );
     // invert z-axis, and scale to greyscale. Therefore alternates with
     // pointing up and down with background colour
    
    fprintf(fo,"</svg>\n");

    fclose(fo);
}


