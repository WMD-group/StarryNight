/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#include <stdbool.h>

#define X 32 // Malloc is for losers.
#define Y 32 
#define Z 1 

int DIM=2; // if DIM==2, the dipoles are constrained to the XY plane
// i.e. a model for dipoles in the Tetragonal phase of MAPI near the
// Orthorhombic-Tetagonal phase transition, where they are constrained to be
// within the tetragonal plane
// if DIM==3, full freedom over the sphere is considered
// i.e. a model for the Tetra phase near the end of its second order transition
// to a fully cubic phase

int T; //global variable so accessible to analysis routines

unsigned long ACCEPT=0; //counters for MC moves
unsigned long REJECT=0;

// CUSTOM STRUCTURES
struct dipole
{
    float x,y,z;
    float length; //length of dipole, to allow for solid state mixture (MA, FA, Ammonia, etc.)
} lattice[X][Y][Z];
// Nb: would be trivial to make this a runtime option with memory allocation;
// BUT - this way the X,Y,Z are constants and can be inspected

struct mixture
{
    float length;
    float prevalence;
} dipoles[10];
int dipolecount=0;

// Prototypes...
static float dot(struct dipole *a, struct dipole *b);
static void random_sphere_point(struct dipole *p);
static void random_X_point(struct dipole *p);
void initialise_lattice();
void initialise_lattice_wall();
void initialise_lattice_slip();
void initialise_lattice_spectrum();
void initialise_lattice_buckled();
void initialise_lattice_slab_delete();

// SIMULATION PARAMETERS
// NB: These are defaults - most are now read from config file

double beta=1.0;  // beta=1/T  T=temperature of the lattice, in units of k_B

struct dipole Efield; //now a vector, still k_B.T units per lattice unit

double K=1.0; //elastic coupling constant for dipole moving within cage
//double Dipole=1.0; //units of k_B.T for spacing = 1 lattice unit
double CageStrain=1.0; // as above

//double dipole_fraction=1.0; //fraction of sites to be occupied by dipoles
int DipoleCutOff=3;

// These variables control the number of loops
int MCMegaSteps=400;
int MCEqmSteps=10;
double MCMegaMultiplier=1.0;
long long int MCMinorSteps=0;
char const *LOGFILE = NULL; //for output filenames

// Simulation display / calculation flags
// False = 0 ; True = 1
int DisplayDumbTerminal=true;
int CalculateRecombination=true;
int CalculateRadialOrderParameter=false;

int ConstrainToX=false;

int CalculatePotential=false;
int CalculateEfield=false;

int SaveDipolesXYZ=false;
int SaveDipolesPNG=false;
int SaveDipolesSVG=false;

//END OF SIMULATION PARAMETERS
// {{ Except for the ones hardcoded into the algorithm :^) }}

void load_config()
{
    int i,j,k, x,y; //for loop iterators

    config_t cfg, *cf; //libconfig config structure
    const config_setting_t *setting;
    double tmp;

    char name[100];
    // Yes, I know, 50 chars are enough for any segfault ^_^

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
        exit(EXIT_FAILURE);
    }

    config_lookup_string(cf,"LOGFILE",&LOGFILE); //library does its own dynamic allocation

    config_lookup_int(cf,"T",&T);

    config_lookup_float(cf,"Efield.x",&tmp);  Efield.x=(float)tmp;
    config_lookup_float(cf,"Efield.y",&tmp);  Efield.y=(float)tmp;
    config_lookup_float(cf,"Efield.z",&tmp);  Efield.z=(float)tmp;

    fprintf(stderr,"Efield: x %f y %f z %f\n",Efield.x,Efield.y,Efield.z);

    //    config_lookup_float(cf,"Eangle",&Eangle);

    config_lookup_float(cf,"K",&K);
//    config_lookup_float(cf,"Dipole",&Dipole); DEPRECATED - see below for list
    config_lookup_float(cf,"CageStrain",&CageStrain);

    
    // read in list of dipoles + prevalence for solid mixture
    setting = config_lookup(cf, "Dipoles");
    dipolecount   = config_setting_length(setting);
    for (i=0;i<dipolecount;i++)
    dipoles[i].length=config_setting_get_float_elem(setting,i);
    setting = config_lookup(cf, "Prevalence");
    dipolecount   = config_setting_length(setting);
    for (i=0;i<dipolecount;i++)
    dipoles[i].prevalence=config_setting_get_float_elem(setting,i);

    // stderr printf to check we read correctly
    for (i=0;i<dipolecount;i++)
    fprintf(stderr,"Dipole: %d Length: %f Prevalence: %f\n",i,dipoles[i].length, dipoles[i].prevalence);
    
    // above doesn't do anything currently - not sure whether I lost the code
    // at some point?

    config_lookup_bool(cf,"ConstrainToX",&ConstrainToX);

//    config_lookup_float(cf,"DipoleFraction",&dipole_fraction); DEPRECATED
//    - now done with Dipoles above

    config_lookup_int(cf,"DipoleCutOff",&DipoleCutOff);

    config_lookup_int(cf,"MCEqmSteps",&MCEqmSteps);

    config_lookup_int(cf,"MCMegaSteps",&MCMegaSteps); 
    config_lookup_float(cf,"MCMoves",&MCMegaMultiplier);

    MCMinorSteps=(int)((float)X*(float)Y*(float)Z*MCMegaMultiplier);

// Simulation display / calculation flags
    config_lookup_bool(cf,"DisplayDumbTerminal",&DisplayDumbTerminal);
    config_lookup_bool(cf,"CalculateRecombination",&CalculateRecombination);
    config_lookup_bool(cf,"CalculateRadialOrderParameter",&CalculateRadialOrderParameter);
      
    config_lookup_bool(cf,"CalculatePotential",&CalculatePotential);
    config_lookup_bool(cf,"CalculateEfield",&CalculateEfield);
    
    config_lookup_bool(cf,"SaveDipolesSVG",&SaveDipolesSVG);
    config_lookup_bool(cf,"SaveDipolesPNG",&SaveDipolesPNG);
    config_lookup_bool(cf,"SaveDipolesXYZ",&SaveDipolesXYZ);

    fprintf(stderr,"Config loaded. \n");
}

// 3-Vector dot-product... hand coded, should probably validate against
// a proper linear albegra library
// + also check generated machine code that it unrolls nicely + etc.
static float dot(struct dipole *a, struct dipole *b)
{
    int D;
    float sum=0.0;

    sum+=a->x*b->x;
    sum+=a->y*b->y;
    sum+=a->z*b->z;

    return(sum);
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

// choose random vector on cubic face (X) <001>
static void random_X_point(struct dipole *p)
{
    int i;
    int x=0,y=0,z=0;

    i=rand_int(6);

    switch(i)
    {
        case 0:
            x=1;
        break;
        case 1:
            x=-1;
        break;
        case 2:
            y=1;
        break;
        case 3:
            y=-1;
        break;
        case 4:
            z=1;
        break;
        case 5:
            z=-1;
        break;
    }

    // convert to floating point + pack back into passed structure
    p->x = (float) x;
    p->y = (float) y;
    p->z = (float) z;
}


