/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
static int rand_int(int SPAN);
static void gen_neighbour();
static double site_energy(int x, int y, int z, struct dipole *newdipole, struct dipole *olddipole);
static void MC_moves(int moves);
static void MC_move();
static void MC_move_openmp();

static int rand_int(int SPAN) // TODO: profile this to make sure it runs at an OK speed.
{
    return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}


// The following code builds a neighbour list (of the delta dx,dy,dzs) for
// speedy evaluation of energy; results in a speedup as it avoids the for loops
// + 'ifs' during Monte Carlo; instead you just pull the deltas from the lookup
// table

enum {MAXNEIGHBOURS=10000};
struct {
    int dx;
    int dy;
    int dz;
    float d;
} neighbours[MAXNEIGHBOURS];
int neighbour=0; //count of neighbours

static void gen_neighbour()
{

    int dx,dy,dz=0;
    float d;

    for (dx=-DipoleCutOff;dx<=DipoleCutOff;dx++)
        for (dy=-DipoleCutOff;dy<=DipoleCutOff;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-DipoleCutOff;dz<=DipoleCutOff;dz++) //NB: conditional zDipoleCutOff to allow for 2D version
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                d=sqrt((float) dx*dx + dy*dy + dz*dz); //that old chestnut; distance in Euler space

                if (d>(float)DipoleCutOff) continue; // Cutoff in d

                // store precomputed list of neighbours
                neighbours[neighbour].dx=dx; neighbours[neighbour].dy=dy; neighbours[neighbour].dz=dz;
                neighbours[neighbour].d=d;
                neighbour++;

                if (neighbour>MAXNEIGHBOURS) // bounds check
                {
                    fprintf(stderr,"Run out of space for the neighbour list with %d neighbours. FAILING TO EXIT!\n\n", neighbour);
                    exit(-1);
                }
            }
    fprintf(stderr,"\nNeighbour list generated: %d neighbours found with DipoleCutOff=%d.\n",neighbour,DipoleCutOff);
}


// Calculate change in site energy of changing from olddipole -> newdipole
static double site_energy(int x, int y, int z, struct dipole *newdipole, struct dipole *olddipole)
{
    int dx,dy,dz=0;
    float d;
    double dE=0.0;
    struct dipole *testdipole, n;

    // This now iterates over the neighbour list of neighbours[0..neighbour]
    // Which contains all the precomputed dx,dy,dz for a spherical cutoff, and
    // the sqrt distances etc.

    // Sum over near neighbours for dipole-dipole interaction
    int i;
//#pragma omp parallel for private(dx,dy,dz,d,n) reduction(+:dE) schedule(static,1)
    // NB: Works, but only modest speed gains!
    for (i=0;i<neighbour;i++)
    {
        // read in dirn to neighbours + precomputed values
        dx=neighbours[i].dx; dy=neighbours[i].dy; dz=neighbours[i].dz;
        d=neighbours[i].d;

        testdipole=& lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];

        n.x=(float)dx/d; n.y=(float)dy/d; n.z=(float)dz/d; //normalised diff. vector

        //True dipole like
        dE+= (olddipole->length * testdipole->length) * 
            (
             ( dot(newdipole,testdipole) - 3*dot(&n,newdipole)*dot(&n,testdipole) ) -
             ( dot(olddipole,testdipole) - 3*dot(&n,olddipole)*dot(&n,testdipole) ) 
            ) / (d*d*d); 

        // Ferroelectric / Potts model - vector form
        //            dE+= - Dipole * dot(newdipole,testdipole) / (d*d*d)
        //                + Dipole * dot(olddipole,testdipole) / (d*d*d);

        // Now reborn as our cage-strain term!
        if ((dx*dx+dy*dy+dz*dz)==1) //only nearest neighbour
            dE+= - CageStrain* dot(newdipole,testdipole)
                + CageStrain * dot(olddipole,testdipole); // signs to energetic drive alignment of vectors (dot product = more +ve, dE = -ve)

    }

    // Interaction of dipole with (unshielded) E-field
    dE+= + dot(newdipole, & Efield)
        - dot(olddipole, & Efield);
    //fprintf(stderr,"%f\n",dot(newdipole, & Efield));

    if (K>0.0) // rarely used anymore; 2D lattice epitaxial strain term
    {
        // along .x projection, squared
        n.x=1.0; n.y=0.0; n.z=0.0;
        dE +=   - K*fabs(dot(newdipole,&n))
            + K*fabs(dot(olddipole,&n));
        // along .y projection, squared
        n.x=0.0; n.y=1.0; n.z=0.0;
        dE +=   - K*fabs(dot(newdipole,&n))
            + K*fabs(dot(olddipole,&n));
    }

    // point charge at centre of space
    //    n.x=x-(X/2); n.y=y-(Y/2); n.z=z-(Z/2);
    //    dE += 1.0 * (dot(newdipole,&n) - dot(olddipole,&n) ) / ((x-X/2)^2 - (y-Y/2)^2 - (z-Z/2)^2);

    return(dE); 
}

static void MC_moves(int moves)
{
    int i;
    //moves/=8; //hard coded domain decomp.
    for (i=0;i<moves;i++)
        MC_move();
}

static void MC_move()
{
    int x, y, z;
    float dE=0.0;
    struct dipole newdipole, *olddipole;

    // Choose random dipole / lattice location

    x=rand_int(X);
    y=rand_int(Y);
    z=rand_int(Z);

    if (lattice[x][y][z].length==0.0) return; //dipole zero length .'. not present

    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well
    if (ConstrainToX)
        random_X_point(& newdipole); //consider any <100> vector
    else
        random_sphere_point(& newdipole);    

    newdipole.length = lattice[x][y][z].length; // preserve length / i.d. of dipole
    olddipole=& lattice[x][y][z];

    //calc site energy
    dE=site_energy(x,y,z, & newdipole,olddipole);

    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x][y][z].x=newdipole.x;
        lattice[x][y][z].y=newdipole.y;
        lattice[x][y][z].z=newdipole.z;
        //      lattice[x][y][z].length=newdipole.length; // never changes with current
        //      algorithms.

        ACCEPT++;
    }
    else
        REJECT++;
}

