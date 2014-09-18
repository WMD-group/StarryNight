/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
static void lattice_angle_log(FILE *log);
static double polarisation();
static double dipole_potential(int x, int y, int z);
static void lattice_potential_log(FILE *log);
void lattice_potential_XY(char * filename);
void lattice_potential_XYZ(char * filename);
static double lattice_energy_log(FILE *log);
double landau_order();

void outputpotential_png(char * filename);
void outputlattice_pnm(char * filename);
void outputlattice_ppm_hsv(char * filename);
void outputlattice_svg(char * filename);
void outputlattice_xyz(char * filename);
void outputlattice_xyz_overprint(char * filename);
void outputlattice_pymol_cgo(char * filename);
void outputlattice_dumb_terminal();

static void lattice_angle_log(FILE *log) // nb: 2D angle
{
    int x,y;
    double angle;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
        {
            angle=atan2(lattice[x][y][0].y, lattice[x][y][0].x);
            fprintf(log,"%f\n",angle);
        }
}

// takes advantage of the fact that the integral of the polarisation is equal
// to the total dipole moment of the dielectric
static double polarisation()
{
    double P=0.0;
    int x,y,z;
    struct dipole n;

    n.x=1.0; n.y=0.0; n.z=0.0;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
                P+=dot(&lattice[x][y][z],&n); //dipole response in direction of Efield

    return(P);
}

//Calculate dipole potential at specific location
static double dipole_potential(int x, int y, int z) 
{
    int dx,dy,dz=0;
    int MAX=10;
    double pot=0.0;
    float d;
    struct dipole r;

    for (dx=-MAX;dx<MAX;dx++)
        for (dy=-MAX;dy<MAX;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-MAX;dz<MAX;dz++)
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                r.x=(float)(dx); r.y=(float)(dy); r.z=(float)(dz);

                d=sqrt((float) r.x*r.x + r.y*r.y + r.z*r.z); //that old chestnut

                if (d>(float)MAX) continue; // Cutoff in d

                // pot(r) = 1/4PiEpsilon * p.r / r^3
                // Electric dipole potential
                pot+=dot(& lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z] ,& r)/(d*d*d);
            }
    return(pot);
}

//Calculate dipole potential at specific location
/*static double dipole_potential(int x, int y, int z) 
  {
  int dx,dy,dz;
  int MAX=10;
  double pot=0.0;
  float d;
  struct dipole r;

  for (dx=0;dx<X;dx++)
  for (dy=0;dy<X;dy++)
  for (dz=0;dz<Z;dz++)
  {
  if (x-dx==0 && y-dy==0 && z-dz==0)
  continue; //no infinities / self interactions please!

  r.x=(float)(x-dx); r.y=(float)(y-dy); r.z=(float)(z-dz);

  d=sqrt((float) r.x*r.x + r.y*r.y + r.z*r.z); //that old chestnut

//            if (d>(float)MAX) continue; // Cutoff in d

// pot(r) = 1/4PiEpsilon * p.r / r^3
// Electric dipole potential
pot+=dot(& lattice[dx][dy][dz],& r)/(d*d*d);
}
return(pot);
}
*/

//Calculates dipole potential along trace of lattice
static void lattice_potential_log(FILE *log)
{
    int x,y,z;
    double pot;

    y=Y/2; //trace across centre of material. I know, I know, PBCs.
    z=0;
    for (x=0;x<X;x++)
    {
        pot=0.0;
        for (y=0;y<Y;y++)
            pot+=dipole_potential(x,y,z);
        fprintf(log,"%d %f %f\n",x,pot/(double)Y,dipole_potential(x,Y/2,z));
    }

}

//Calculates dipole potential across XY lattice
void lattice_potential_XY(char * filename)
{
    int x,y;
    double pot;
    FILE *fo;
    fo=fopen(filename,"w");

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            fprintf(fo,"%d %d %f\n",x,y,dipole_potential(x,y,0));
}

//Calculates dipole potential across XYZ volume
void lattice_potential_XYZ(char * filename)
{
    int x,y,z;
    double pot;
    FILE *fo;
    fo=fopen(filename,"w");

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
                fprintf(fo,"%d %d %f\n",x,y,dipole_potential(x,y,z));
}



void outputpotential_png(char * filename)
{
    int i,k,pixel;
    FILE *fo;
    fo=fopen(filename,"w");

    fprintf (fo,"P2\n%d %d\n%d\n", X, Y, SHRT_MAX);

    for (i=0;i<X;i++)
    {
        for (k=0;k<Y;k++)
        {
            pixel=SHRT_MAX/2+(int)(SHRT_MAX*0.1*dipole_potential(i,k,0));

            // Bounds checking :^)
            if (pixel<0) pixel=0;
            if (pixel>SHRT_MAX) pixel=SHRT_MAX;

            fprintf(fo,"%d ",pixel);
        }
        fprintf(fo,"\n");
    }

}

/* This whole function defunct - no longer have angles...
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

// TODO: Calculate lattice electric field profile as a result of
// dipoles. Integrate out to full size of lattice? Seems a bit
// heavy handed. Same cut-offs as used in dipole calculation??


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
*/

double landau_order()
{
    int x,y,z;
    double landau=0.0;
    struct dipole orientation;

    orientation.x=0.0; orientation.y=0.0; orientation.z=0.0;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            { 
                orientation.x+=lattice[x][y][z].x;
                orientation.y+=lattice[x][y][z].y;
                orientation.z+=lattice[x][y][z].z;
            }

    landau=dot(&orientation,&orientation) / (double)(X*Y*Z * X*Y*Z); // u.u = |u|^2, 
    // so need to divide by N*N to put Landau parameter on range [0;1]
    return(landau);
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
            fprintf(fo,"%d ",(int)(SHRT_MAX*atan2(lattice[i][k][0].y,lattice[i][k][0].x)/(2*M_PI)));
        fprintf(fo,"\n");
    }

}

// Outputs a PPM bitmap of lattice dipole orientation on a HSV colourwheel
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
            h=M_PI+atan2(lattice[i][k][0].y,lattice[i][k][0].x); //Nb: assumes 0->2PI interval!
            v=0.5+0.4*lattice[i][k][0].z; //darken towards the south (-z) pole
            s=0.6-0.6*fabs(lattice[i][k][0].z); //desaturate towards the poles

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

            if (lattice[i][k][0].x == 0.0 && lattice[i][k][0].y == 0.0 && lattice[i][k][0].z == 0.0)
            { r=0.0; g=0.0; b=0.0; } // #FADE TO BLACK
            //zero length dipoles, i.e. absent ones - appear as black pixels

            fprintf(fo,"%c%c%c",(char)(254.0*r),(char)(254.0*g),(char)(254.0*b));
        }
    fclose(fo); //don't forget :^)
}

//Outputs an SVG file of pointing lattice dipoles; designed to overlay with PPM
//routine above
void outputlattice_svg(char * filename)
{
    int x,y;

    FILE *fo;
    fo=fopen(filename,"w");

    fprintf(fo,"<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" height=\"%d\" width=\"%d\">\n",X,Y);

    //our arrow marker...
    fprintf(fo," <marker id=\"triangle\" viewBox=\"0 0 10 10\" refX=\"7\" refY=\"5\" markerUnits=\"strokeWidth\" markerWidth=\"2\" markerHeight=\"2\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" /></marker>\n");

    //No markers...  marker-end=\"url(#triangle)\"

    for (x=0;x<X;x++) // care with X&Y - non-intuitive to get agreement with outputlattice_ppm_hsv for post-production overlaying
        for (y=0;y<Y;y++)
            fprintf(fo," <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:rgb(%d,%d,%d);stroke-width:0.17\" marker-end=\"url(#triangle)\" />\n",
                    y+0.5 + 0.4*lattice[x][y][0].y, 
                    x+0.5 + 0.4*lattice[x][y][0].x,
                    y+0.5 - 0.4*lattice[x][y][0].y,
                    x+0.5 - 0.4*lattice[x][y][0].x,
                    (int)((-lattice[x][y][0].z+1.0)*127.0),
                    (int)((-lattice[x][y][0].z+1.0)*127.0),
                    (int)((-lattice[x][y][0].z+1.0)*127.0)
                   );
    // invert z-axis, and scale to greyscale. Therefore alternates with
    // pointing up and down with background colour

    fprintf(fo,"</svg>\n");

    fclose(fo);
}

#define ZSCALE 5.0 // Scales Z-axis in Pymol xyz / CGO outputs

void outputlattice_xyz(char * filename)
{
    int x,y,z;
    float r=1.6/2; // half length of C-N molecule
    float d=4.0; // lattice size - for placing molecule
    // artificially small - to make molecules relatively bigger!
    // Nb: set to 3.0 to get pymol to draw bonds between aligned MA    
    FILE *fo;
    fo=fopen(filename,"w");
    fprintf(fo,"%d\n\n",X*Y*Z*2); //number of atoms... i.e. lattice sites times 2

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                fprintf(fo,"C %f %f %f\n",d*x+r*lattice[x][y][z].x, d*y+r*lattice[x][y][z].y, ZSCALE*(d*z)+r*lattice[x][y][z].z);
                fprintf(fo,"N %f %f %f\n",d*x-r*lattice[x][y][z].x, d*y-r*lattice[x][y][z].y, ZSCALE*(d*z)-r*lattice[x][y][z].z);
            }
}

void outputlattice_xyz_overprint(char * filename)
{
    int x,y,z;
    float r=6.0; // half length of C-N molecule
    float d=6.4; // lattice size - for placing molecule

    FILE *fo;
    fo=fopen(filename,"w");
    fprintf(fo,"%d\n\nC 0.0 0.0 0.0\n",1+(X*Y*Z)); //number of atoms...

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
                fprintf(fo,"N %f %f %f\n",r*lattice[x][y][z].x, r*lattice[x][y][z].y, r*lattice[x][y][z].z);
}

// Outputs Pymol CGO sphere primitives of lattice dipole orientation on a HSV colourwheel
void outputlattice_pymol_cgo(char * filename)
{
    float d=4.0; //to agree with XYZ file
    float a=1.4; // radius of sphere within above
    int x,y,z;
    float angle;

    float r,g,b; // RGB
    float h,s,v; // HSV
    float p,t,q,f; // intemediates for HSV->RGB conversion
    int hp;

    FILE *fo;
    fo=fopen(filename,"w");
    fprintf(fo,"from pymol.cgo import *\nfrom pymol import cmd\n");
    fprintf(fo,"obj = [ ALPHA, 0.7"); // Alpha = degree of transparency (1.0 = opaque, 0.0 = transparent)

    //Set Saturation + Value, vary hue
    s=0.6; v=0.8;


    for (x=0;x<X;x++) //force same ordering as SVG...
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                h=M_PI+atan2(lattice[x][y][z].y,lattice[x][y][z].x); //Nb: assumes 0->2PI interval!
                v=0.5+0.4*lattice[x][y][z].z; //darken towards the south (-z) pole
                s=0.6-0.6*fabs(lattice[x][y][z].z); //desaturate towards the poles

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

                if (lattice[x][y][z].x == 0.0 && lattice[x][y][z].y == 0.0 && lattice[x][y][z].z == 0.0)
                { r=0.0; g=0.0; b=0.0; } // #FADE TO BLACK
                //zero length dipoles, i.e. absent ones - appear as black pixels

                fprintf(fo,",COLOR, %f, %f, %f,\n",r,g,b);
                // sphere of colour centered at the site
                fprintf(fo,"SPHERE, %f, %f, %f, %f\n",x*d,y*d,ZSCALE*z*d,a*d/2);

                // makes a square (plane) of colour at the site
                //  - presently half works - odd triangle taken out of the square!
                /*            fprintf(fo,"BEGIN, TRIANGLE_STRIP, NORMAL, 0.0, 0.0, 1.0,\n");
                              fprintf(fo,"VERTEX, %f, %f, %f,",d*(x-0.5), d*(y-0.5),  ZSCALE*(d*z));
                              fprintf(fo,"VERTEX, %f, %f, %f,",d*(x+0.5), d*(y-0.5),  ZSCALE*(d*z));
                              fprintf(fo,"VERTEX, %f, %f, %f,",d*(x+0.5), d*(y+0.5),  ZSCALE*(d*z));
                              fprintf(fo,"VERTEX, %f, %f, %f,",d*(x-0.5), d*(y+0.5),  ZSCALE*(d*z));
                              fprintf(fo,"END");*/
            }
    fprintf(fo,"]\n");
    fprintf(fo,"cmd.load_cgo(obj,'battenberg')");

    fclose(fo); //don't forget :^)
}

float DMAX=55.0; //sensible starting value...

void outputlattice_dumb_terminal()
{
    const char * arrows="-\\|/-\\|/"; // "Dancing at angles"
    int x,y;
    float a;
    int z=0;
    float new_DMAX=0.0; //used to calibrate next colour scale, based on present maxima of data
    float variance=0.0; // sum of potential^2
    float mean=0.0;

    fprintf(stderr,"%*s%*s\n",X+3, "DIPOLES", (2*X)+4,"POTENTIAL"); //padded labels

    for (y=0;y<Y;y++)
    {
        for (x=0;x<X;x++)
        {
            a=atan2(lattice[x][y][z].y,lattice[x][y][z].x);
            a=a/(M_PI); //fraction of circle
            a=a+1.0; //map to [0,2]
            a=a+0.125; //I could tell you what this magic number is, but then I'd have to kill you.
            // OK -seriously, it's 45degrees/2 in our current basis, so that
            // the colours + text are centered _AROUND_ the cardinal
            // directions, not oscillating either side of N,NE,E... etc.
            if (a>2.0) a=a-2.0; //wrap around so values always show.
            a*=4; //pieces of eight
            fprintf (stderr,"%c[%d",27,31+((int)a)%8 ); // Sets colour of output routine
            if (a<4.0)                                  // makes colour bold / normal depending on arrow orientation
                fprintf(stderr,";7");
            char arrow=arrows[(int)a];
            if (lattice[x][y][z].z> sqrt(2)/2.0) arrow='o';
            if (lattice[x][y][z].z<-sqrt(2)/2.0) arrow='x';

            if (lattice[x][y][z].x==0.0 && lattice[x][y][z].y==0.0 && lattice[x][y][z].z==0.0) arrow='*'; 

            fprintf(stderr,"m%c %c[0m",arrow,27);  // prints arrow
            fprintf(stderr,"%c[37m%c[0m",27,27); //RESET

            //            fprintf(stderr,"%c ",arrows[(int)a]); // dumb - just black 'n'
            //            white
        }

        // OK - now potential plot :^)
        float potential;
        //        const char * density=".,:;o*O#"; //increasing potential density
        const char * density="012345689";
        fprintf(stderr,"    ");
        for (x=0;x<X;x++)
        {
            potential=dipole_potential(x,y,z);
            variance+=potential*potential;
            mean+=potential;

            if (fabs(potential)>new_DMAX)
                new_DMAX=fabs(potential); // used to calibrate scale - technically this changes
            //printf("%f\t",potential); //debug routine to get scale

            //fprintf(stderr,"%c[%d",27,31+((int)(8.0*fabs(potential)/DMAX))%8); //8 colours
            //fprintf(stderr,"%c[48;5;%d",27,17+(int)(214.0*fabs(potential)/DMAX)); // Xterm 256 color map - (16..231)
            fprintf(stderr,"%c[48;5;%d",27,232+12+(int)(12.0*potential/DMAX)); // Xterm 256 color map - shades of grey (232..255)
            // https://code.google.com/p/conemu-maximus5/wiki/AnsiEscapeCodes#xterm_256_color_processing_requirements

            //if (potential<0.0) // if negative
            //    fprintf(stderr,";7"); // bold

            a=atan2(lattice[x][y][z].y,lattice[x][y][z].x);
            a=a/(M_PI); //fraction of circle
            a=a+1.0; //map to [0,2]
            a=a+0.125; //I could tell you what this magic number is, but then I'd have to kill you.
            // OK -seriously, it's 45degrees/2 in our current basis, so that
            // the colours + text are centered _AROUND_ the cardinal
            // directions, not oscillating either side of N,NE,E... etc.
            if (a>2.0) a=a-2.0; //wrap around so values always show.
            a*=4; //pieces of eight
            char arrow=arrows[(int)a];  // selectss arrow
            if (lattice[x][y][z].z> sqrt(2)/2.0) arrow='o'; // override for 'up' (towards you - physics arrow style 'o')
            if (lattice[x][y][z].z<-sqrt(2)/2.0) arrow='x'; // and 'down' (away from you, physics arrow style 'x')


            fprintf(stderr,"m%c%c%c[0m",density[(int)(8.0*fabs(potential)/DMAX)],arrow,27);
        }

        fprintf(stderr,"\n");
    }
    mean=mean/(X*Y);
    variance=variance/(X*Y); 
    fprintf(stderr,"T: %d DMAX: %f new_DMAX: %f variance: %f mean: %f\n",T,DMAX,new_DMAX,variance,mean);
    fprintf(stdout,"T: %d DMAX: %f new_DMAX: %f variance: %f mean: %f\n",T,DMAX,new_DMAX,variance,mean);
    DMAX=(new_DMAX+DMAX)/2.0; // mean of old and new (sampled, noisy) value
    DMAX=new_DMAX; // infinite fast following - but leads to fluctuations at steady state
}
