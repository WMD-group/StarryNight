SRCs= src/mt19937ar-cok.c                src/starrynight-config.c           src/starrynight-main.c             src/xorshift1024star.c \
	  src/starrynight-analysis.c         src/starrynight-lattice.c          src/starrynight-montecarlo-core.c  src/xorshift128plus.c

# Code compilation
starrynight: $(SRCs) 
	gcc -O4 -lm -lconfig -o starrynight src/starrynight-main.c

starrynight-openmp: ${SRCs} 
	gcc -O4 -lm -lconfig -fopenmp -o starrynight src/starrynight-main.c

starrynight-mac-openmp: ${SRCs}
	/usr/local/bin/gcc-4.8 -O4 -lm -lconfig -fopenmp -lgomp -o starrynight src/starrynight-main.c

profile: ${SRCs} 
	gcc -lm -lconfig -o starrynight src/starrynight-main.c -pg

all: starrynight

clean:
	rm starrynight *.pnm *.jpg *.gif *.avi *.svg 
	rm *.png 
	rm *.log 
	rm *.dat 
	rm *.xyz
# Make file magics to assist running jobs 

cx1:
	    # Local version of libconfig, within starrynight directory
	    gcc -Llibconfig-1.5/lib -Ilibconfig-1.5/lib \
			-O4 -lm -o starrynight src/starrynight-main.c libconfig-1.5/lib/.libs/libconfig.a
# module load intel-suite
cx1-icc: 
	icc -Llibconfig-1.5/lib -Ilibconfig-1.5/lib \
	-O4 -o starrynight src/starrynight-main.c libconfig-1.5/lib/.libs/libconfig.a -lm

parallel: starrynight
	seq 0 10 1000 | parallel  ./starrynight {}  

superparallel: starrynight
	awk 'BEGIN{for (i=0;i<1000;i=i+20) { for (j=0.0;j<=4.0;j=j+1.0) printf ("%f %f\n",i,j); }}' \
		| parallel --colsep ' ' ./starrynight {1} {2}  > aggregate.dat

parallel-annamaria: starrynight
	seq 0.9 0.02 1.0 | caffeinate parallel ./starrynight {} | sort -k2 -g > variance.dat

parallel-CageStrain: starrynight
	seq 0 0.5 3.0 | caffeinate parallel ./starrynight {} > landau.dat

parallel-T: starrynight
	seq 75 75 600 | parallel ./starrynight {} > starrynight-parallel-T.log
