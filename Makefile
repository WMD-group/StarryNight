# Code compilation
starrynight: starrynight-analysis.c   starrynight-config.c  starrynight-lattice.c  starrynight-main.c
	gcc -O4 -lm -lconfig -o starrynight starrynight-main.c

starrynight-openmp: starrynight-analysis.c   starrynight-config.c  starrynight-lattice.c  starrynight-main.c
	gcc -O4 -lm -lconfig -fopenmp -o starrynight starrynight-main.c

starrynight-mac-openmp: starrynight-analysis.c   starrynight-config.c  starrynight-lattice.c  starrynight-main.c 
	/usr/local/bin/gcc-4.8 -O4 -lm -lconfig -fopenmp -lgomp -o starrynight starrynight-main.c

profile: starrynight-analysis.c   starrynight-config.c  starrynight-lattice.c  starrynight-main.c
	gcc -lm -lconfig -o starrynight starrynight-main.c -pg

all: starrynight

clean:
	rm starrynight *.pnm *.jpg *.gif *.avi *.svg *.png
# Make file magics to assist running jobs 

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
	seq 0 75 600 | parallel ./starrynight {} > starrynight-parallel-T.log

