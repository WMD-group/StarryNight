starrynight: starrynight.c
	gcc -O4 -lm -lconfig -o starrynight starrynight.c

starrynight-mac-openmp: starrynight.c
	/usr/local/bin/gcc-4.8 -O4 -lm -lconfig -fopenmp -lgomp -o starrynight starrynight.c

parallel: starrynight
	seq 0 20 1000 | parallel  ./starrynight {}  | sort -k2 -g > T-dep.dat

all: starrynight

clean:
	rm starrynight *.pnm *.jpg *.gif *.avi *.svg *.png
