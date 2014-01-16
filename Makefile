starrynight: starrynight.c
	gcc -O3 -lm -o starrynight starrynight.c

all: starrynight

clean:
	rm starrynight
