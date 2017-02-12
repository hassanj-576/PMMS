FLAGS = -std=c99


all:fail.o input.o output.o

input:
	gcc -c input.o -o input.o $(FLAGS)

output:
	gcc -c output.c -o output.o $(FLAGS)

fail:
	gcc -c fail.c -o fail.o $(FLAGS)


clean:
	rm *.o
