CC := gcc
SRC := ../src/
OUT := main
FLAGS := -std=c99 -mavx

all:
	make -C ../src
	$(CC) -c main.c -o main.o  $(FLAGS)
	$(CC) main.o $(SRC)output.o $(SRC)input.o $(SRC)fail.o -o $(OUT) $(FLAGS)
	
debug:
	$(CC) -c main.c -o main.o
	$(CC) -g main.o $(SRC)input.o $(SRC)fail.o -o $(OUT)

clean:
	rm *.o
	rm $(OUT)
