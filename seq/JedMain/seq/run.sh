make 
#./main -t pattern_100x150.pgm -c pattern_100x150.pgm -n 7 -m 24 -L 0 -H 15
#valgrind --leak-check=yes ./main -m 100 -n 150 -i 900 -e 0.1 -k 10 -L 0 -H 100 -c pat1_100x150.pgm -t pat1_100x150.pgm
#./main -m 100 -n 150 -i 900 -e 0.1 -k 10 -L 0 -H 100 -c pat1_100x150.pgm -t pat1_100x150.pgm
./main -m 100 -n 150 -i 900 -e 0 -k 100 -L 0 -H 100 -c pat1_100x150.pgm -t pat1_100x150.pgm
#./main -n 150 -m 100 -i 42 -e 0.0001 -k 50 -c pat1_100x150.pgm -t pat1_100x150.pgm -r 1 -k 15 -L 0 -H 100:
#./main
rm main
