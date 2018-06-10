all: binom binom_old monte

binom: hw4_binom.o
	g++ -o $@ $^ -std=c++11 -O2

binom_old: hw4_binom_old.o
	g++ -o $@ $^ -std=c++11 -O2

monte: hw4_monte.o
	g++ -o $@ $^ -std=c++11 -O2

hw4_binom.o: hw4_binom.cpp
	g++ -c hw4_binom.cpp -std=c++11 -O2

hw4_binom_old.o: hw4_binom_old.cpp
	g++ -c hw4_binom_old.cpp -std=c++11 -O2

hw4_monte.o: hw4_monte.cpp
	g++ -c hw4_monte.cpp -std=c++11 -O2

clean:
	rm *.o binom binom_old monte
