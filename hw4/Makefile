all: binom binom2

binom: hw4_binom.o
	g++ -o $@ $^ -std=c++11 -O2

binom2: hw4_binom2.o
	g++ -o $@ $^ -std=c++11 -O2

hw4_binom.o: hw4_binom.cpp
	g++ -c hw4_binom.cpp -std=c++11 -O2

hw4_binom2.o: hw4_binom2.cpp
	g++ -c hw4_binom2.cpp -std=c++11 -O2

clean:
	rm *.o
