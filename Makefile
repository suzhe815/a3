CXX = g++ -g
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic

.PHONY: clean

meina3: meina3.o vector.o matrix.o unit.o
	$(CXX) $(CXXFLAGS) -o meina3 meina3.o vector.o matrix.o unit.o
 
meina3.o: meina3.cpp matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c meina3.cpp
	
test.o: vector.o matrix.o test.o
	$(CXX) $(CXXFLAGS) -o test vector.o matrix.o test.o
 
matrix.o: matrix.cpp matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp

vector.o: vector.cpp vector.h
	$(CXX) $(CXXFLAGS) -c vector.cpp

clean:
	rm -f  meina3.o vector.o matrix.o
