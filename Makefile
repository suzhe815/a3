CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic

.PHONY: clean

meina3: meina3.cpp vector.o matrix.o unit.o
	$(CXX) $(CXXFLAGS) -o meina3 vector.o matrix.o unit.o
	
matrix.o: matrix.cpp matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp

vector.o: vector.cpp vector.h
	$(CXX) $(CXXFLAGS) -c vector.cpp

clean:
	rm -f  test vector.o matrix.o
