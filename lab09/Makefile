CXX = g++
CXXFLAGS = -std=c++14 -Wall -pedantic -g -ggdb3 -O2 -lgsl -lgslcblas -lm

SOURCES = $(wildcard *.cpp)
OBJ = $(SOURCES:.cpp=.o)
HEADERS = $(wildcard *.h)

main: $(OBJ) $(HEADERS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJ) 
%.o: %.cpp $(HEADERS)
	$(CXX) -c $< $(CXXFLAGS) 
run:
	./main
plots: 
	gnuplot plot.sh
	
all:  main  run plots
check:
	valgrind ./main 
.PHONY: clean
clean:
	rm -f *.o main
	rm -f *.dat
	rm -f *.png
	rm -f *.gif
