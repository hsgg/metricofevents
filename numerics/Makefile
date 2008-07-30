#CXX=/usr/lib/colorgcc/bin/g++
#CC=/usr/lib/colorgcc/bin/gcc
CXX=g++
CC=gcc
FLAGS_QUICKCOMPILE=-Wall -pedantic
FLAGS=$(FLAGS_QUICKCOMPILE)
#FLAGS+=-O2 -march=pentium4 -mfpmath=sse
LIBS=-lginac -lcln

EXT=exe


# Metrics:
METRIC=../metrics/schwarzschild.metric
#METRIC=../metrics/robertson-walker-schwarzschild.metric
#METRIC=../metrics/kerr-newman.metric


CS=../calculus/coordinate_system.cpp
CS_H=$(CS:.cpp=.h)
CS_O=$(CS:.cpp=.o)


all: programs datafiles plot


plot: gnuplot.gplt plot.dat sphere.dat sphere2.dat
	gnuplot -persist gnuplot.gplt

wrong: plot gnuplot-wrong.gplt
	gnuplot -persist gnuplot-wrong.gplt


programs: algebra2src.$(EXT) numerics.$(EXT) sphere.$(EXT)

datafiles: plot.dat sphere.dat sphere2.dat


# algebra2src
algebra2src.$(EXT): algebra2src.o $(CS_O)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

algebra2src.o: algebra2src.cpp $(CS_H) Makefile
	$(CXX) $(FLAGS) -c -o $@ $<

$(CS_O): $(CS) $(CS_H)
	@make -C ../calculus


# numerics
numerics.$(EXT): numerics.o init.o metric.o particle.o emfield.o
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)
numerics.o: numerics.cpp global.h init.h metric.h myfloat.h particle.h
	$(CXX) $(FLAGS) -c -o $@ $<

# init.o
init.o: init.cpp init.h myfloat.h
	$(CC) $(FLAGS) -c -o $@ $<

# particle.o
particle.o: particle.cpp particle.h global.h myfloat.h
	$(CXX) $(FLAGS) -c -o $@ $<

# metric.o
metric.o: metric.cpp metric.h global.h myfloat.h
	$(CXX) $(FLAGS_QUICKCOMPILE) -c -o $@ $<
metric.h: algebra2src.$(EXT) $(METRIC)
	./$< "$(METRIC)"
metric.cpp: algebra2src.$(EXT) $(METRIC)
	./$< "$(METRIC)"

# emfield
emfield.o: emfield.cpp emfield.h metric.h myfloat.h
	$(CXX) $(FLAGS) -c -o $@ $<

# sphere
sphere.$(EXT): sphere.o init.o
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)
sphere.o: sphere.cpp global.h init.h myfloat.h
	$(CXX) $(FLAGS) -c -o $@ $<


# plot.dat
plot.dat: numerics.$(EXT)
	./$<

# sphere.dat
sphere.dat: sphere.$(EXT)
	./$<

# sphere2.dat
sphere2.dat: sphere.$(EXT)
	./$<




# cleaners
clean:
	rm -f *.o
	rm -f *.dat

realclean: clean
	rm -f algebra2src numerics sphere
	rm -f algebra2src.exe numerics.exe sphere.exe
	rm -f metric.h metric.cpp
	rm -f plot.dat sphere.dat sphere2.dat
	@make -C ../calculus clean
