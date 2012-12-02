#CXX=/usr/lib/colorgcc/bin/g++
#CC=/usr/lib/colorgcc/bin/gcc
CXX=g++
CC=gcc
FLAGS_QUICKCOMPILE=-Wall -pedantic -O2 -march=native
CFLAGS += $(FLAGS_QUICKCOMPILE)
#CFLAGS+=-O2 -march=pentium4 -mfpmath=sse
LDFLAGS += -lginac -lcln

EXT=exe


# Metrics:
#METRIC=../metrics/schwarzschild.metric
#METRIC=../metrics/robertson-walker-schwarzschild.metric
METRIC=../metrics/kerr-newman.metric
#METRIC=../metrics/schwarzschild-isotropic.metric


CS=../calculus/coordinate_system.cpp
CS_H=$(CS:.cpp=.h)
CS_O=$(CS:.cpp=.o)


.PHONY: all plot wrong spec programs datafiles globalplot globalspec
all: programs globalplot globalspec


plot: gnuplot.gplt plot.dat sphere.dat sphere2.dat
	gnuplot -persist gnuplot.gplt

wrong: plot gnuplot-wrong.gplt
	gnuplot -persist gnuplot-wrong.gplt

spec: plot gnuplot-spec.gplt
	gnuplot -persist gnuplot-spec.gplt

globalplot: gnuplot-globalplot.gplt globalplot.dat
	gnuplot -persist $<

globalspec: gnuplot-globalspec.gplt globalspec.dat
	gnuplot -persist $<


programs: algebra2src.$(EXT) single.$(EXT) sphere.$(EXT) many.$(EXT)

datafiles: plot.dat sphere.dat sphere2.dat


# algebra2src
algebra2src.$(EXT): algebra2src.o $(CS_O)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

algebra2src.o: algebra2src.cpp $(CS_H) Makefile
	$(CXX) $(CFLAGS) -c -o $@ $<

$(CS_O): $(CS) $(CS_H)
	$(MAKE) -C ../calculus obj CXXFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"


# single
single.$(EXT): single.o numerics.o init.o metric.o particle.o emfield.o spectrum.o misclib.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

single.o: single.cpp numerics.h global.h init.h metric.h myfloat.h particle.h misclib.h
	$(CXX) $(CFLAGS) -c -o $@ $<

# many
many.$(EXT): many.o numerics.o init.o metric.o particle.o emfield.o spectrum.o misclib.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

many.o: many.cpp numerics.h global.h init.h metric.h myfloat.h particle.h misclib.h
	$(CXX) $(CFLAGS) -c -o $@ $<

# *.o
numerics.o: numerics.cpp numerics.h global.h metric.h myfloat.h particle.h
	$(CXX) $(CFLAGS) -c -o $@ $<

init.o: init.cpp init.h myfloat.h
	$(CXX) $(CFLAGS) -c -o $@ $<

spectrum.o: spectrum.cpp spectrum.h myfloat.h misclib.h
	$(CC) $(CFLAGS) -c -o $@ $<

particle.o: particle.cpp particle.h global.h myfloat.h
	$(CXX) $(CFLAGS) -c -o $@ $<

metric.o: metric.cpp metric.h global.h myfloat.h
	$(CXX) $(FLAGS_QUICKCOMPILE) -c -o $@ $<
metric.h metric.cpp: algebra2src.$(EXT) $(METRIC)
	./$< "$(METRIC)"

emfield.o: emfield.cpp emfield.h metric.h myfloat.h
	$(CXX) $(CFLAGS) -c -o $@ $<

misclib.o: misclib.cpp misclib.h metric.h global.h myfloat.h init.h
	$(CXX) $(CFLAGS) -c -o $@ $<

# sphere
sphere.$(EXT): sphere.o init.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
sphere.o: sphere.cpp global.h init.h myfloat.h
	$(CXX) $(CFLAGS) -c -o $@ $<


# plot.dat
plot.dat: single.$(EXT)
	time ./$< input

# sphere.dat, sphere2.dat
sphere.dat sphere2.dat: sphere.$(EXT)
	./$<

# globalspec.dat
globalspec.dat: many.$(EXT)
	rm globalplot.dat/plot*
	./$<

.PHONY: globalplot.dat
globalplot.dat: many.$(EXT) globalspec.dat
	mkdir -p $@

gnuplot-globalplot.gplt: globalplot.dat
	echo "set title \"Whohoo\"" >$@
	echo "set polar" >>$@
	echo "set grid polar" >>$@
	echo -n "plot " >>$@
	for dat in globalplot.dat/*.dat; do \
		echo -n "\"$${dat}\"," >>$@; \
	done
	echo "" >>$@
	sed -i -e 's/,$$//' $@

# cleaners
.PHONY: clean realclean
clean:
	rm -f *.o

realclean: clean
	rm -f algebra2src.exe single.exe sphere.exe many.exe
	rm -f metric.h metric.cpp
	rm -f dtau.dat plot.dat sphere.dat sphere2.dat wrong.dat spec.dat
	@make -C ../calculus clean
