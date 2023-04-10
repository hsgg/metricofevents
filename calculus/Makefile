CXXFLAGS += -g -Wall -pedantic
CFLAGS += $(CXXFLAGS)
LDLIBS += -lginac -lcln

CC = $(CXX)

PDFLATEX = pdflatex
PDFVIEWER = evince


HDR_CPP = coordinate_system.h \
	  functions.h
HDR_C = mytime.h
SRC_CPP = metric2energy.cpp \
	  $(HDR_CPP:.h=.cpp)
SRC_C = $(HDR_C:.h=.c)

SRC = $(SRC_CPP) $(SRC_C)
OBJ = $(SRC_CPP:.cpp=.o) $(SRC_C:.c=.o)
PROG = metric2energy

DEPFILE = depend
PROGDEPFILE = progdepend


# all
.PHONY: all
all: $(PROGDEPFILE) $(PROG)

obj: $(OBJ)


# dependecies
$(DEPFILE): $(SRC)
	$(CXX) -MM $(SRC) > $(DEPFILE)

$(PROGDEPFILE): $(DEPFILE)
	# Technically, the dependecies-only should be removed.
	sed -e "s/\.o//g" -e "s/\.cpp/.o/g" -e "s/\.c/.o/g" -e "s/\.h/.o/g" \
		$(DEPFILE) > $(PROGDEPFILE)

-include $(PROGDEPFILE)
-include $(DEPFILE)


# cleaners
.PHONY: clean
clean:
	rm -f $(OBJ)
	rm -f $(PROG)
	rm -f $(DEPFILE) $(PROGDEPFILE)
