IDIR =.

CC=gcc
CPP=g++

CCFLAGS=-I$(IDIR)
CPPFLAGS=-I$(IDIR)
LDFLAGS=

ifeq ($(NO_OMP), 1)
	CCFLAGS+=-DNO_OMP
	CPPFLAGS+=-DNO_OMP
else
	CCFLAGS+=-fopenmp
	CPPFLAGS+=-fopenmp
endif

ifeq ($(NO_GL), 1)
	CCFLAGS+=-DNO_GL
	CPPFLAGS+=-DNO_GL
else
	ifeq ($(shell uname),Darwin)
		LDFLAGS+=-framework Carbon -framework OpenGL -framework GLUT
	else
		LDFLAGS+=-lGL -lglut -lGLU
	endif
endif

ifeq ($(NO_CL), 1)
	CCFLAGS+=-DNO_CL
	CPPFLAGS+=-DNO_CL
else
	ifeq ($(shell uname),Darwin)
		LDFLAGS+=-framework OpenCL
	else
		LDFLAGS+=-lOpenCL
	endif
endif

#LDFLAGS+=-lavformat.52.84.0 -lavcodec.52.93.0
LDFLAGS+=-lavformat -lavcodec

ifeq ($(DEBUG), 1)
	CCFLAGS+=-g -DDEBUG
	CPPFLAGS+=-g -DDEBUG
	OBJDIR=objdir-debug
else
	OBJDIR=objdir-release
endif

_DEPS =	NDimensionalDisplayInterface.h \
	BaseNddiDisplay.h \
	GlNddiDisplay.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =	BaseNddiDisplay.o \
	GlNddiDisplay.o

OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))


$(OBJDIR)/%.o: %.cpp $(DEPS)
	mkdir -p $(dir $@)
	$(CPP) -c -o $@ $< $(CPPFLAGS)

$(OBJDIR)/%.o: %.c $(DEPS)
	mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(CCFLAGS)

all:	pixelbridge

pixelbridge: $(OBJDIR)/PixelBridgeMain.o $(OBJ)
	$(CPP) $(LDFLAGS) -o $(OBJDIR)/$@ $^ $(CPPFLAGS)

.PHONY: clean

clean:
	rm -Rf objdir-* *~ core $(INCDIR)/*~ 

debug:
	$(MAKE) $(MAKEFILE) DEBUG=1
