IDIR =.

CC=gcc
CPP=g++

CCFLAGS=-I$(IDIR)
CPPFLAGS=-I$(IDIR)
#LDFLAGS=--static
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
		CPPFLAGS+=-DGL_GLEXT_PROTOTYPES
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

CPPFLAGS+=-D__STDC_CONSTANT_MACROS

LDFLAGS+=-lz -lavutil -lavformat -lavcodec -lswscale

ifeq ($(DEBUG), 1)
	CCFLAGS+=-g -DDEBUG
	CPPFLAGS+=-g -DDEBUG
	OBJDIR=../build/objdir-debug
else
	OBJDIR=../build/objdir-release
endif

_DEPS =	GlNddiDisplay.h \
	InputVector.h \
	NDimensionalDisplayInterfaceExtended.h \
	NDimensionalDisplayInterface.h \
	Rewinder.h \
	Tiler.h \
	BaseNddiDisplay.h \
	BlendingGlNddiDisplay.h \
	CachedTiler.h \
	ClCoefficientPlane.h \
	ClFrameVolume.h \
	ClInputVector.h \
	ClNddiDisplay.h \
	CoefficientMatrix.h \
	CoefficientPlane.h \
	CostModel.h \
	FfmpegPlayer.h \
	FlatTiler.h \
	FrameVolume.h \
	CostModel.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =	GlNddiDisplay.o \
	PixelBridgeMain.o \
	Rewinder.o \
	BaseNddiDisplay.o \
	BlendingGlNddiDisplay.o \
	CachedTiler.o \
	ClNddiDisplay.o \
	FfmpegPlayer.o \
	FlatTiler.o

OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))


$(OBJDIR)/%.o: %.cpp $(DEPS)
	mkdir -p $(dir $@)
	$(CPP) -c -o $@ $< $(CPPFLAGS)

$(OBJDIR)/%.o: %.c $(DEPS)
	mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(CCFLAGS)

all:	pixelbridge

pixelbridge: $(OBJDIR)/PixelBridgeMain.o $(OBJ)
	$(CPP) -o $(OBJDIR)/$@ $^ $(CPPFLAGS) $(LDFLAGS) 

.PHONY: clean

clean:
	rm -Rf ../build *~ core $(INCDIR)/*~ 

debug:
	$(MAKE) $(MAKEFILE) DEBUG=1
