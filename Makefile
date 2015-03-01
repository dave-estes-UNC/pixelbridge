SRCDIR=./src
IDIR=./src

CC=gcc
CPP=g++

CCFLAGS=-I$(IDIR)
CPPFLAGS=-I$(IDIR)
#LDFLAGS=--static
LDFLAGS=

ifneq ($(NO_HACKS), 1)
	CCFLAGS+=-DSUPRESS_EXCESS_RENDERING -DSKIP_COMPUTE_WHEN_SCALER_ZERO
	CPPFLAGS+=-DSUPRESS_EXCESS_RENDERING -DSKIP_COMPUTE_WHEN_SCALER_ZERO
endif

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

ifeq ($(PROFILE), 1)
	CCFLAGS+=-pg
	CPPFLAGS+=-pg
	LDFLAGS+=-pg
endif

CPPFLAGS+=-D__STDC_CONSTANT_MACROS

LDFLAGS+=-lz -lavutil -lavformat -lavcodec -lswscale `pkg-config --libs opencv`

ifeq ($(DEBUG), 1)
	CCFLAGS+=-g -DDEBUG
	CPPFLAGS+=-g -DDEBUG
	OBJDIR=build/objdir-debug
else
	CCFLAGS+= -O3 -ffast-math
	CPPFLAGS+= -O3 -ffast-math
	OBJDIR=build/objdir-release
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
	ClCoefficientPlanes.h \
	ClFrameVolume.h \
	ClInputVector.h \
	ClNddiDisplay.h \
	CoefficientMatrix.h \
	CoefficientPlanes.h \
	CostModel.h \
	DctTiler.h \
	ItTiler.h \
	FfmpegPlayer.h \
	FlatTiler.h \
	FrameVolume.h \
	CostModel.h \
	PixelBridgeFeatures.h \
	Configuration.h \
	Player.h \
	RandomPlayer.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = \
	BaseNddiDisplay.o \
	BlendingGlNddiDisplay.o \
	CachedTiler.o \
	DctTiler.o \
	ItTiler.o \
	FfmpegPlayer.o \
	RandomPlayer.o \
	FlatTiler.o \
	GlNddiDisplay.o \
	PixelBridgeMain.o \
	Rewinder.o

ifneq ($(NO_CL), 1)
_OBJ += \
	ClNddiDisplay.o
endif

OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))


$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	mkdir -p $(dir $@)
	$(CPP) -c -o $@ $< $(CPPFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(CCFLAGS)

all:	pixelbridge

pixelbridge: $(OBJDIR)/PixelBridgeMain.o $(OBJ)
	$(CPP) -o $(OBJDIR)/$@ $^ $(CPPFLAGS) $(LDFLAGS) 

.PHONY: clean

clean:
	rm -Rf build *~ core

debug:
	$(MAKE) $(MAKEFILE) DEBUG=1
