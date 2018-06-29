SRCDIR=./src
IDIR=./src

CC=gcc
CPP=g++

CCFLAGS=-std=c++11 -I$(IDIR)
CPPFLAGS=-std=c++11 -I$(IDIR) -I$(IDIR)/cereal/include
#LDFLAGS=--static
LDFLAGS=

ifneq ($(NO_HACKS), 1)
	CCFLAGS+=-DSUPRESS_EXCESS_RENDERING -DSKIP_COMPUTE_WHEN_SCALER_ZERO
	CPPFLAGS+=-DSUPRESS_EXCESS_RENDERING -DSKIP_COMPUTE_WHEN_SCALER_ZERO
endif

ifneq ($(NO_OMP), 1)
	CCFLAGS+=-fopenmp -DUSE_OMP
	CPPFLAGS+=-fopenmp -DUSE_OMP
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
		CCFLAGS+=-I/opt/AMDAPP/include/
		CPPFLAGS+=-I/opt/AMDAPP/include/
		LDFLAGS+=-L/opt/AMDAPP/lib/x86_64 -L/opt/AMDAPP/lib/x86 -lOpenCL
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
	CCFLAGS+= -O3 -ffast-math -DNDEBUG
	CPPFLAGS+= -O3 -ffast-math -DNDEBUG
	OBJDIR=build/objdir-release
endif

_DEPS = \
	CachedTiler.h \
	Configuration.h \
	DctTiler.h \
	FfmpegPlayer.h \
	FlatTiler.h \
	ItTiler.h \
	MultiDctTiler.h \
	PixelBridgeFeatures.h \
	Player.h \
	Rewinder.h \
	ScaledDctTiler.h \
	RandomPlayer.h \
	Tiler.h \
	nddi/BaseNddiDisplay.h \
	nddi/BlendingGlNddiDisplay.h \
	nddi/ClCoefficientPlanes.h \
	nddi/ClFrameVolume.h \
	nddi/ClInputVector.h \
	nddi/ClNddiDisplay.h \
	nddi/CoefficientMatrix.h \
	nddi/CoefficientPlanes.h \
	nddi/CostModel.h \
	nddi/Features.h \
	nddi/FrameVolume.h \
	nddi/GlNddiDisplay.h \
	nddi/InputVector.h \
	nddi/NDimensionalDisplayInterface.h \
	nddi/NDimensionalDisplayInterfaceExtended.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = \
	CachedTiler.o \
	DctTiler.o \
	FfmpegPlayer.o \
	FlatTiler.o \
	ItTiler.o \
	MultiDctTiler.o \
	PixelBridgeMain.o \
	RandomPlayer.o \
	Rewinder.o \
	ScaledDctTiler.o \
	nddi/BaseNddiDisplay.o \
	nddi/BlendingGlNddiDisplay.o \
	nddi/GlNddiDisplay.o

ifneq ($(NO_CL), 1)
_OBJ += \
	nddi/ClNddiDisplay.o
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
