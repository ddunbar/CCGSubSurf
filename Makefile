ifndef SRCROOT
  SRCROOT := $(shell pwd)
endif
ifndef OBJROOT
  OBJROOT := $(SRCROOT)/obj
endif

CC := clang
LD := clang

# Default compile options
ifndef CFLAGS
CFLAGS := -arch x86_64
endif
ifndef OPTFLAGS
OPTFLAGS := -O3 -fstrict-aliasing
endif

# Compiler
ExtraFlags := \
	-pedantic -Wall -Wextra \
	-Wno-unused-parameter -Wno-missing-field-initializers
CommonFlags := $(CFLAGS) $(OPTFLAGS) $(ExtraFlags)
CompileOnlyFlags := -std=c99 -I$(SRCROOT)/src
CompileFlags := $(CompileOnlyFlags) $(CommonFlags)
LinkFlags := $(CommonFlags)
LinkFlags.glutTest := $(LinkFlags) \
	-framework OpenGL -framework GLUT
LinkFlags.perfTest := $(LinkFlags)

# Source files

Headers := $(wildcard $(SRCROOT)/*.h)
Sources.Lib := CCGSubSurf.c
Sources.Tests := $(Sources.Lib) QMesh.c

Sources.glutTest := $(Sources.Tests) glutTest.c
Objects.glutTest := $(Sources.glutTest:%.c=$(OBJROOT)/%.o)

Sources.perfTest := $(Sources.Tests) perfTest.c
Objects.perfTest := $(Sources.perfTest:%.c=$(OBJROOT)/%.o)

all: $(OBJROOT)/glutTest $(OBJROOT)/perfTest

$(OBJROOT)/glutTest: $(Objects.glutTest)
	$(LD) -o $@ $(Objects.glutTest) $(LinkFlags.glutTest)

$(OBJROOT)/perfTest: $(Objects.perfTest)
	$(LD) -o $@ $(Objects.perfTest) $(LinkFlags.perfTest)

$(OBJROOT)/%.o: $(SRCROOT)/tests/%.c $(Headers) Makefile $(OBJROOT)/.dir
	$(CC) -c -o $@ $< $(CompileFlags)
$(OBJROOT)/%.o: $(SRCROOT)/src/%.c $(Headers) Makefile $(OBJROOT)/.dir
	$(CC) -c -o $@ $< $(CompileFlags)

%/.dir:
	mkdir -p $* > /dev/null
	echo > $@
