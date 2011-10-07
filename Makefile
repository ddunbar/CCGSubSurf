ifndef SRCROOT
  SRCROOT := $(shell pwd)
endif
ifndef OBJROOT
  OBJROOT := $(SRCROOT)/obj
endif

CC := clang
LD := clang

# Compiler
ExtraFlags := -Wall -Wextra \
	-Wno-unused-parameter -Wno-missing-field-initializers
CompileFlags := $(CFLAGS) $(OPTFLGS) $(ExtraFlags)
LinkFlags := $(CompileFlags) \
	-framework OpenGL -framework GLUT

# Source files

Sources.Lib := CCGSubSurf.c
Sources.Tests := $(Sources.Lib) QMesh.c
Sources.glutTest := $(Sources.Tests) glutTest.c
Objects.glutTest := $(Sources.glutTest:%.c=$(OBJROOT)/%.o)

all: $(OBJROOT)/glutTest

$(OBJROOT)/glutTest: $(Objects.glutTest)
	$(LD) -o $@ $(Objects.glutTest) $(LinkFlags)

$(OBJROOT)/%.o: $(SRCROOT)/%.c $(OBJROOT)/.dir
	$(CC) -c -o $@ $< $(CompileFlags)

%/.dir:
	mkdir -p $* > /dev/null
	echo > $@
