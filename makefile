# For the Upgraded main.cpp in the folder v0.1.10/
# cpp compiler to use
CXX= g++

# source files for the program, along with object (.o) files (ext change only)
SRCS := main.cpp NodeTree.cpp Particle.cpp simulator.cpp \
        commandline.cpp gui.cpp OrbitalMechanics.cpp
# SRC := main.cpp
# SRCS := $(wildcard *.cpp)
SRCDIR := src
OBJDIR := obj
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS.o=.d)
BINS := $(SRCS:.cpp=)

# Compiler flags meant for all files
CPPFLAGS += -std=c++11 -fmax-errors=2 -Werror
# This line lets us choose which file needs specific flags
CPPFLAGS += $(CPPFLAGS-$@)
# adding a special flag for the executable
CPPFLAGS-primus += -lsfml-graphics -lsfml-system -lsfml-window -lboost_program_options

ifdef fast
CPPFLAGS += -O3
endif

all: primus

EXEOBJS = $(addprefix $(OBJDIR)/, $(OBJS))

# main program to build, and dependencies, along with compile statement
primus: $(EXEOBJS)
	$(CXX) -o primus $(EXEOBJS) $(CPPFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
		$(CXX) $(CPPFLAGS) -c $< -o $@

# depend: .depend
#
# .depend: $(SRCDIR)/$(SRCS)
# 	$(RM) ./.depend
# 	$(CXX) $(CPPFLAGS) -MM src/$^>>./.depend;

# make clean removes executable program and all generated .o files.
.PHONY: clean
clean:
			 $(RM) primus $(EXEOBJS)

# I broke the dependency creator bit when I moved to the multifile system.
# include .depend
# DO NOT DELETE
