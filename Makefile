CXX ?= g++

PROJECT_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CPPDIR       := $(PROJECT_ROOT)/cpp
BINDIR       := $(PROJECT_ROOT)/bin
OBJDIR       := $(PROJECT_ROOT)/build/obj

CPLEX_STUDIO_DIR ?= /opt/ibm/ILOG/CPLEX_Studio201
CPLEXDIR          := $(CPLEX_STUDIO_DIR)/cplex
CONCERTDIR        := $(CPLEX_STUDIO_DIR)/concert
SYSTEM            ?= x86-64_linux
LIBFORMAT         ?= static_pic

CPPFLAGS := -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include
CXXFLAGS ?= -O2
CXXFLAGS += -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG
LDFLAGS  += -L$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) \
            -L$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LDLIBS   := -lconcert -lilocplex -lcplex -lm -lpthread -ldl

OBJECTS := $(OBJDIR)/cmpuf.o \
           $(OBJDIR)/modele.o \
           $(OBJDIR)/tri.o \
           $(OBJDIR)/pmedian.o \
           $(OBJDIR)/transport.o

.PHONY: all clean cmpuf

all: cmpuf

cmpuf: $(BINDIR)/cmpuf

$(BINDIR)/cmpuf: $(OBJECTS) | $(BINDIR)
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(OBJDIR)/%.o: $(CPPDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJDIR) $(BINDIR):
	mkdir -p $@

clean:
	$(RM) -r $(PROJECT_ROOT)/build

-include $(wildcard $(OBJDIR)/*.d)
