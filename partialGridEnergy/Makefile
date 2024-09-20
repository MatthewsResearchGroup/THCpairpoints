CXX=g++
CXXFLAGS=-fopenmp -std=c++17 -g -O0 -DDEBUG -Wno-unused-result -I../docopt.cpp
BLAS=-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
LAPACK=
LDFLAGS=-fopenmp
LIBS=-lm -ldl ../docopt.cpp/libdocopt.a 

#
# Use optimized flags and libraries on M2
#
#ifeq ($(patsubst m3login%,m3login,$(shell hostname)),m3login)
#	CXX=g++
#	CXXFLAGS=-fopenmp -std=c++17 -march=haswell -g -O0 -DDEBUG 
#	BLAS=-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
#	LAPACK=
#	TORCH_DIR=/work/group/matthewsgroup/apps/libtorch
#else
#
# Otherwise find libtorch wherever it may be
#
#	TORCH_DIR=$(HOME)/libtorch
#endif

CXXFLAGS+=-I../common
CXXFLAGS+=-I../docopt.cpp
CXXFLAGS+=-I../randutils

#Hack for macOS "g++"
IS_APPLE=$(shell g++ --version 2>/dev/null | grep -q Apple && echo yes || echo no)
ifeq ($(IS_APPLE),yes)
	CXXFLAGS:=$(subst -fopenmp,,$(CXXFLAGS))
	LDFLAGS:=$(subst -fopenmp,,$(LDFLAGS))
endif

exe:=../bin/partialGridEnergy
src:=$(shell find . -iname "*.cxx") $(shell find ../common -iname "*.cxx")
objs:=$(addsuffix .o, $(basename $(src)))
deps:=$(join $(addsuffix .deps/,$(dir $(objs))),$(notdir $(objs:.o=.d)))

.PHONY: all
all: $(exe)

-include $(deps)

.PHONY: clean
clean:
	rm -f $(objs) $(exe)

$(exe): $(objs)
	@mkdir -p $(dir $(exe))
	$(CXX) $(LDFLAGS) -o $(exe) $(objs) $(BLAS) $(LAPACK) $(LIBS)

%.o: %.cxx
	@mkdir -p $(dir $@).deps
	$(CXX) $(CXXFLAGS) -MT $@ -MP -MMD -MF $(dir $@).deps/$(notdir $(@:.o=.d)) -o $@ -c $<

