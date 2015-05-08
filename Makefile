#
# Copyright (c) 2012-2014 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

##################################################
# flags / settings for gfortran >= 4.6 compiler

DEBUG = -O1 -g -Wall -Wextra -fcheck=all -fno-realloc-lhs -std=f2008 -fall-intrinsics

# for MS-Windows comment out the following line, since OpenMP doesn't typically work with
# minGW. I think it might work with cygwin, though.
OMP = -fopenmp

PERF = -O3 -march=native #-static

# either put gfortran in your PATH variable, or specify the path here
F90 = gfortran

# these are different if using Intel Fortran compiler, etc. Should work for gfortran
CPP = -cpp
FREE = -free

##################################################
# user typically shouldn't have to ever modify anything below this point
##################################################

PERFLDFLAGS = $(PERF)

EXTERNAL = cbessel.o
HILEV = time.o laplace_hankel_solutions.o driver_io.o integration.o
OBJS =  constants.o $(EXTERNAL) types.o utility.o invlap.o  $(HILEV)

MAIN = driver.o

OPTOBJS = $(patsubst %.o,%.opt.o,$(OBJS) $(MAIN))
DEBUGOBJS = $(patsubst %.o,%.debug.o,$(OBJS) $(MAIN))

F90SRC=$(patsubst %.o,%.f90,$(OBJS) $(MAIN))

OUT = unconfined
DEBUGOUT = debug_$(OUT)

LD = $(F90) -fbacktrace

# can't normally link statically when using OpenMP (Linux/MacOS only)
# following line is an "unsupported" workaround
#LD += -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

####### default optimized (no debugging) target ##########################
# use OMP in optimized link step
driver: $(OPTOBJS)
	$(LD)  $(PERFLDFLAGS) $(OMP) -o $(OUT) $(OPTOBJS)

####### compiler debugging ###
debug_driver: $(DEBUGOBJS)
	$(LD) -o $(DEBUGOUT) $(DEBUGOBJS)

complex_bessel.mod:cbessel.f90

cbessel.debug.o: cbessel.f90
	$(F90) -c $(PERF) -o cbessel.debug.o cbessel.f90
cbessel.opt.o: cbessel.f90
	$(F90) -c $(PERF) -o cbessel.opt.o cbessel.f90

# explicityly include openMP in driver routine
driver.opt.o: driver.f90
	$(F90) -c $(PERF) $(DEFAULTS) $(OMP) -o $@ $<

####### rule for making optimized object files ############
%.opt.o: %.f90
	$(F90) -c $(PERF) $(DEFAULTS) -o $@ $<

####### rule for making debugging object files ############
%.debug.o: %.f90
	$(F90) -c $(DEBUG) $(DEFAULTS) -o $@ $<

# optimized objects
constants.opt.o constants.mod: constants.f90
types.opt.o types.mod: types.f90 constants.mod
invlap.opt.o invlap.mod: invlap.f90 constants.mod types.mod
utility.opt.o utility.mod: utility.f90 constants.mod
time.opt.o time.mod: time.f90 constants.mod types.mod
laplace_hankel_solutions.opt.o laplace_hankel_solutions.mod: \
 laplace_hankel_solutions.f90 constants.mod types.mod time.mod \
 utility.mod cbessel.mod
driver_io.opt.o driver_io.mod: driver_io.f90 constants.mod types.mod \
 utility.mod
integration.opt.o integration.mod: integration.f90 constants.mod types.mod
driver.opt.o: driver.f90 types.mod driver_io.mod constants.mod \
 laplace_hankel_solutions.mod invlap.mod integration.mod

# debug objects
constants.debug.o constants.mod: constants.f90
types.debug.o types.mod: types.f90 constants.mod
invlap.debug.o invlap.mod: invlap.f90 constants.mod types.mod
utility.debug.o utility.mod: utility.f90 constants.mod
time.debug.o time.mod: time.f90 constants.mod types.mod
laplace_hankel_solutions.debug.o laplace_hankel_solutions.mod: \
 laplace_hankel_solutions.f90 constants.mod types.mod time.mod \
 utility.mod cbessel.mod
driver_io.debug.o driver_io.mod: driver_io.f90 constants.mod types.mod \
 utility.mod
integration.debug.o integration.mod: integration.f90 constants.mod types.mod
driver.debug.o: driver.f90 types.mod driver_io.mod constants.mod \
 laplace_hankel_solutions.mod invlap.mod integration.mod


###### clean up #################################
clean:
	rm -f *.o *.mod $(OUT) $(DEBUGOUT) $(MATOUT)

