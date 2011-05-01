##################################################
### flags / settings for intel >= 12.0 compiler
### version 12 still doesn't seem to handle complex extended-precision
### hyperbolic trig functions (although these are part of F2008)
##
##DEBUG = -O0 -g -warn all -check all -traceback
##OMP = -openmp
##PERF = -O3 -xHOST -ipo 
##REAL = -real-size 64
##F90 = ifort
##CPP = -cpp
##FREE = -free
##PERFLDFLAGS = ${PERF}
####################################################


##################################################
# flags / settings for gfortran >= 4.6 compiler

DEBUG = -O0 -g -Wall -Wextra -fbacktrace -fwhole-file
DEBUG += -frange-check -fcheck=all #-finit-integer=-999 -finit-real=snan ## -ffpe-trap=invalid
OMP = -fopenmp
PERF = -O2 -march=native -fwhole-file 
#PERF += $(OMP)
DEFAULTS = -fdefault-real-8 -fdefault-integer-8 # for constants like 1.0, 2.0, etc.
F90 = gfortran-4.7
CPP = -cpp
FREE = -free
PERFLDFLAGS = $(PERF)
##################################################

EXTERNAL = cbessel.o
HILEV = time.o laplace_hankel_solutions.o driver_io.o integration.o
OBJS = $(EXTERNAL) constants.o types.o invlap.o utility.o $(HILEV)

MAIN = driver.o

OPTOBJS = $(patsubst %.o,%.opt.o,$(OBJS) $(MAIN))
DEBUGOBJS = $(patsubst %.o,%.debug.o,$(OBJS) $(MAIN))

F90SRC=$(patsubst %.o,%.f90,$(OBJS) $(MAIN))

OUT = unconfined
DEBUGOUT = debug_$(OUT)

LD = $(F90)

####### default optimized (no debugging) target ##########################
# only the driver routine has openmp directives
driver: $(OPTOBJS)
	$(LD)  $(PERFLDFLAGS) $(OMP) -o $(OUT) $(OPTOBJS)


####### compiler debugging ### 
##(no optimization, checks for out-of-bounds arrays and gives more warninngs, but still runs to completion)
debug_driver: $(DEBUGOBJS)
	$(LD) -o $(DEBUGOUT) $(DEBUGOBJS)


# always compile external libray with debugging off (and no default-real-8)
cbessel.debug.o: cbessel.f90
	$(F90) -c $(PERF) -o cbessel.debug.o cbessel.f90
cbessel.opt.o: cbessel.f90
	$(F90) -c $(PERF) -o cbessel.opt.o cbessel.f90


####### rule for making optimized object files ############
%.opt.o: %.f90
	$(F90) -c $(PERF) $(DEFAULTS) -o $@ $<

####### rule for making debugging object files ############
%.debug.o: %.f90
	$(F90) -c $(DEBUG) $(DEFAULTS) -o $@ $<

constants.opt.o constants.mod : constants.f90
integration.opt.o integration.mod : integration.f90 constants.mod types.mod
invlap.opt.o inverse_laplace_transform.mod: invlap.f90 constants.mod types.mod 
utility.opt.o utilities.mod: utility.f90 constants.mod 
laplace_hankel_solutions.opt.o laplace_hankel_solution.mod: \
 laplace_hankel_solutions.f90 constants.mod types.mod time_mod.mod
driver.opt.o: laplace_hankel_solution.mod constants.mod utilities.mod io.mod driver.f90
time.opt.o time_mod.mod: time.f90 constants.mod types.mod
types.opt.o types.mod: types.f90 constants.mod
driver_io.opt.o io.mod: driver_io.f90 constants.mod types.mod utilities.mod 
driver.opt.o: driver.f90 types.mod constants.mod integration.mod laplace_hankel_solution.mod \
 inverse_laplace_transform.mod 

constants.debug.o constants.mod : constants.f90
integration.debug.o integration.mod : integration.f90 constants.mod types.mod
invlap.debug.o inverse_laplace_transform.mod: invlap.f90 constants.mod types.mod 
utility.debug.o utilities.mod: utility.f90 constants.mod 
laplace_hankel_solutions.debug.o laplace_hankel_solution.mod: \
 laplace_hankel_solutions.f90 constants.mod types.mod time_mod.mod
driver.debug.o: laplace_hankel_solution.mod constants.mod utilities.mod io.mod driver.f90
time.debug.o time_mod.mod: time.f90 constants.mod types.mod
types.debug.o types.mod: types.f90 constants.mod
driver_io.debug.o io.mod: driver_io.f90 constants.mod types.mod utilities.mod 
driver.debug.o: driver.f90 types.mod constants.mod integration.mod laplace_hankel_solution.mod \
 inverse_laplace_transform.mod 

###### clean up #################################
clean:
	rm -f *.o *.mod $(OUT) $(DEBUGOUT) $(MATOUT)

