
##################################################
# flags / settings for gfortran >= 4.6 compiler

DEBUG = -O0 -g -Wall -Wextra -fwhole-file
DEBUG += -frange-check -fcheck=all 
OMP = -fopenmp
PERF = -Ofast -march=native 
PERF += $(OMP)
DEFAULTS = -fdefault-real-8 -fdefault-integer-8 
F90 = gfortran-4.7
CPP = -cpp
FREE = -free
PERFLDFLAGS = $(PERF)
##################################################

EXTERNAL = bessel.o
HILEV = time.o laplace_hankel_solutions.o driver_io.o integration.o
OBJS =  constants.o $(EXTERNAL) types.o invlap.o utility.o $(HILEV)

MAIN = driver.o

OPTOBJS = $(patsubst %.o,%.opt.o,$(OBJS) $(MAIN))
DEBUGOBJS = $(patsubst %.o,%.debug.o,$(OBJS) $(MAIN))

F90SRC=$(patsubst %.o,%.f90,$(OBJS) $(MAIN))

OUT = unconfined
DEBUGOUT = debug_$(OUT)

LD = $(F90) -fbacktrace

####### default optimized (no debugging) target ##########################
# only the driver routine has openmp directives
driver: $(OPTOBJS)
	$(LD)  $(PERFLDFLAGS) $(OMP) -o $(OUT) $(OPTOBJS)


####### compiler debugging ### 
##(no optimization, checks for out-of-bounds arrays and gives more warninngs, but still runs to completion)
debug_driver: $(DEBUGOBJS)
	$(LD) -o $(DEBUGOUT) $(DEBUGOBJS)


# always compile external libray with debugging off (and no default-real-8)
bessel.debug.o: bessel.f90 constants.mod
	$(F90) -c $(PERF) -o bessel.debug.o bessel.f90
bessel.opt.o: bessel.f90 constants.mod
	$(F90) -c $(PERF) -o bessel.opt.o bessel.f90


####### rule for making optimized object files ############
%.opt.o: %.f90
	$(F90) -c $(PERF) $(DEFAULTS) -o $@ $<

####### rule for making debugging object files ############
%.debug.o: %.f90
	$(F90) -c $(DEBUG) $(DEFAULTS) -o $@ $<

constants.opt.o constants.mod: constants.f90
types.opt.o types.mod: types.f90 constants.mod
invlap.opt.o invlap.mod: invlap.f90 constants.mod types.mod constants.mod \
 types.mod constants.mod types.mod
utility.opt.o utility.mod: utility.f90 constants.mod constants.mod \
 constants.mod constants.mod constants.mod constants.mod constants.mod
time.opt.o time.mod: time.f90 constants.mod types.mod
laplace_hankel_solutions.opt.o laplace_hankel_solutions.mod: \
 laplace_hankel_solutions.f90 constants.mod types.mod time.mod \
 utility.mod constants.mod constants.mod types.mod time.mod utility.mod \
 constants.mod types.mod time.mod utility.mod constants.mod \
 types.mod utility.mod complex_bessel.mod
driver_io.opt.o driver_io.mod: driver_io.f90 constants.mod types.mod \
 utility.mod constants.mod types.mod constants.mod types.mod
integration.opt.o integration.mod: integration.f90 constants.mod types.mod \
 constants.mod types.mod constants.mod utility.mod constants.mod
driver.opt.o: driver.f90 types.mod driver_io.mod constants.mod \
 laplace_hankel_solutions.mod invlap.mod integration.mod

constants.debug.o constants.mod: constants.f90
types.debug.o types.mod: types.f90 constants.mod
invlap.debug.o invlap.mod: invlap.f90 constants.mod types.mod constants.mod \
 types.mod constants.mod types.mod
utility.debug.o utility.mod: utility.f90 constants.mod constants.mod \
 constants.mod constants.mod constants.mod constants.mod constants.mod
time.debug.o time.mod: time.f90 constants.mod types.mod
laplace_hankel_solutions.debug.o laplace_hankel_solutions.mod: \
 laplace_hankel_solutions.f90 constants.mod types.mod time.mod \
 utility.mod constants.mod constants.mod types.mod time.mod utility.mod \
 constants.mod types.mod time.mod utility.mod constants.mod \
 types.mod utility.mod complex_bessel.mod
driver_io.debug.o driver_io.mod: driver_io.f90 constants.mod types.mod \
 utility.mod constants.mod types.mod constants.mod types.mod
integration.debug.o integration.mod: integration.f90 constants.mod types.mod \
 constants.mod types.mod constants.mod utility.mod constants.mod
driver.debug.o: driver.f90 types.mod driver_io.mod constants.mod \
 laplace_hankel_solutions.mod invlap.mod integration.mod


###### clean up #################################
clean:
	rm -f *.o *.mod $(OUT) $(DEBUGOUT) $(MATOUT)

