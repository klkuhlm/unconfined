##################################################
# flags / settings for gfortran >= 4.6 compiler

DEBUG = -O0 -g -Wall -Wextra -fcheck=all
OMP = -fopenmp
PERF = -Ofast -mtune=native
DEFAULTS = -fdefault-real-8 -fdefault-integer-8
F90 = gfortran-4.7
CPP = -cpp
FREE = -free
PERFLDFLAGS = $(PERF)
##################################################

EXTERNAL = cbessel.o
HILEV = time.o laplace_hankel_solutions.o driver_io.o integration.o
OBJS =  constants.o $(EXTERNAL) types.o invlap.o utility.o $(HILEV)

MAIN = driver.o

OPTOBJS = $(patsubst %.o,%.opt.o,$(OBJS) $(MAIN))
DEBUGOBJS = $(patsubst %.o,%.debug.o,$(OBJS) $(MAIN))

F90SRC=$(patsubst %.o,%.f90,$(OBJS) $(MAIN))

OUT = unconfined
DEBUGOUT = debug_$(OUT)

# can't link statically when using OpenMP!!!
LD = $(F90) -fbacktrace # -static

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

