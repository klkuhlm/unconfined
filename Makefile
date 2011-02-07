
EXTERNAL = cbessel.o  
KKFILES = laplace_hankel_solutions.o  
OBJS = $(EXTERNAL) constants.o invlap.o utility.o $(KKFILES)

MAIN = driver.o

OPTOBJS = $(patsubst %.o,%.opt.o,$(OBJS) $(MAIN))
DEBUGOBJS = $(patsubst %.o,%.debug.o,$(OBJS) $(MAIN))

F90SRC=$(patsubst %.o,%.f90,$(OBJS) $(MAIN))

OUT = slug
DEBUGOUT = debug_slug

LD = $(F90)

####### default optimized (no debugging) target ##########################
driver: $(OPTOBJS)
	$(LD)  $(PERFLDFLAGS) -o $(OUT) $(OPTOBJS)

####### compiler debugging ### 
##(no optimization, checks for out-of-bounds arrays and gives more warninngs, but still runs to completion)
debug_driver: $(DEBUGOBJS)
	$(LD) -o $(DEBUGOUT) $(DEBUGOBJS)

####### rule for making optimized object files ############
%.opt.o: %.f90
	$(F90) -c -cpp $(INTEL) $(PERF) -o $@ $<

####### rule for making debugging object files ############
%.debug.o: %.f90
	$(F90) -c -cpp $(INTEL) $(DEBUG) -o $@ $<

cbessel.opt.o complex_bessel.mod : cbessel.f90
constants.opt.o constants.mod : constants.f90
invlap.opt.o inverse_laplace_transform.mod : constants.mod invlap.f90 
utility.opt.o shared_data.mod utilities.mod: constants.mod inverse_laplace_transform.mod utility.f90
laplace_hankel_solutions.opt.o lap_hank_soln.mod: complex_bessel.mod inverse_laplace_transform.mod constants.mod shared_data.mod utilities.mod laplace_hankel_solutions.f90
driver.opt.o: lap_hank_soln.mod constants.mod shared_data.mod utilities.mod driver.f90

cbessel.debug.o complex_bessel.mod : cbessel.f90
constants.debug.o constants.mod : constants.f90
invlap.debug.o inverse_laplace_transform.mod : constants.mod invlap.f90 
utility.debug.o shared_data.mod utilities.mod: constants.mod inverse_laplace_transform.mod utility.f90
laplace_hankel_solutions.debug.o lap_hank_soln.mod: complex_bessel.mod inverse_laplace_transform.mod constants.mod shared_data.mod utilities.mod laplace_hankel_solutions.f90
driver.debug.o: lap_hank_soln.mod utilities.mod constants.mod shared_data.mod driver.f90


###### clean up #################################
clean:
	rm -f *.o *.mod $(OUT) $(DEBUGOUT) $(MATOUT)

