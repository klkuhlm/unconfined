This is a modern Fortran (2003/2008) well-test simulator.  It performs
numerical Laplace-Hankel inversion, implementing the the main
unconfined approaches still in use today.  The program is free
software (MIT license), which can essentially be used, modified, or
redistributed for any purpose, given the license is left intact.

This code is a command-line utility, which reads a text input file and
writes a text datafile, formatted for simple plotting using available
software obtained elsewhere (e.g., MS-Excel, python matplotlib, or
gnuplot).  The code is accurate and relatively fast, using OpenMP to
execute in parallel on a multi-processor Linux or Mac computer (a
recent Intel compiler is needed to create parallel executables for
MS-Windows).

The input parameters are explained in input-explanation.txt

The solutions implemented include:
------------------------------------------
1) Mishra & Neuman (2010,2011) : Unsaturated/saturated flow to a partially
penetrating well. http://dx.doi.org/10.1029/2009WR008899  http://dx.doi.org/10.1029/2010WR010177

NB: The Mishra & Neuman solutions given in the WRR papers are somewhat
ill-behaved.  My code implements them in three different ways.

1a) One approach to solve M/N follows the Malama (2014)
  http://dx.doi.org/10.1002/2013WR014909 simplified formulation --
  replacing the no-flow boundary condition at the land surface with a
  "finiteness" boundary condition.  This solution has only been
  derived so far for the fully penetrating no wellbore storage case.

1b) A second approach to solve M/N discritizes the vadose zone using
  finite differences (in Laplace-Hankel space).  This approach works and
  can be used as a "check" on the algebra and mathematics in the
  closed-form Laplace-space approaches.  This allows partial
  penetration, but no wellbore storage for now.
  
1c) The third approach to solving M/N implements the solution listed
  in their paper directly (and naively). This approach fails for some
  combinations of parameters, and is often suffers from severe cancell-
  ation in the transition region, between early and late time.

2) Malama (2011) : Alternative linearization of the
moving water table boundary condition.  Basically an improvement on
Neuman (1974). http://dx.doi.org/10.1016/j.jhydrol.2010.11.007

3) Moench (2001,1995) : The hybrid water table boundary condition of
Moench (1995), but including the multiple delayed yield (α)
coefficients, as used in the large Cape Cod, Massachusetts pumping
test in USGS Water Supply Paper 1629. http://dx.doi.org/10.1111/j.1745-6584.1995.tb00293.x  http://pubs.usgs.gov/pp/pp1629/pdf/pp1629ver2.pdf

4) Neuman (1974,1972) : The standard moving water table solution used
by most hydrologists for well-test interpretation. http://dx.doi.org/10.1029/WR008i004p01031  http://dx.doi.org/10.1029/WR010i002p00303

5) Hantush (1961) : The confined solution which includes the effects
of partial penetration, but using a three-layer approach of Malama
(2011), rather than the typical finite cosine transform.

6) Theis (1935) : The confined fully penetrating solution, which all
other solutions build upon.

The code is distributed as a collection of Fortran source files and a
makefile.  On Linux/Unix/Mac platforms this is trivial to turn into a
command-line program, by simply going to the source directory and
typing:

make

An unoptimized debugging version (produces some warnings that can be
safely ignored) can be compiled via:

make debug_driver

On MS-Windows, you will need the mingw compilation environment (OpenMP
doesn't work under mingw, though -- so single thread only) or the
Intel Fortran compiler (which works and provides OpenMP as well).  I
have recently compiled it (Feb 2015) with the free mingw toolchain,
and can either provide assistance setting this up, or provide you with
a binary. There are a few settings in the makefile that must be changed
to get it to compile using mingw (see comments there)

Data
------------------------------------------

I am in the process of collecting unconfined pumping datasets for
benchmarking and conducting a "beauty pagent" between the different 
unconfined models.  I am currently working to include the following 
datasets:

Moench et al., 2001 (Cape Cod, Massachusets): a 320-gpm 72-hour
pumping test in a thick glacial outwash aquifer with many shallow
piezometers and screened observation wells (significan partial
penetration effects) the thinned PEST dataset is available from the
authors electronically, and in a report (USGS Professional Paper 1629
v2). No recovery data.

Wenzel, 1942 (Grand Island, Nebraska): a 540-gpm 48-hour pumping test
available from an old report (USGS Water Supply Paper 887).  This test
also has significant partial penetration effects. I have entered this
data into spreadsheet form. 36 hours of recovery data available. A
large number of observation piezometers and a few screened wells (82 +
pumping well) were monitored. The spreadsheets are available for
anyone to view on Google docs at:

https://docs.google.com/spreadsheet/ccc?key=0AlJMuEYu7Z-5dGJfdzBibk4zNDB4UG9DN1FpQ0FnX1E&usp=sharing

Bevan, 2002 (Borden, Ontario); a 12-gpm 168-hour (week long) test was
obtained electronically from the authors, published in Michael Bevan's
MS thesis. The aquifer was relatively thin, compared to the Cape Cod
and Grand Island tests. Moisture content data and geophysical data
were collected in the vadose zone before, during, and after
testing. Five days of recovery data were collected.

I will make all the data available as they are cleaned/prepared
for use in my inverse modeling exercise.

Kris Kuhlman (klkuhlm _at_ sandia _dot_ gov)
Feb, 2015

