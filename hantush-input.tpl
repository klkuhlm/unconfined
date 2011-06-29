2  1  T  T  T                :: quiet output?, model choice, dimensionless output?, timeseries?, piezometer?
2.0189D-2                      :: Q (volumetric pumping rate) [L^3/T]
6.5D-1   3.5D-1               :: l/b, d/b (normalized __depth__ to bottom/top of screened interval)
2.54D-2    2.54D-2                       :: rw, rc (radius of well casing and tubing)
1.0D+0                        :: gamma (dimensionless wellbore skin)
001  0.0D+0  1.0D+0            :: pumping well time behavior and parameters 
5.2669D+1                        :: b (initial saturated thickness)
1.225D-3  5.288D-1                :: Kr,kappa (radial K and ratio Kz/Kr)
3.766D-6  2.521D-1               :: Ss,Sy
2.0D0                        ::  beta Malama linearization parameter
2.95D+0  3.7D-1   2.0D+0  2.2D-1   2.0D+1      :: Mishra/Neuman 2010; a_c,a_k,  psi_a,psi_k,  L
14  1.0D-8  1.0D-9           :: deHoog invlap;  M,alpha,tol
8  6                         :: tanh-sinh quad;  2^k-1 order, # extrapollation steps
2  2  10  50                 :: G-L quad;  min/max zero split, # zeros to integrate, # abcissa/zero
timedata.dat  15.5            :: file to read time data from (and default t otherwise)
spacedata.dat  5.0          :: file to read space data from (and default r otherwise)
%%%%z%%%%  0.0D0  5  2.54D-2  20.0                  :: relative top obs well screen OR piezometer loc, bottom screen, quadrature order across screen (not for piezometer)
hantush-test-%%%%z%%%%.out
