from mpmath import *

for i in range(8):
    ii = i+1
    z = ii*pi/4.0 + ii*pi/2.0*1j
    nu = (ii+2)*pi/4.0
    print ii,z,nu,besselj(nu,z),bessely(nu,z)

