#include <quadmath.h>

// use arb C library to provide arbitrary precision Bessel
// functions of complex argument and real order
// http://fredrikj.net/arb/acb_hypgeom.html
#include "acb_hypgeom.h"

void arf_set_float128(arf_t res, __float128 x)
{
  // conversion from quad precision float to arf type
  double d1, d2, d3;
  arf_t t;

  // scale results, in case quad-precision value has exponent
  // larger than can be handled by doubles used in conversion
  fmpz_t ef;
  int ei;
  long el;
  x = frexpq(x, &ei);
  el = (long)ei;
  
  d1 = x;
  d2 = (x - d1);
  d3 = ((x - d1) - d2);
  
  arf_init(t);
  arf_set_d(res, d1);
  arf_set_d(t, d2);
  arf_add(res, res, t, ARF_PREC_EXACT, ARF_RND_DOWN);
  arf_set_d(t, d3);
  arf_add(res, res, t, ARF_PREC_EXACT, ARF_RND_DOWN);
  arf_clear(t);
  
  fmpz_init(ef);
  fmpz_set_si(ef, el);
  arf_mul_2exp_fmpz(res, res, ef);

  fmpz_clear(ef);
  
}

__float128 arf_get_float128(const arf_t x)
{
  // conversion from arf type to quad precision float
  arf_t t1, t2, t3;
  double d1, d2, d3;
  __float128 res;
  
  // need to do same frexpq scaling (but on arf_t type -- not sure which function to use)
  
  arf_init(t1);
  arf_init(t2);
  arf_init(t3);
  
  arf_set_round(t1, x, 53, ARF_RND_DOWN);
  arf_sub(t3, x, t1, ARF_PREC_EXACT, ARF_RND_DOWN);
  arf_set_round(t2, t3, 53, ARF_RND_DOWN);
  arf_sub(t3, t3, t2, 53, ARF_RND_DOWN);

  d1 = arf_get_d(t1, ARF_RND_DOWN);
  d2 = arf_get_d(t2, ARF_RND_DOWN);
  d3 = arf_get_d(t3, ARF_RND_DOWN);
  
  res = ((__float128) d1) + ((__float128) d2) + ((__float128) d3);

  arf_clear(t1);
  arf_clear(t2);
  arf_clear(t3);

  // perform an ldexpq here to apply exponent to float128.
  
  return res;
}

__complex128 arb_J(__float128 gcc_nu, __complex128 gcc_z)
{

  // arb complex types
  acb_t acb_z, acb_res, acb_nu;

  // libquadmath complex
  __complex128 gcc_J;
  
  acb_init(acb_res);
  acb_init(acb_nu);
  acb_zero(acb_nu);
  acb_init(acb_z);
  
  // gcc_z -> acb_z (argument)
  arf_set_float128(arb_midref(acb_realref(acb_z)), __real__(gcc_z));
  arf_set_float128(arb_midref(acb_imagref(acb_z)), __imag__(gcc_z));

  // gcc_nu -> acb_nu (order)
  arf_set_float128(arb_midref(acb_realref(acb_nu)), gcc_nu);
  
  // compute bessel function
  acb_hypgeom_bessel_j(acb_res, acb_nu, acb_z, 159); // not sure about this precision

  //  arb_res -> gcc_K1
  __real__(gcc_J) = arf_get_float128(arb_midref(acb_realref(acb_res)));
  __imag__(gcc_J) = arf_get_float128(arb_midref(acb_imagref(acb_res)));
  
  acb_clear(acb_res);
  acb_clear(acb_nu);
  acb_clear(acb_z);

  return gcc_J;
}

__complex128 arb_Y(__float128 gcc_nu, __complex128 gcc_z)
{

  // arb complex types
  acb_t acb_z, acb_res, acb_nu;

  // libquadmath complex
  __complex128 gcc_Y;
  
  acb_init(acb_res);
  acb_init(acb_nu);
  acb_zero(acb_nu);
  acb_init(acb_z);

  // gcc_z -> acb_z (argument)
  arf_set_float128(arb_midref(acb_realref(acb_z)), __real__(gcc_z));
  arf_set_float128(arb_midref(acb_imagref(acb_z)), __imag__(gcc_z));

  // gcc_nu -> acb_nu (order)
  arf_set_float128(arb_midref(acb_realref(acb_nu)), gcc_nu);

  // compute bessel function
  acb_hypgeom_bessel_y(acb_res, acb_nu, acb_z, 159); // not sure about this precsion

  //  arb_res -> gcc_K1
  __real__(gcc_Y) = arf_get_float128(arb_midref(acb_realref(acb_res)));
  __imag__(gcc_Y) = arf_get_float128(arb_midref(acb_imagref(acb_res)));
  
  acb_clear(acb_res);
  acb_clear(acb_nu);
  acb_clear(acb_z);

  return gcc_Y;
}


