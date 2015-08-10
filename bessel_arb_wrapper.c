#include <quadmath.h>  
#include "acb_hypgeom.h"

#include "mpfr.h"

void arf_set_float128(arf_t res, __float128 x)
{
    double d1, d2, d3;
    arf_t t;
    
    //int e;
    //x = frexpq(x, &e);
    
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

    // not sure how to apply the exponent from frexpq down here
    // it seems I need to use : arf_set_fmpz_2exp(arf_t y, const fmpz_t m, const fmpz_t e)
    // but not sure how to get arf and int into fmpz_t form for this.
}

__float128 arf_get_float128(const arf_t x)
{
    arf_t t1, t2, t3;
    double d1, d2, d3;
    __float128 res;

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

    return res;
}

__float128 arf_get_float128_via_mpfr(const arf_t x)
{
  __float128 res;
  mpfr_t t;

  mpfr_init(t);
  arf_get_mpfr(t, x, MPFR_RNDN);  
  res = mpfr_get_float128(t, MPFR_RNDN);
  mpfr_clear(t);
  
  return res;
}

void arf_set_float128_via_mpfr(arf_t res, __float128 x)
{
  mpfr_t t;
  int n;

  mpfr_init(t);
  // n is result of mpfr_check_range()
  n = mpfr_set_float128(t, x, MPFR_RNDN);
  arf_set_mpfr(res, t);
  mpfr_clear(t);
  
}

__complex128 quadprec_to_arb_and_back_K1(__complex128 gcc_z)
{

  // arb complex types
  acb_t acb_z, acb_res, acb_nu;

  // libquadmath complex
  __complex128 gcc_K1;
  
  acb_init(acb_res);
  acb_init(acb_nu);
  acb_one(acb_nu);
  acb_init(acb_z);
  
  // gcc_z -> arb_z
  arf_set_float128(arb_midref(acb_realref(acb_z)), __real__(gcc_z));
  arf_set_float128(arb_midref(acb_imagref(acb_z)), __imag__(gcc_z));

  // compute bessel function
  acb_hypgeom_bessel_k(acb_res, acb_nu, acb_z, 159);

  //  arb_res -> gcc_K1
  __real__(gcc_K1) = arf_get_float128(arb_midref(acb_realref(acb_res)));
  __imag__(gcc_K1) = arf_get_float128(arb_midref(acb_imagref(acb_res)));
  
  acb_clear(acb_res);
  acb_clear(acb_nu);
  acb_clear(acb_z);

  return gcc_K1;
}


