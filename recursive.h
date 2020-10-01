

#ifndef FASTPOLY_RECURSIVE_H
#define FASTPOLY_RECURSIVE_H

#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>

#include <chrono>
#include <sys/resource.h>


void poly_interpolate_zp_recursive(long degree, NTL::ZZ_p *X, NTL::ZZ_p *Y, NTL::ZZ_pX &P);

void poly_evaluate_zp_recursive(long degree, NTL::ZZ_pX &P, NTL::ZZ_p *X, NTL::ZZ_p *Y);



#endif //FASTPOLY_RECURSIVE_H
