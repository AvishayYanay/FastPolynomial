

#ifndef FASTPOLY_RECURSIVE_H
#define FASTPOLY_RECURSIVE_H

#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>

#include <chrono>
#include <sys/resource.h>

using namespace std;
using namespace NTL;
using namespace chrono;


void poly_interpolate_zp_recursive(long degree, ZZ_p *X, ZZ_p *Y, ZZ_pX &P);

void poly_evaluate_zp_recursive(long degree, ZZ_pX &P, ZZ_p *X, ZZ_p *Y);



#endif //FASTPOLY_RECURSIVE_H
