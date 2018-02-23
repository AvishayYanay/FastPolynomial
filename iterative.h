//
// Created by bush on 22/02/18.
//

#ifndef FASTPOLY_ITERATIVE_H
#define FASTPOLY_ITERATIVE_H

#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>

#include <chrono>
#include <sys/resource.h>

using namespace std;
using namespace NTL;
using namespace chrono;

void poly_interpolate_zp_iterative(long degree, ZZ_p *X, ZZ_p *Y, ZZ_pX &P);

void poly_evaluate_zp_iterative(long degree, ZZ_pX &P, ZZ_p *X, ZZ_p *Y);

#endif //FASTPOLY_ITERATIVE_H
