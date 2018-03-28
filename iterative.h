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

void poly_interpolate_zp_iterative(long degree, NTL::ZZ_p *X, NTL::ZZ_p *Y, NTL::ZZ_pX &P);

void poly_evaluate_zp_iterative(long degree, NTL::ZZ_pX &P, NTL::ZZ_p *X, NTL::ZZ_p *Y);

#endif //FASTPOLY_ITERATIVE_H
