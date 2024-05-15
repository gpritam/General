#pragma once

#ifndef AllPeriodic1D_H
#define AllPeriodic1D_H

#include "General.h"
#include "Fourier_Quadrature.h"
#include <bits/stdc++.h>

#define TECPLOT
//#define VISIT

extern int Nx;

extern double xl, xr, Lambda, Lx;

extern double *x, *rhsPDE;

void CreatePDEPlan(const int Nx, const double xl, const double xr);
void DestroyPDEPlan();

void SolvePDE(double *&U, double (*ConstructOperator)(const int));

void ComputeError(double *U, double (*ExactFunction)(const double));
void WriteFile(double *U, double (*ExactFunction)(const double));
#endif
