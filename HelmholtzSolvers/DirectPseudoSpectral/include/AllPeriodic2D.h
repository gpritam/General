#pragma once

#ifndef AllPeriodic2D_H
#define AllPeriodic2D_H

#include "General.h"
#include "Fourier_Quadrature.h"
#include <bits/stdc++.h>

#define TECPLOT
//#define VISIT

extern int Nx, Ny;

extern double xl, xr, yl, yr, Lx, Ly, Lambda;

extern double *x, *y, **rhsPDE, Operator;

void CreatePDEPlan(const int Nx, const int Ny, const double xl, const double xr, const double yl, const double yr);
void DestroyPDEPlan();

void SolvePDE(double **&U, double (*ConstructOperator)(const int, const int));

void ComputeError(double **U, double (*ExactFunction)(const double, const double));
void WriteFile(double **U, double (*ExactFunction)(const double, const double));
#endif
