#pragma once

#ifndef AllPeriodic3D_H
#define AllPeriodic3D_H

#include "General.h"
#include "Fourier_Quadrature.h"
#include <bits/stdc++.h>

#define TECPLOT
//#define VISIT

extern int Nx, Ny, Nz;

extern double xl, xr, yl, yr, zl, zr, Lx, Ly, Lz, Lambda;

extern double *x, *y, *z, ***rhsPDE, Operator;

void CreatePDEPlan(const int Nx, const int Ny, const int Nz, const double xl, const double xr, const double yl, const double yr, const double zl, const double zr);
void DestroyPDEPlan();

void SolvePDE(double ***&U, double (*ConstructOperator)(const int, const int, const int));

void ComputeError(double ***U, double (*ExactFunction)(const double, const double, const double));
void WriteFile(double ***U, double (*ExactFunction)(const double, const double, const double));
#endif
