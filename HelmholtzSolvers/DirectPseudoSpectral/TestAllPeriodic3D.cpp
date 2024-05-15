//_______________________________________________________________________________
//_______________________________________________________________________________
// Program for 3D Helmholtz equation with periodic boundary condition.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 20.03.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

//_______________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=TestAllPeriodic3D.cpp
// make run
//_______________________________________________________________________________

#include "AllPeriodic3D.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x, 
                       const double y, 
                       const double z )
{
    return sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x, 
             const double y, 
             const double z )
{
    return -(64.0*PI*PI/(Lx*Lx) + 36.0*PI*PI/(Ly*Ly) + 16.0*PI*PI/(Lz*Lz) + Lambda)*sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
double ConstructOperator ( const int iw, 
                           const int jw, 
                           const int kw )
{
    return (-4.0*PI*PI*(iw*iw/(Lx*Lx) + jw*jw/(Ly*Ly) + kw*kw/(Lz*Lz)) - Lambda);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main ()
{
    std::cout.flags(std::ios::dec);
    std::cout.precision(8);    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // Create plan
    CreatePDEPlan(64,32,256,0.0,1.0,0.0,1.0,0.0,1.0);
    
    // Solve first PDE
    double ***V;
    
    Allocate(V,2*Nx,2*Ny,2*Nz);
    
    auto start = high_resolution_clock::now();
    
    Lambda = 1.0;
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            for (int k{}; k < Nz; ++k)
                    rhsPDE[i][j][k] = RHS(x[i],y[j],z[k]);
    
    SolvePDE(V,ConstructOperator);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " milliseconds." << std::endl << std::endl;
    
    WriteFile(V,ExactFunction);
    
    ComputeError(V,ExactFunction);
    
    // Deatroy plan    
    Deallocate(V,2*Nx,2*Ny,2*Nz);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    return 0;
}
