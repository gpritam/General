//_______________________________________________________________________________
//_______________________________________________________________________________
// Program for 2D Helmholtz equation with periodic boundary condition.
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
// make TARGET=TestAllPeriodic2D.cpp
// make run
//_______________________________________________________________________________

#include "AllPeriodic2D.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x, 
                       const double y )
{
    return sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x, 
             const double y )
{
    return -(64.0*PI*PI/(Lx*Lx) + 36.0*PI*PI/(Ly*Ly) + Lambda)*sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
double ConstructOperator ( const int iw, 
                           const int jw )
{
    return (-4.0*PI*PI*(iw*iw/(Lx*Lx) + jw*jw/(Ly*Ly)) - Lambda);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main ()
{
    std::cout.flags( std::ios::dec );
    std::cout.precision(8);    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // Create plan
    CreatePDEPlan(512,512,0.0,1.0,0.0,1.0);
    
    // Solve first PDE
    double **V;
    
    Allocate(V,2*Nx,2*Ny);
    
    auto start = high_resolution_clock::now();
    
    Lambda = 1.0;
    
    for (int j{}; j < Ny; ++j)
        for (int i{}; i < Nx; ++i)
            rhsPDE[i][j] = RHS(x[i],y[j]);
    
    SolvePDE(V,ConstructOperator);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " milliseconds." << std::endl << std::endl;
    
    WriteFile(V,ExactFunction);
    
    ComputeError(V,ExactFunction);
    
    // Deatroy plan    
    Deallocate(V,2*Nx,2*Ny);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    return 0;
}
