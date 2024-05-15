//_______________________________________________________________________________
//_______________________________________________________________________________
// Test program for 1D Helmholtz equation with periodic boundary condition.
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
// make TARGET=TestAllPeriodic1D.cpp
// make run
//_______________________________________________________________________________

#include "AllPeriodic1D.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x )
{
    return sin(8.0*PI*(x-xl)/Lx);
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x )
{
    return -(64.0*PI*PI/(Lx*Lx) + Lambda)*sin(8.0*PI*(x-xl)/Lx);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
double ConstructOperator ( const int iw )
{
    return (-4.0*PI*PI*iw*iw/(Lx*Lx) - Lambda);
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
    CreatePDEPlan(512,0.0,1.0);
    
    // Solve first PDE
    double *V;
    
    Allocate(V,2*Nx);
    
    auto start = high_resolution_clock::now();
    
    Lambda = 1.0;
    
    for (int i{}; i < Nx; ++i)
        rhsPDE[i] = RHS(x[i]);
    
    SolvePDE(V,ConstructOperator);
    
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " microseconds." << std::endl << std::endl;
    
    WriteFile(V,ExactFunction);
    
    ComputeError(V,ExactFunction);
    
    // Deatroy plan    
    Deallocate(V,2*Nx);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    return 0;
}
