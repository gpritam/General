//_______________________________________________________________________________
//_______________________________________________________________________________
// Program for 1D Helmholtz equation with periodic boundary condition.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 15.03.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

#include "AllPeriodic1D.h"

int Nx;

double xl, xr, Lambda, Lx;

double *x, *rhsPDE, Operator;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( const int Nx, 
                     const double xl, 
                     const double xr )
{
    CreateFFTPlan(Nx);
    
    ::Nx = Nx;
    ::xl = xl;
    ::xr = xr;
    
    Lx = (xr-xl);
    
    Allocate(x,Nx);
    Allocate(rhsPDE,2*Nx);
    
    for (int i{}; i < Nx; ++i)
        x[i] = xl + i*(Lx/Nx);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    DestroyFFTPlan(Nx);
    
    Deallocate(x,Nx);
    Deallocate(rhsPDE,2*Nx);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double InverseOperator ( const int i, 
                         double *&U )
{
    U[i] = (Absolute(Operator) < 1.0E-15 ? 0.0 : rhsPDE[i]/Operator);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SolvePDE ( double *&U, 
                double (*ConstructOperator)(const int) )
{
    FourierTransformX(rhsPDE,1,Nx);
    
    int imax = Nx/2;
    
    for (int iw{}; iw <= imax; ++iw)
    {
        Operator = ConstructOperator(iw);
        
        if (iw == 0)
        {
            InverseOperator(0,U);
            InverseOperator(1,U);
        }
        else if (iw == imax)
        {
            InverseOperator(Nx,U);
            InverseOperator(Nx+1,U);
        }
        else
        {
            InverseOperator(2*iw,U);
            InverseOperator(2*iw+1,U);
            InverseOperator(2*(Nx-iw),U);
            InverseOperator(2*(Nx-iw)+1,U);
        }
    }
    
    FourierTransformX(U,-1,Nx);
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( double *U, 
                    double (*ExactFunction)(const double) )
{
    double Error = 0.0;
    
    for (int i{}; i < Nx; ++i)
        Error = Maximum(Error,Absolute(U[i]-ExactFunction(x[i])));
    
    std::cout << "Error = " << Error << std::endl << std::endl;
}

//_______________________________________________________________________________
// Write file in Tecplot format
//_______________________________________________________________________________
void WriteFile ( double *U, 
                 double (*ExactFunction)(const double) )
{
    char *s;
    
    Allocate(s,200);
    
    #ifdef TECPLOT
    sprintf(s,"Output/Solution.tec");
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    for (int i{}; i < Nx; ++i)
        FileWrite << x[i] << "\t" << U[i] << "\t" << ExactFunction(x[i]) << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    sprintf(s,"Output/Solution.curve");
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    FileWrite << "#Computed" << std::endl;
    
    for (int i{}; i < Nx; ++i)
        FileWrite << x[i] << "\t" << U[i] << std::endl;
    
    FileWrite.close();
    #endif
    
    Deallocate(s,200);
}
