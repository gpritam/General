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

#include "AllPeriodic2D.h"

int Nx, Ny;

double xl, xr, yl, yr, Lx, Ly, Lambda;

double *x, *y, **rhsPDE, Operator;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( const int Nx, 
                     const int Ny, 
                     const double xl, 
                     const double xr, 
                     const double yl, 
                     const double yr )
{
    CreateFFTPlan(Maximum(Nx,Ny));
    
    ::Nx = Nx;
    ::Ny = Ny;
    ::xl = xl;
    ::xr = xr;
    ::yl = yl;
    ::yr = yr;
    
    Lx = (xr-xl);
    Ly = (yr-yl);
    
    Allocate(x,Nx);
    Allocate(y,Ny);
    Allocate(rhsPDE,2*Nx,2*Ny);
    
    for (int i{}; i < Nx; ++i)
        x[i] = xl + i*(Lx/Nx);
    
    for (int j{}; j < Ny; ++j)
        y[j] = yl + j*(Ly/Ny);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    DestroyFFTPlan(Maximum(Nx,Ny));
    
    Deallocate(x,Nx);
    Deallocate(y,Ny);
    Deallocate(rhsPDE,2*Nx,2*Ny);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double InverseOperator ( const int i, 
                         const int j, 
                         double **&U )
{
    U[i][j] = (Absolute(Operator) < 1.0E-15 ? 0.0 : rhsPDE[i][j]/Operator);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SolvePDE ( double **&U, 
                double (*ConstructOperator)(const int, const int) )
{
    FourierTransformXY(rhsPDE,1,Nx,Ny);
    
    int imax = Nx/2, jmax = Ny/2;
    
    for (int iw{}; iw <= imax; ++iw)
    {
        for (int jw{}; jw <= jmax; ++jw)
        {
            Operator = ConstructOperator(iw,jw);
            
            if (iw == 0)
            {
                if (jw == 0)
                {
                    InverseOperator(0,0,U);
                    InverseOperator(0,1,U);
                }
                else if (jw == jmax)
                {
                    InverseOperator(0,Ny,U);
                    InverseOperator(0,Ny+1,U);
                }
                else
                {
                    InverseOperator(0,2*jw,U);
                    InverseOperator(0,2*jw+1,U);
                    InverseOperator(0,2*(Ny-jw),U);
                    InverseOperator(0,2*(Ny-jw)+1,U);
                }
            }
            else if (iw == imax)
            {
                if (jw == 0)
                {
                    InverseOperator(imax,0,U);
                    InverseOperator(imax,1,U);
                }
                else if (jw == jmax)
                {
                    InverseOperator(imax,Ny,U);
                    InverseOperator(imax,Ny+1,U);
                }
                else
                {
                    InverseOperator(imax,2*jw,U);
                    InverseOperator(imax,2*jw+1,U);
                    InverseOperator(imax,2*(Ny-jw),U);
                    InverseOperator(imax,2*(Ny-jw)+1,U);
                }
            }
            else
            {
                if (jw == 0)
                {
                    InverseOperator(iw,0,U);
                    InverseOperator(iw,1,U);
                    InverseOperator((Nx-iw),0,U);
                    InverseOperator((Nx-iw),1,U);
                }
                else if (jw == jmax)
                {
                    InverseOperator(iw,Ny,U);
                    InverseOperator(iw,Ny+1,U);
                    InverseOperator((Nx-iw),Ny,U);
                    InverseOperator((Nx-iw),Ny+1,U);
                }
                else
                {
                    InverseOperator(iw,2*jw,U);
                    InverseOperator(iw,2*jw+1,U);
                    InverseOperator(iw,2*(Ny-jw),U);
                    InverseOperator(iw,2*(Ny-jw)+1,U);
                    InverseOperator((Nx-iw),2*jw,U);
                    InverseOperator((Nx-iw),2*jw+1,U);
                    InverseOperator((Nx-iw),2*(Ny-jw),U);
                    InverseOperator((Nx-iw),2*(Ny-jw)+1,U);
                }
            }
        }
    }
    
    FourierTransformXY(U,-1,Nx,Ny);
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( double **U, 
                    double (*ExactFunction)(const double, const double) )
{
    double Error = 0.0;
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            Error = Maximum(Error,Absolute(U[i][j]-ExactFunction(x[i],y[j])));
    
    std::cout << "Error = " << Error << std::endl << std::endl;
}

//_______________________________________________________________________________
// Write file in Tecplot format
//_______________________________________________________________________________
void WriteFile ( double **U, 
                 double (*ExactFunction)(const double, const double) )
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
    
    FileWrite << "TITLE = \"Helmholtz solver\"" << std::endl;
    FileWrite << "Variables = \"X\",\"Y\",\"U\",\"Exact\"" << std::endl;
    FileWrite << "Zone I = " << Nx << ", J = " << Ny << ", DATAPACKING=POINT" << std::endl;
    
    for (int j{}; j < Ny; ++j)
        for (int i{}; i < Nx; ++i)
            FileWrite << x[i] << "\t" << y[j] << "\t" << U[i][j] << "\t" << ExactFunction(x[i],y[j]) << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    sprintf(s,"Output/Solution.vtk");
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    FileWrite << "# vtk DataFile Version 3.1" << std::endl;
    FileWrite << "Helmholtz solver" << std::endl;
    FileWrite << "ASCII" << std::endl;
    FileWrite << "DATASET STRUCTURED_GRID" << std::endl;
    FileWrite << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
    FileWrite << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
    
    for (int j{}; j < Ny; ++j)
        for (int i{}; i < Nx; ++i)
            FileWrite << x[i] << "\t" << y[j] << "\t" << 0.0 << std::endl;
    
    FileWrite << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
    FileWrite << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
    
    for (int j{}; j < Ny; ++j)
        for (int i{}; i < Nx; ++i)
            FileWrite << U[i][j] << std::endl;
    
    FileWrite.close();
    #endif
    
    Deallocate(s,200);
}
