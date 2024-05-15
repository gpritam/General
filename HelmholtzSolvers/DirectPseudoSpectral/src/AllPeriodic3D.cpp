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

#include "AllPeriodic3D.h"

int Nx, Ny, Nz;

double xl, xr, yl, yr, zl, zr, Lx, Ly, Lz, Lambda;

double *x, *y, *z, ***rhsPDE, Operator;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( const int Nx, 
                     const int Ny, 
                     const int Nz, 
                     const double xl, 
                     const double xr, 
                     const double yl, 
                     const double yr, 
                     const double zl, 
                     const double zr )
{
    CreateFFTPlan(Maximum(Nx,Ny,Nz));
    
    ::Nx = Nx;
    ::Ny = Ny;
    ::Nz = Nz;
    ::xl = xl;
    ::xr = xr;
    ::yl = yl;
    ::yr = yr;
    ::zl = zl;
    ::zr = zr;
    
    Lx = (xr-xl);
    Ly = (yr-yl);
    Lz = (zr-zl);
    
    Allocate(x,Nx);
    Allocate(y,Ny);
    Allocate(z,Nz);
    Allocate(rhsPDE,2*Nx,2*Ny,2*Nz);
    
    for (int i{}; i < Nx; ++i)
        x[i] = xl + i*(Lx/Nx);
    
    for (int j{}; j < Ny; ++j)
        y[j] = yl + j*(Ly/Ny);
    
    for (int k{}; k < Nz; ++k)
        z[k] = zl + k*(Lz/Nz);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    DestroyFFTPlan(Maximum(Nx,Ny,Nz));
    
    Deallocate(x,Nx);
    Deallocate(y,Ny);
    Deallocate(z,Nz);
    Deallocate(rhsPDE,2*Nx,2*Ny,2*Nz);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double InverseOperator ( const int i, 
                         const int j, 
                         const int k, 
                         double ***&U )
{
    U[i][j][k] = (Absolute(Operator) < 1.0E-15 ? 0.0 : rhsPDE[i][j][k]/Operator);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SolvePDE ( double ***&U, 
                double (*ConstructOperator)(const int, const int, const int) )
{
    FourierTransformXYZ(rhsPDE,1,Nx,Ny,Nz);
    
    int imax = Nx/2, jmax = Ny/2, kmax = Nz/2;
    
    for (int iw{}; iw <= imax; ++iw)
    {
        for (int jw{}; jw <= jmax; ++jw)
        {
            for (int kw{}; kw <= kmax; ++kw)
            {
                Operator = ConstructOperator(iw,jw,kw);
                
                if (iw == 0)
                {
                    if (jw == 0)
                    {
                        if (kw == 0)
                        {
                            InverseOperator(0,0,0,U);
                            InverseOperator(0,0,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(0,0,Nz,U);
                            InverseOperator(0,0,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(0,0,2*kw,U);
                            InverseOperator(0,0,2*kw+1,U);
                            InverseOperator(0,0,2*(Nz-kw),U);
                            InverseOperator(0,0,2*(Nz-kw)+1,U);
                        }
                    }
                    else if (jw == jmax)
                    {
                        if (kw == 0)
                        {
                            InverseOperator(0,jmax,0,U);
                            InverseOperator(0,jmax,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(0,jmax,Nz,U);
                            InverseOperator(0,jmax,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(0,jmax,2*kw,U);
                            InverseOperator(0,jmax,2*kw+1,U);
                            InverseOperator(0,jmax,2*(Nz-kw),U);
                            InverseOperator(0,jmax,2*(Nz-kw)+1,U);
                        }
                    }
                    else
                    {
                        if (kw == 0)
                        {
                            InverseOperator(0,jw,0,U);
                            InverseOperator(0,jw,1,U);
                            InverseOperator(0,Ny-jw,0,U);
                            InverseOperator(0,Ny-jw,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(0,jw,Nz,U);
                            InverseOperator(0,jw,Nz+1,U);
                            InverseOperator(0,Ny-jw,Nz,U);
                            InverseOperator(0,Ny-jw,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(0,jw,2*kw,U);
                            InverseOperator(0,jw,2*kw+1,U);
                            InverseOperator(0,Ny-jw,2*kw,U);
                            InverseOperator(0,Ny-jw,2*kw+1,U);
                            InverseOperator(0,jw,2*(Nz-kw),U);
                            InverseOperator(0,jw,2*(Nz-kw)+1,U);
                            InverseOperator(0,Ny-jw,2*(Nz-kw),U);
                            InverseOperator(0,Ny-jw,2*(Nz-kw)+1,U);
                        }
                    }
                }
                else if (iw == imax)
                {
                    if (jw == 0)
                    {
                        if (kw == 0)
                        {
                            InverseOperator(imax,0,0,U);
                            InverseOperator(imax,0,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(imax,0,Nz,U);
                            InverseOperator(imax,0,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(imax,0,2*kw,U);
                            InverseOperator(imax,0,2*kw+1,U);
                            InverseOperator(imax,0,2*(Nz-kw),U);
                            InverseOperator(imax,0,2*(Nz-kw)+1,U);
                        }
                    }
                    else if (jw == jmax)
                    {
                        if (kw == 0)
                        {
                            InverseOperator(imax,jmax,0,U);
                            InverseOperator(imax,jmax,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(imax,jmax,Nz,U);
                            InverseOperator(imax,jmax,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(imax,jmax,2*kw,U);
                            InverseOperator(imax,jmax,2*kw+1,U);
                            InverseOperator(imax,jmax,2*(Nz-kw),U);
                            InverseOperator(imax,jmax,2*(Nz-kw)+1,U);
                        }
                    }
                    else
                    {
                        if (kw == 0)
                        {
                            InverseOperator(imax,jw,0,U);
                            InverseOperator(imax,jw,1,U);
                            InverseOperator(imax,Ny-jw,0,U);
                            InverseOperator(imax,Ny-jw,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(imax,jw,Nz,U);
                            InverseOperator(imax,jw,Nz+1,U);
                            InverseOperator(imax,Ny-jw,Nz,U);
                            InverseOperator(imax,Ny-jw,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(imax,jw,2*kw,U);
                            InverseOperator(imax,jw,2*kw+1,U);
                            InverseOperator(imax,Ny-jw,2*kw,U);
                            InverseOperator(imax,Ny-jw,2*kw+1,U);
                            InverseOperator(imax,jw,2*(Nz-kw),U);
                            InverseOperator(imax,jw,2*(Nz-kw)+1,U);
                            InverseOperator(imax,Ny-jw,2*(Nz-kw),U);
                            InverseOperator(imax,Ny-jw,2*(Nz-kw)+1,U);
                        }
                    }
                }
                else
                {
                    if (jw == 0)
                    {
                        if (kw == 0)
                        {
                            InverseOperator(iw,0,0,U);
                            InverseOperator(iw,0,1,U);
                            InverseOperator(Nx-iw,0,0,U);
                            InverseOperator(Nx-iw,0,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(iw,0,Nz,U);
                            InverseOperator(iw,0,Nz+1,U);
                            InverseOperator(Nx-iw,0,Nz,U);
                            InverseOperator(Nx-iw,0,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(iw,0,2*kw,U);
                            InverseOperator(iw,0,2*kw+1,U);
                            InverseOperator(Nx-iw,0,2*kw,U);
                            InverseOperator(Nx-iw,0,2*kw+1,U);
                            InverseOperator(iw,0,2*(Nz-kw),U);
                            InverseOperator(iw,0,2*(Nz-kw)+1,U);
                            InverseOperator(Nx-iw,0,2*(Nz-kw),U);
                            InverseOperator(Nx-iw,0,2*(Nz-kw)+1,U);
                        }
                    }
                    else if (jw == jmax)
                    {
                        if (kw == 0)
                        {
                            InverseOperator(iw,jmax,0,U);
                            InverseOperator(iw,jmax,1,U);
                            InverseOperator(Nx-iw,jmax,0,U);
                            InverseOperator(Nx-iw,jmax,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(iw,jmax,Nz,U);
                            InverseOperator(iw,jmax,Nz+1,U);
                            InverseOperator(Nx-iw,jmax,Nz,U);
                            InverseOperator(Nx-iw,jmax,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(iw,jmax,2*kw,U);
                            InverseOperator(iw,jmax,2*kw+1,U);
                            InverseOperator(Nx-iw,jmax,2*kw,U);
                            InverseOperator(Nx-iw,jmax,2*kw+1,U);
                            InverseOperator(iw,jmax,2*(Nz-kw),U);
                            InverseOperator(iw,jmax,2*(Nz-kw)+1,U);
                            InverseOperator(Nx-iw,jmax,2*(Nz-kw),U);
                            InverseOperator(Nx-iw,jmax,2*(Nz-kw)+1,U);
                        }
                    }
                    else
                    {
                        if (kw == 0)
                        {
                            InverseOperator(iw,jw,0,U);
                            InverseOperator(iw,jw,1,U);
                            InverseOperator(Nx-iw,jw,0,U);
                            InverseOperator(Nx-iw,jw,1,U);
                            InverseOperator(iw,Ny-jw,0,U);
                            InverseOperator(iw,Ny-jw,1,U);
                            InverseOperator(Nx-iw,Ny-jw,0,U);
                            InverseOperator(Nx-iw,Ny-jw,1,U);
                        }
                        else if (kw == kmax)
                        {
                            InverseOperator(iw,jw,Nz,U);
                            InverseOperator(iw,jw,Nz+1,U);
                            InverseOperator(Nx-iw,jw,Nz,U);
                            InverseOperator(Nx-iw,jw,Nz+1,U);
                            InverseOperator(iw,Ny-jw,Nz,U);
                            InverseOperator(iw,Ny-jw,Nz+1,U);
                            InverseOperator(Nx-iw,Ny-jw,Nz,U);
                            InverseOperator(Nx-iw,Ny-jw,Nz+1,U);
                        }
                        else
                        {
                            InverseOperator(iw,jw,2*kw,U);
                            InverseOperator(iw,jw,2*kw+1,U);
                            InverseOperator(Nx-iw,jw,2*kw,U);
                            InverseOperator(Nx-iw,jw,2*kw+1,U);
                            InverseOperator(iw,Ny-jw,2*kw,U);
                            InverseOperator(iw,Ny-jw,2*kw+1,U);
                            InverseOperator(Nx-iw,Ny-jw,2*kw,U);
                            InverseOperator(Nx-iw,Ny-jw,2*kw+1,U);
                            InverseOperator(iw,jw,2*(Nz-kw),U);
                            InverseOperator(iw,jw,2*(Nz-kw)+1,U);
                            InverseOperator(Nx-iw,jw,2*(Nz-kw),U);
                            InverseOperator(Nx-iw,jw,2*(Nz-kw)+1,U);
                            InverseOperator(iw,Ny-jw,2*(Nz-kw),U);
                            InverseOperator(iw,Ny-jw,2*(Nz-kw)+1,U);
                            InverseOperator(Nx-iw,Ny-jw,2*(Nz-kw),U);
                            InverseOperator(Nx-iw,Ny-jw,2*(Nz-kw)+1,U);
                        }
                    }
                }
            }
        }
    }
    
    FourierTransformXYZ(U,-1,Nx,Ny,Nz);
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( double ***U, 
                    double (*ExactFunction)(const double, const double, const double) )
{
    double Error = 0.0;
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            for (int k{}; k < Nz; ++k)
                Error = Maximum(Error,Absolute(U[i][j][k]-ExactFunction(x[i],y[j],z[k])));
    
    std::cout << "Error = " << Error << std::endl << std::endl;
}

//_______________________________________________________________________________
// Write file in Tecplot format
//_______________________________________________________________________________
void WriteFile ( double ***U, 
                 double (*ExactFunction)(const double, const double, const double) )
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
    FileWrite << "Variables = \"X\",\"Y\",\"Z\",\"U\",\"Exact\"" << std::endl;
    FileWrite << "Zone I = " << Nx << ", J = " << Ny << ", K = " << Nz << ", DATAPACKING=POINT" << std::endl;
    
    for (int k{}; k < Nz; ++k)
        for (int j{}; j < Ny; ++j)
            for (int i{}; i < Nx; ++i)
                FileWrite << x[i] << "\t" << y[j] << "\t" << z[k] << "\t" << U[i][j][k] << "\t" << ExactFunction(x[i],y[j],z[k]) << std::endl;
    
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
    FileWrite << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << std::endl; 
    FileWrite << "POINTS " << Nx*Ny*Nz << " FLOAT" << std::endl;
    
    for (int k{}; k < Nz; ++k)
        for (int j{}; j < Ny; ++j)
            for (int i{}; i < Nx; ++i)
                FileWrite << x[i] << "\t" << y[j] << "\t" << z[k] << std::endl;
    
    FileWrite << std::endl << "POINT_DATA " << Nx*Ny*Nz << std::endl;
    FileWrite << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
    
    for (int k{}; k < Nz; ++k)
        for (int j{}; j < Ny; ++j)
            for (int i{}; i < Nx; ++i)
                FileWrite << U[i][j][k] << std::endl;
    
    FileWrite.close();
    #endif
    
    Deallocate(s,200);
}
