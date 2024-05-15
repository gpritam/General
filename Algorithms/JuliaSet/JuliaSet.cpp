//_______________________________________________________________________________
//_______________________________________________________________________________
// Program for plotting Julia sets when c is lying on the periphery of the unit 
// circle having center z=0 on the argand plane.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 15.08.2021
// IIT Kanpur
//_______________________________________________________________________________
//_______________________________________________________________________________

//_______________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=JuliaSet.cpp
// make run
//_______________________________________________________________________________

#include <bits/stdc++.h>

using namespace std::chrono;

#define TECPLOT
//#define VISIT

const double PI = 4.0*atan(1.0);

//_______________________________________________________________________________
// Determines if a point in the argand plane is included in the Julia set or, not 
// for a fixed c.
// 
// z_n = z_(n-1)^2 + c
// 
// z_(n-1) = zx + izy and, c = cx + icy
//_______________________________________________________________________________
double IsInJuliaSet ( double zx,
                      double zy,
                      const double cx,
                      const double cy )
{
    const int Imax = 1000;
    
    double Zx, Zy;
    
    int k = 0;
    
    while (((zx*zx + zy*zy) < 100.0) && (k < Imax))
    {
        Zx = zx*zx - zy*zy + cx;
        Zy = 2.0*zx*zy + cy;
        
        zx = Zx;
        zy = Zy;
        
        k++;
    }
    
    return double(k)/Imax;
} 

//_______________________________________________________________________________
// This is the main code
//_______________________________________________________________________________
int main() 
{
    std::cout.flags( std::ios::dec );
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    const double xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5;
    
    const int Nx = 512, Ny = Nx*(ymax-ymin)/(xmax-xmin), Nf = 100;
    
    double zx, zy, cx, cy, alpha;
    
    char *s = new char[200];
    
    // Plot Julia sets
    #ifdef TECPLOT
    for (int k{}; k < Nf; ++k)
    {
        alpha = 2.0*k*PI/Nf;
        
        cx = 0.7885*cos(alpha); cy = 0.7885*sin(alpha);
        
        sprintf(s,"Output/JuliaSet-%d.tec",k);
        
        std::ofstream FileWrite(s, std::ios::out);
        FileWrite.flags( std::ios::dec | std::ios::scientific);
        FileWrite.precision(8);
        
        if (!FileWrite)
        {
            std::cout << "Output file couldnot be opened." << std::endl;
            exit(1);
        }
        
        FileWrite << "TITLE = \"Julia set\"" << std::endl;
        FileWrite << "VARIABLES = \"x\", \"y\", \"fractal\""<< std::endl;
        FileWrite << "Zone I = " << Ny << ", J = " << Nx << ", DATAPACKING=POINT" << std::endl;
        
        for (int i{}; i < Nx; ++i)
        {
            zx = xmin + i*(xmax-xmin)/(Nx-1.0);
            
            for (int j = 0; j < Ny; ++j)
            {
                zy = ymin + j*(ymax-ymin)/(Ny-1.0);
                
                FileWrite << zx << "\t" << zy << "\t" << IsInJuliaSet(zx,zy,cx,cy) << std::endl;
            }
        }
        
        FileWrite.close();
        
        std::cout << "Julia set: " << k << std::endl;
    }
    #endif
    
    #ifdef VISIT
    for (int k{}; k < Nf; ++k)
    {
        alpha = 2.0*k*PI/Nf;
        
        cx = 0.7885*cos(alpha); cy = 0.7885*sin(alpha);
        
        sprintf(s,"Output/JuliaSet-%d.vtk",k);
        
        std::ofstream FileWriteVisIt(s, std::ios::out);
        FileWriteVisIt.flags( std::ios::dec | std::ios::scientific);
        FileWriteVisIt.precision(8);
        
        if (!FileWriteVisIt)
        {
            std::cout << "Output file couldnot be opened." << std::endl;
            exit(1);
        }
        
        FileWriteVisIt << "# vtk DataFile Version 3.1" << std::endl;
        FileWriteVisIt << "Julia set" << std::endl;
        FileWriteVisIt << "ASCII" << std::endl;
        FileWriteVisIt << "DATASET STRUCTURED_POINTS" << std::endl;
        FileWriteVisIt << "DIMENSIONS " << Nx << " " << Ny << " 1" << std::endl;
        FileWriteVisIt << "ORIGIN " << xmin << " " << ymin << " 0" << std::endl;
        FileWriteVisIt << "SPACING " << (xmax-xmin)/(Nx-1.0) << " " << (ymax-ymin)/(Ny-1.0) << " 1" << std::endl << std::endl;
        
        FileWriteVisIt << "POINT_DATA " << Nx*Ny << std::endl;
        FileWriteVisIt << "SCALARS fractal float" << std::endl;
        FileWriteVisIt << "LOOKUP_TABLE default" << std::endl;
        
        for (int j{}; j < Ny; ++j)
        {
            zy = ymin + j*(ymax-ymin)/(Ny-1.0);
            
            for (int i{}; i < Nx; ++i)
            {
                zx = xmin + i*(xmax-xmin)/(Nx-1.0);
                
                FileWriteVisIt << IsInJuliaSet(zx,zy,cx,cy) << std::endl;
            }
        }
        
        FileWriteVisIt.close();
        
        std::cout << "Julia set: " << k << std::endl;
    }
    #endif
    
    delete [] s;
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
