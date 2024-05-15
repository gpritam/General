//_______________________________________________________________________________
//_______________________________________________________________________________
// Program for plotting Mandlebrot set.
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
// make TARGET=MandelbrotSet.cpp
// make run
//_______________________________________________________________________________

#include <bits/stdc++.h>

using namespace std::chrono;

#define TECPLOT
//#define VISIT

int main()
{
    std::cout.flags( std::ios::dec );
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    const double xmin = -2.5, xmax = 1.0, ymin = -1.0, ymax = 1.0;
    
    const int Imax = 1000, Nx = 1024, Ny = Nx*(ymax-ymin)/(xmax-xmin);
    
    double dx = (xmax-xmin)/Nx, dy = (ymax-ymin)/Ny, cx, cy, zx, zy, Zx, Zy;
    
    int k;
    
    // Plot Mandelbrot set
    #ifdef TECPLOT
    std::ofstream FileWrite("Output/Mandelbrot.tec", std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::scientific);
    FileWrite.precision(8);
    
    if (!FileWrite)
    {
        std::cout << "Output file couldnot be opened." << std::endl;
        exit(1);
    }
    
    FileWrite << "TITLE = \"Mandelbrot set\"" << std::endl;
    FileWrite << "VARIABLES = \"x\", \"y\", \"fractal\""<< std::endl;
    FileWrite << "Zone I = " << Ny << ", J = " << Nx << ", DATAPACKING=POINT" << std::endl;
    
    for (int i{}; i < Nx; ++i)
    {
        cx = xmin + i*dx;
        
        for (int j{}; j < Ny; ++j)
        {
            cy = ymax - j*dy;
            
            zx = 0.0;
            zy = 0.0;
            
            for (k = 1; (k < Imax) && ((zx*zx + zy*zy) < 4.0); ++k)
            {
                Zx = zx*zx - zy*zy + cx;
                Zy = 2.0*zx*zy + cy;
                
                zx = Zx;
                zy = Zy;
            }
            
            FileWrite << cx << "\t" << cy << "\t" << double(k)/Imax << std::endl;
        }
    }
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    std::ofstream FileWriteVisIt("Output/Mandelbrot.vtk", std::ios::out);
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
    FileWriteVisIt << "SPACING " << (xmax-xmin)/Nx << " " << (ymax-ymin)/Ny << " 1" << std::endl << std::endl;
    
    FileWriteVisIt << "POINT_DATA " << Nx*Ny << std::endl;
    FileWriteVisIt << "SCALARS fractal float" << std::endl;
    FileWriteVisIt << "LOOKUP_TABLE default" << std::endl;
    
    for (int j{}; j < Ny; ++j)
    {
        cy = ymax - j*dy;
        
        for (int i{}; i < Nx; ++i)
        {
            cx = xmin + i*dx;
            
            zx = 0.0;
            zy = 0.0;
            
            for (k = 1; (k < Imax) && ((zx*zx + zy*zy) < 4.0); ++k)
            {
                Zx = zx*zx - zy*zy + cx;
                Zy = 2.0*zx*zy + cy;
                
                zx = Zx;
                zy = Zy;
            }
            
            FileWriteVisIt << double(k)/Imax << std::endl;
        }
    }
    
    FileWriteVisIt.close();
    #endif
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    
    return 0;
}
