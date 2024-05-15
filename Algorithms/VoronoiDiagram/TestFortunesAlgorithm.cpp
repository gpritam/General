//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Serial implementation of Fortune's sweep algorithm
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 27.09.2022
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=TestFortunesAlgorithm.cpp
// make run
//________________________________________________________________________________________________

#include "FortunesAlgorithm.h"

using namespace std::chrono;

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main()
{
    std::cout.flags(std::ios::dec);
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	const int N = 380;
	
	const double xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0, CircleEventAccuracy = 1.0E-8;
	
	double **Vertices;
    
	Allocate(Vertices,N,2);
	
	GeneratePoints(Vertices,xmin,xmax,ymin,ymax,N);
	
	GenerateVoronoiDiagram(Vertices,xmin,xmax,ymin,ymax,N,CircleEventAccuracy);
	
	Deallocate(Vertices,N,2);
	
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    std::cout << std::endl << "Total time elapsed : " << duration.count() << " milliseconds." << std::endl << std::endl;
    
    return 0;
}
