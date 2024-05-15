//_______________________________________________________________________________
//_______________________________________________________________________________
// Developed by: Dr. Pritam Giri
// Date : 1.09.2022
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

//_______________________________________________________________________________
// To run this code, issue the following commands
// 
// g++ -O3 -o run A00_Bubble_Sort.cpp
// ./run
//_______________________________________________________________________________

#include <bits/stdc++.h>

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
template <class T> void Swap ( T &a, T &b )
{
	if (a != b)
	{
		a = a + b;
		b = a - b;
		a = a - b;
	}
}

//_______________________________________________________________________________
// Print an error message and exit the program
//_______________________________________________________________________________
void ErrorMessage ( const char* message )
{
	std::cout << message << std::endl;
	
	exit(1);
}

//_______________________________________________________________________________
// Allocate an 1D array in contiguous memory locations
// First element can be accessed by A[-offset]
//_______________________________________________________________________________
template <class T> void Allocate ( T *&A,
								   const int n,
								   const int offset = 0 )
{
	A = new (std::nothrow) T [n];
	
	if (A == 0)
		ErrorMessage("Error: Memory can not be allocated!");
	
	for (int i{}; i < n; i++)
		A[i] = T(0.0);
	
	A += offset;
}

//_______________________________________________________________________________
// Deallocate an 1D array
//_______________________________________________________________________________
template <class T> void Deallocate ( T *A,
									 const int n,
									 const int offset = 0 )
{
	A -= offset;
	
	delete [] A;
	
	A = nullptr;
}

//_______________________________________________________________________________
// Fill a vector with elements between 0 and 100
//_______________________________________________________________________________
template <class T> void FillArray ( T *&A, 
				  					const int n )
{
	for (int i{}; i < n; ++i)
		A[i] = 100.0*double(rand())/double(RAND_MAX);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
template <class T> void BubbleSort ( T *A, 
									 const int n )
{
	bool sorted;
	
	for (int i{}; i < (n-1); ++i)
	{
		sorted = true;
		
		for (int j{}; j < (n-i-1); ++j)
		{
			if (A[j] > A[j+1])
			{
				Swap(A[j],A[j+1]);
				
				sorted = false;
			}
		}
		
		if (sorted)
			break;
	}
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main(int argc, char *argv[])
{
	std::cout.flags(std::ios::fixed);
	std::cout.precision(4);
	
	int *A;
	const int n = 10;
	
	Allocate(A,n);
	FillArray(A,n);
	
	BubbleSort(A,n);
	
	for (int i{}; i < n; ++i)
		std::cout << A[i] << std::endl;
	
	Deallocate(A,n);
	
	return 0;
}
