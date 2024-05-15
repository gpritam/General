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
// g++ -O3 -o run A04_Quick_Sort.cpp
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
		A[i] = ((double)rand()/((double)(RAND_MAX)))*100.0;
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
template <class T> int  Partition ( T *A, 
									const int l, 
									const int r )
{
	int j = (l-1);
	
	T pivot = A[r];
	
	for (int i{l}; i < r; ++i)
		if (A[i] <= pivot)
			Swap(A[++j],A[i]);
	
	Swap(A[++j],A[r]);
	
	return j;
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
template <class T> void QuickSort ( T *A, 
									const int l, 
									const int r )
{
	if (l < r)
	{
		int m = Partition(A,l,r);
		
		QuickSort(A,l,m-1);
		QuickSort(A,m+1,r);
	}
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
template <class T> void QuickSort ( T *A, 
									const int n )
{
	QuickSort(A,0,n-1);
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
	
	QuickSort(A,n);
	
	for (int i{}; i < n; ++i)
		std::cout << A[i] << std::endl;
	
	Deallocate(A,n);
	
	return 0;
}
