#include "SpaceFillingCurves.h"

//_______________________________________________________________________________
// Operator precedence
// 
// 1. ()  	2. ~  	3. * / %  	4. + -  	5. << >>  	6. &  	7. ^  	8. | 
//_______________________________________________________________________________

//_______________________________________________________________________________
// Count number of bits
//_______________________________________________________________________________
Uint CountBits ( Int x )
{
	Uint Nb = 0;
	
	while (x)
	{
		Nb++;
		
		x >>= 1;
	}
	
	return Nb;
}

//_______________________________________________________________________________
// x ^= (x/2)
//_______________________________________________________________________________
void BinaryToGray ( Int &x )
{
	x ^= (x >> 1);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void GrayToBinary ( Int &x )
{
	Int y = x;
	
	while (y)
	{
		y >>= 1;
		x ^= y;
	}
}

//_______________________________________________________________________________
// Swap two unsigned integers
//_______________________________________________________________________________
void SwapIntegers ( Int &x, 
					Int &y )
{
	x ^= y;
	y ^= x;
	x ^= y;
}

//_______________________________________________________________________________
// Get bit of a particular position (returns true if the bit is 1)
//_______________________________________________________________________________
bool GetBitAtPosition ( Int x,
						Uint n )
{
	if (n >= sizeof(x)*CHAR_BIT)
	{
		std::cout << "Enter a valid bit position!" << std::endl;
		
		exit(1);
	}
	
	Int M = 1 << n;
	
	return ( (x & M) ? true : false );
}

//_______________________________________________________________________________
// Set bit at a particular position
//_______________________________________________________________________________
void SetBitAtPosition ( Int &x,
						Uint n )
{
	if (n >= sizeof(x)*CHAR_BIT)
	{
		std::cout << "Enter a valid bit position!" << std::endl;
		
		exit(1);
	}
	
	Int M = 1 << n;
	
	x |= M;
}

//_______________________________________________________________________________
// Clear bit at a particular position
//_______________________________________________________________________________
void ClearBitAtPosition ( Int &x,
						  Uint n )
{
	if (n >= sizeof(x)*CHAR_BIT)
	{
		std::cout << "Enter a valid bit position!" << std::endl;
		
		exit(1);
	}
	
	Int M = 1 << n;
	
	x &= ~M;
}

//_______________________________________________________________________________
// Toggle bit at a particular position
//_______________________________________________________________________________
void ToggleBitAtPosition ( Int &x,
						   Uint n )
{
	if (n >= sizeof(x)*CHAR_BIT)
	{
		std::cout << "Enter a valid bit position!" << std::endl;
		
		exit(1);
	}
	
	Int M = 1 << n;
	
	x ^= M;
}

//_______________________________________________________________________________
// This program prints the binary representation of x
//_______________________________________________________________________________
void PrintBits ( Int x )
{
	//Uint Nb = sizeof(x)*CHAR_BIT;
	
	Uint Nb = CountBits(x);
	
	for (Uint i = 0; i < Nb; ++i)
	{
		Int M = 1 << (Nb-i-1);
		
		std::cout << ( (x & M) ? 1 : 0 );
	}
	
	std::cout << std::endl;
}

//_______________________________________________________________________________
// Interleave bits of 2 integers i and j: (i,j) -> I
// Nb : It can be (Number of bits of I)/2 atmost
//_______________________________________________________________________________
void Mortonize ( Int i, 
				 Int j,  
				 Int &I, 
				 Uint Nb )
{
	if (!Nb)
	{
		Nb = sizeof(I)*CHAR_BIT >> 1;
	}
	else if ((sizeof(I)*CHAR_BIT >> 1) < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving bits!" << std::endl;
		
		exit(1);
	}
	
	I = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (2*n);
		
		if ( j & M0 )
			I |= M1;
		
		M1 <<= 1;
		
		if ( i & M0 )
			I |= M1;
	}
}

//_______________________________________________________________________________
// Interleave bits of 3 integers i, j, and k: (i,j,k) -> I
// Nb : It can be floor((Number of bits of I)/3) atmost
//_______________________________________________________________________________
void Mortonize ( Int i, 
				 Int j,  
				 Int k, 
				 Int &I, 
				 Uint Nb )
{
	if (!Nb)
	{
		Nb = (sizeof(I)*CHAR_BIT)/3;
	}
	else if ((sizeof(I)*CHAR_BIT)/3 < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving bits!" << std::endl;
		
		exit(1);
	}
	
	I = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (3*n);
		
		if ( k & M0 )
			I |= M1;
		
		M1 <<= 1;
		
		if ( j & M0 )
			I |= M1;
		
		M1 <<= 1;
		
		if ( i & M0 )
			I |= M1;
	}
}

//_______________________________________________________________________________
// De-interleave bits for 2 integers i and j: I -> (i,j)
// Nb : It can be (Number of bits of I)/2 atmost
//_______________________________________________________________________________
void DeMortonize ( Int &i, 
				   Int &j,  
				   Int I, 
				   Uint Nb )
{
	if (!Nb)
	{
		Nb = sizeof(I)*CHAR_BIT >> 1;
	}
	else if ((sizeof(I)*CHAR_BIT >> 1) < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving Nb bits!" << std::endl;
		
		exit(1);
	}
	
	i = j = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (2*n);
		
		if ( I & M1 )
			j |= M0;
		
		M1 <<= 1;
		
		if ( I & M1 )
			i |= M0;
	}
}

//_______________________________________________________________________________
// De-interleave bits for 3 integers i, j, and k: I -> (i,j,k)
// Nb : It can be floor((Number of bits of I)/3) atmost
//_______________________________________________________________________________
void DeMortonize ( Int &i, 
				   Int &j, 
				   Int &k,  
				   Int I, 
				   Uint Nb )
{
	if (!Nb)
	{
		Nb = (sizeof(I)*CHAR_BIT)/3;
	}
	else if ((sizeof(I)*CHAR_BIT)/3 < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving Nb bits!" << std::endl;
		
		exit(1);
	}
	
	i = j = k = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (3*n);
		
		if ( I & M1 )
			k |= M0;
		
		M1 <<= 1;
		
		if ( I & M1 )
			j |= M0;
		
		M1 <<= 1;
		
		if ( I & M1 )
			i |= M0;
	}
}

//_______________________________________________________________________________
// (i,j) -> I
// Nb : It can be (Number of bits of I)/2 atmost
//_______________________________________________________________________________
void Hilbertize ( Int i, 
				  Int j,  
				  Int &I, 
				  Uint Nb )
{
	if (!Nb)
	{
		Nb = sizeof(I)*CHAR_BIT >> 1;
	}
	else if ((sizeof(I)*CHAR_BIT >> 1) < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving bits!" << std::endl;
		
		exit(1);
	}
	
	Int P, Q, t;
	
	// Redo excess work
	for (Uint n = (Nb-1); n > 0; --n)
	{
		Q = 1 << n;
		P = Q - 1;
		
		if (i & Q)
			i ^= P;
		
		if (j & Q)
			i ^= P;
		else
		{
			t = (i ^ j) & P;
			i ^= t;
			j ^= t;
		}
	}
	
	// Gray encode
	j ^= i;
	
	t = 0;
	
	for (Q = 1 << (Nb-1); Q > 1; Q >>= 1)
		if (j & Q)
			t ^= Q-1;
	
	i ^= t;
	j ^= t;
	
	// Mortonize
	I = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (2*n);
		
		if ( j & M0 )
			I |= M1;
		
		M1 <<= 1;
		
		if ( i & M0 )
			I |= M1;
	}
}

//_______________________________________________________________________________
// (i,j,k) -> I
// Nb : It can be floor((Number of bits of I)/3) atmost
//_______________________________________________________________________________
void Hilbertize ( Int i, 
				  Int j,  
				  Int k, 
				  Int &I, 
				  Uint Nb )
{
	if (!Nb)
	{
		Nb = (sizeof(I)*CHAR_BIT)/3;
	}
	else if ((sizeof(I)*CHAR_BIT)/3 < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving bits!" << std::endl;
		
		exit(1);
	}
	
	Int P, Q, t;
	
	for (Uint n = (Nb-1); n > 0; --n)
	{
		Q = 1 << n;
		P = Q - 1;
		
		if (i & Q)
			i ^= P;
		
		if (j & Q)
			i ^= P;
		else
		{
			t = (i ^ j) & P;
			i ^= t;
			j ^= t;
		}
		
		if (k & Q)
			i ^= P;
		else
		{
			t = (i ^ k) & P;
			i ^= t;
			k ^= t;
		}
	}
	
	// Gray encode
	j ^= i;
	k ^= j;
	
	t = 0;
	
	for (Q = 1 << (Nb-1); Q > 1; Q >>= 1)
		if (k & Q)
			t ^= Q-1;
	
	i ^= t;
	j ^= t;
	k ^= t;
	
	// Mortonize
	I = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (3*n);
		
		if ( k & M0 )
			I |= M1;
		
		M1 <<= 1;
		
		if ( j & M0 )
			I |= M1;
		
		M1 <<= 1;
		
		if ( i & M0 )
			I |= M1;
	}
}

//_______________________________________________________________________________
// I -> (i,j)
// Nb : It can be (Number of bits of I)/2 atmost
//_______________________________________________________________________________
void DeHilbertize ( Int &i, 
				    Int &j,  
				    Int I, 
				    Uint Nb )
{
	if (!Nb)
	{
		Nb = sizeof(I)*CHAR_BIT >> 1;
	}
	else if ((sizeof(I)*CHAR_BIT >> 1) < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving Nb bits!" << std::endl;
		
		exit(1);
	}
	
	// De-Mortonize
	i = j = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (2*n);
		
		if ( I & M1 )
			j |= M0;
		
		M1 <<= 1;
		
		if ( I & M1 )
			i |= M0;
	}
	
	Int P, Q, t;
	
	t = j >> 1;
	
	j ^= i;
	i ^= t;
	
	for (Uint n{1}; n < Nb; ++n)
	{
		Q = 1 << n;
		P = Q - 1;
		
		if (j & Q)
			i ^= P;
		else
		{
			t = (i ^ j) & P;
			i ^= t;
			j ^= t;
		}
		
		if (i & Q)
			i ^= P;
	}
}

//_______________________________________________________________________________
// I -> (i,j,k)
// Nb : It can be floor((Number of bits of I)/3) atmost
//_______________________________________________________________________________
void DeHilbertize ( Int &i, 
				    Int &j, 
					Int &k,  
				    Int I, 
				    Uint Nb )
{
	if (!Nb)
	{
		Nb = (sizeof(I)*CHAR_BIT)/3;
	}
	else if ((sizeof(I)*CHAR_BIT)/3 < Nb)
	{
		std::cout << "Int data type cannot accommodate the number generated by interleaving Nb bits!" << std::endl;
		
		exit(1);
	}
	
	// De-Mortonize
	i = j = k = 0;
	
	Int M0, M1;
	
	for (Uint n{}; n < Nb; ++n)
	{
		M0 = 1 << n;
		M1 = 1 << (3*n);
		
		if ( I & M1 )
			k |= M0;
		
		M1 <<= 1;
		
		if ( I & M1 )
			j |= M0;
		
		M1 <<= 1;
		
		if ( I & M1 )
			i |= M0;
	}
	
	Int P, Q, t;
	
	t = k >> 1;
	
	k ^= j;
	j ^= i;
	
	i ^= t;
	
	for (Uint n{1}; n < Nb; ++n)
	{
		Q = 1 << n;
		P = Q - 1;
		
		if (k & Q)
			i ^= P;
		else
		{
			t = (i ^ k) & P;
			i ^= t;
			k ^= t;
		}
		
		if (j & Q)
			i ^= P;
		else
		{
			t = (i ^ j) & P;
			i ^= t;
			j ^= t;
		}
		
		if (i & Q)
			i ^= P;
	}
}

//_______________________________________________________________________________
// Plot 2D Morton curve in html
//_______________________________________________________________________________
void DrawMorton2D ( const Int N )
{
	std::ofstream WriteFile("Output/figure.js", std::ios::out);
	WriteFile.flags( std::ios::dec );
	WriteFile.precision(4);
	
	if ( !WriteFile )
	{
		std::cout << "Output file couldnot be opened!" << std::endl;
		
		exit(1);
	}
	
	Int ie, je, is = 0, js = 0;
	
	Uint Nb = CountBits(N);
	
	double h = 500.0/N;
	
	WriteFile << "function draw()" << std::endl << "{" << std::endl << "\tconst canvas = document.querySelector('#canvas');" << std::endl << std::endl;
	
    WriteFile << "\tif (!canvas.getContext)" << std::endl << "\t\treturn;" << std::endl << std::endl;
	
    WriteFile << "\tconst ctx = canvas.getContext('2d');" << std::endl;
    WriteFile << "\tctx.strokeStyle = 'black';" << std::endl;
    WriteFile << "\tctx.lineWidth = 1;" << std::endl << std::endl;
	
	for (Int ii{}; ii <= N; ++ii)
		WriteFile << "\tctx.beginPath(); ctx.moveTo(" << ii*h << "," << 0 << "); ctx.lineTo(" << ii*h << "," << 500.0 << "); ctx.stroke();" << std::endl;
	
	for (Int ii{}; ii <= N; ++ii)
		WriteFile << "\tctx.beginPath(); ctx.moveTo(" << 0 << "," << ii*h << "); ctx.lineTo(" << 500.0 << "," << ii*h << "); ctx.stroke();" << std::endl;
	
	WriteFile << std::endl << "\tctx.lineWidth = 2;" << std::endl;	
	WriteFile << "\tctx.strokeStyle = 'red';" << std::endl << std::endl;
	
	for (Int I{1}; I < N*N; ++I)
	{
		DeMortonize(ie,je,I,Nb);
		
		WriteFile << "\tctx.beginPath(); ctx.moveTo(" << (is+0.5)*h << "," << (js+0.5)*h << "); ctx.lineTo(" << (ie+0.5)*h << "," << (je+0.5)*h << "); ctx.stroke();" << std::endl;
		
		is = ie;
		js = je;
	}
	
	WriteFile << "}" << std::endl << std::endl << "draw();" << std::endl;
	
	WriteFile.close();
}

//_______________________________________________________________________________
// Plot 3D Morton curve in Tecplot
//_______________________________________________________________________________
void DrawMorton3D ( const Int N )
{
	Int i, j, k;
	
	Uint Nb = CountBits(N);
	
	double h = 1.0/N;
	
	std::ofstream WriteFile("Output/Morton3D.tec", std::ios::out);
	WriteFile.flags( std::ios::fixed );
	WriteFile.precision(4);
	
	if ( !WriteFile )
	{
		std::cout << "Output file couldnot be opened!" << std::endl;
		
		exit(1);
	}
	
    WriteFile << "TITLE = \"Morton curve\"" << std::endl;
    WriteFile << "Variables = \"X\",\"Y\",\"Z\"" << std::endl;
    WriteFile << "Zone N = " << (N*N*N) << ", E = " << (N*N*N-1) << ", DATAPACKING = POINT, ZONETYPE = FELINESEG" << std::endl;
	
	for (Int I{}; I < (N*N*N); ++I)
	{
		DeMortonize(i,j,k,I,Nb);
		
		WriteFile << i*h << "\t" << j*h << "\t" << k*h << std::endl;
	}
	
	for (Int I{1}; I < (N*N*N); ++I)
		WriteFile << I << " " << (I+1) << std::endl;
	
	WriteFile.close();
}

//_______________________________________________________________________________
// Plot 2D Hilbert curve in html
//_______________________________________________________________________________
void DrawHilbert2D ( const Int N )
{
	std::ofstream WriteFile("Output/figure.js", std::ios::out);
	WriteFile.flags( std::ios::dec );
	WriteFile.precision(4);
	
	if ( !WriteFile )
	{
		std::cout << "Output file couldnot be opened!" << std::endl;
		
		exit(1);
	}
	
	Int ie, je, is = 0, js = 0;
	
	Uint Nb = CountBits(N);
	
	double h = 500.0/N;
	
	WriteFile << "function draw()" << std::endl << "{" << std::endl << "\tconst canvas = document.querySelector('#canvas');" << std::endl << std::endl;
	
    WriteFile << "\tif (!canvas.getContext)" << std::endl << "\t\treturn;" << std::endl << std::endl;
	
    WriteFile << "\tconst ctx = canvas.getContext('2d');" << std::endl;
    WriteFile << "\tctx.strokeStyle = 'black';" << std::endl;
    WriteFile << "\tctx.lineWidth = 1;" << std::endl << std::endl;
	
	for (Int ii = 0; ii <= N; ++ii)
		WriteFile << "\tctx.beginPath(); ctx.moveTo(" << ii*h << "," << 0 << "); ctx.lineTo(" << ii*h << "," << 500.0 << "); ctx.stroke();" << std::endl;
	
	for (Int ii = 0; ii <= N; ++ii)
		WriteFile << "\tctx.beginPath(); ctx.moveTo(" << 0 << "," << ii*h << "); ctx.lineTo(" << 500.0 << "," << ii*h << "); ctx.stroke();" << std::endl;
	
	WriteFile << std::endl << "\tctx.lineWidth = 2;" << std::endl;	
	WriteFile << "\tctx.strokeStyle = 'red';" << std::endl << std::endl;
	
	for (Int I = 1; I < N*N; ++I)
	{
		DeHilbertize(ie,je,I,Nb);
		
		WriteFile << "\tctx.beginPath(); ctx.moveTo(" << (is+0.5)*h << "," << (js+0.5)*h << "); ctx.lineTo(" << (ie+0.5)*h << "," << (je+0.5)*h << "); ctx.stroke();" << std::endl;
		
		is = ie;
		js = je;
	}
	
	WriteFile << "}" << std::endl << std::endl << "draw();" << std::endl;
	
	WriteFile.close();
}

//_______________________________________________________________________________
// Plot 3D Hilbert curve in Tecplot
//_______________________________________________________________________________
void DrawHilbert3D ( const Int N )
{
	Int i, j, k;
	
	Uint Nb = CountBits(N);
	
	double h = 1.0/N;
	
	std::ofstream WriteFile("Output/Hilbert3D.tec", std::ios::out);
	WriteFile.flags( std::ios::fixed );
	WriteFile.precision(4);
	
	if ( !WriteFile )
	{
		std::cout << "Output file couldnot be opened!" << std::endl;
		
		exit(1);
	}
	
    WriteFile << "TITLE = \"Hilbert curve\"" << std::endl;
    WriteFile << "Variables = \"X\",\"Y\",\"Z\"" << std::endl;
    WriteFile << "Zone N = " << (N*N*N) << ", E = " << (N*N*N-1) << ", DATAPACKING = POINT, ZONETYPE = FELINESEG" << std::endl;
	
	for (Int I{}; I < (N*N*N); ++I)
	{
		DeHilbertize(i,j,k,I,Nb);
		
		WriteFile << i*h << "\t" << j*h << "\t" << k*h << std::endl;
	}
	
	for (Int I{1}; I < (N*N*N); ++I)
		WriteFile << I << " " << (I+1) << std::endl;
	
	WriteFile.close();
}
