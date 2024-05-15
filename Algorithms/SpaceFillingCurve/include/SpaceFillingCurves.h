#pragma once

#ifndef SpaceFillingCurves_H
#define SpaceFillingCurves_H

#include <bits/stdc++.h>

using Int  = unsigned int;
using Uint  = short unsigned int;

Uint CountBits(Int x);
void BinaryToGray(Int &x);
void GrayToBinary(Int &x);
void SwapIntegers(Int &x, Int &y);
bool GetBitAtPosition(Int x, Uint n);
void SetBitAtPosition(Int &x, Uint n);
void ClearBitAtPosition(Int &x, Uint n);
void ToggleBitAtPosition(Int &x, Uint n);
void PrintBits(Int x);
void Mortonize(Int i, Int j, Int &I, Uint Nb = 0);
void Mortonize(Int i, Int j, Int k, Int &I, Uint Nb = 0);
void DeMortonize(Int &i, Int &j, Int I, Uint Nb = 0);
void DeMortonize(Int &i, Int &j, Int &k, Int I, Uint Nb = 0);
void Hilbertize(Int i, Int j, Int &I, Uint Nb = 0);
void Hilbertize(Int i, Int j, Int k, Int &I, Uint Nb = 0);
void DeHilbertize(Int &i, Int &j, Int I, Uint Nb = 0);
void DeHilbertize(Int &i, Int &j, Int &k, Int I, Uint Nb = 0);
void DrawMorton2D(const Int N);
void DrawMorton3D(const Int N);
void DrawHilbert2D(const Int N);
void DrawHilbert3D(const Int N);
#endif
