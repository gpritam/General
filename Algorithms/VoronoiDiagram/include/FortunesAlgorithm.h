#pragma once

#ifndef FortunesAlgorithm_H
#define FortunesAlgorithm_H

#include <bits/stdc++.h>
#include "General.h"

//_______________________________________________________________________________
// Data structures
//_______________________________________________________________________________
struct Edge
{
	double X[2][2];
	
	Edge *PreviousEdge, *NextEdge;
};

struct Site
{
	double x, y;
	
	Site *NextSite, *PreviousSite;
	
	int N;
	
	double **X;
	
	Edge *EdgeHead, *EdgeLeaf;
};

struct Event;

struct Arc
{
	double xl, yl, xr, yr;
	
	Site *S;
	
	Arc *PreviousArc, *NextArc;
	
	bool CircleEvent;
	
	Event *CircleEventPointer;
};

struct Event
{
	double y;
	
	bool IsCircleEvent;
	
	Site *S;
	Arc *A;
	Event *PreviousEvent, *NextEvent;
};

void GeneratePoints(double **&Vertices, const double xmin, const double xmax, const double ymin, const double ymax, const int N);
double TriangleArea(double x0, double y0, double x1, double y1, double x2, double y2);
bool IsOnParabola(const double x, const double y, const double xf, const double yf, const double yd);
void IntersectionParabola(Arc *A0, Arc *A1, const double yd, double &x, double &y, const double xmin, const double xmax, const double ymin, const double ymax);
void CheckCircleEvent(Arc **A, const double yd, const double xmin, const double xmax, const double ymin, const double ymax, const double CircleEventAccuracy = 1.0E-8);
void SiteEvent(const double xmin, const double xmax, const double ymin, const double ymax, const double CircleEventAccuracy = 1.0E-8);
void CircleEvent(const double xmin, const double xmax, const double ymin, const double ymax, const double CircleEventAccuracy = 1.0E-8);
void GenerateVoronoiDiagram(double **&Vertices, const double xmin, const double xmax, const double ymin, const double ymax, const int N,    const double CircleEventAccuracy = 1.0E-8);
#endif
