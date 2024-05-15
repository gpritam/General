#include "FortunesAlgorithm.h"

Event *EventHead = nullptr, *EventTail, *TemporaryEvent;
Arc *ArcHead, *ArcTail, *TemporaryArc;
Site *SiteHead = nullptr, *SiteTail, *TemporarySite;

int EventNumber = 0;

double ylast;

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void GeneratePoints ( double **&Vertices, 
					  const double xmin, 
					  const double xmax, 
					  const double ymin, 
					  const double ymax, 
					  const int N )
{
	double marginx = 0.005*(xmax-xmin), marginy = 0.005*(ymax-ymin);
	
	for (int i{}; i < N; ++i)
	{
		Vertices[i][0] = xmin + (xmax-xmin-2.0*marginx)*RandomNumber() + marginx;
		Vertices[i][1] = ymin + (ymax-ymin-2.0*marginy)*RandomNumber() + marginy;
	}
}

//_______________________________________________________________________________
// This function returns the area of a triangle having vertices 
// (x0,y0), (x1,y1), (x2,y2)
//_______________________________________________________________________________
double TriangleArea ( double x0,
					  double y0,
					  double x1,
					  double y1,
					  double x2,
					  double y2 )
{
	return 0.5*((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
}

//_______________________________________________________________________________
// This program checks if (x,y) is on a parabola whose focus is (xf,yf) and 
// directrix is y = yd 
//_______________________________________________________________________________
bool IsOnParabola ( const double x, 
					const double y, 
					const double xf, 
					const double yf, 
					const double yd )
{
	return (Absolute(0.5*(xf-x)*(xf-x)/(yf-yd) + 0.5*(yf+yd) - y) < 1.0E-12 ? true : false);
}

//_______________________________________________________________________________
// This function gives the intersection point between two parabolas A0 and A1 with 
// common directrix y=yd 
//_______________________________________________________________________________
void IntersectionParabola ( Arc *A0,
							Arc *A1,
							const double yd,
							double &x,
							double &y, 
							const double xmin, 
							const double xmax, 
							const double ymin, 
							const double ymax )
{
	const double epsilon = 1.0E-16;
	
	if ( (A0 != nullptr) && (A1 != nullptr) )
	{
		double xfl = A0->S->x, yfl = A0->S->y;
		double xfr = A1->S->x, yfr = A1->S->y;
		
		if ( (Absolute(yfl-yd) < epsilon) && (Absolute(yfr-yd) > epsilon) )
		{
			x = xfl;
			y = 0.5*(x-xfr)*(x-xfr)/(yfr-yd) + 0.5*(yfr+yd);
		}
		else if ( (Absolute(yfl-yd) > epsilon) && (Absolute(yfr-yd) < epsilon) )
		{
			x = xfr;
			y = 0.5*(xfl-x)*(xfl-x)/(yfl-yd) + 0.5*(yfl+yd);
		}
		else if ( (Absolute(yfl-yd) > epsilon) && (Absolute(yfr-yd) > epsilon) )
		{
			double a = yfl-yfr;
			double b = 2.0*( xfl*(yfr-yd) - xfr*(yfl-yd) );
			double c = (yfr-yfl)*(yfl-yd)*(yfr-yd) + xfl*xfl*(yd-yfr) - xfr*xfr*(yd-yfl);
			
			if (Absolute(a) < epsilon)
				x = -c/b;
			else
			{
				if ( (b*b-4.0*a*c) < epsilon)
					x = -0.5*b/a;
				else
				{
					double xfirst = 0.5*(-b+sqrt(b*b-4.0*a*c))/a;
					double xsecond = 0.5*(-b-sqrt(b*b-4.0*a*c))/a;
					
					x =  (yfl < yfr ? Maximum(xfirst,xsecond) : Minimum(xfirst,xsecond));
				}
			}
			
			y = 0.5*(xfl-x)*(xfl-x)/(yfl-yd) + 0.5*(yfl+yd);
		}
		else if ( (Absolute(yfl-yd) < epsilon) && (Absolute(yfr-yd) < epsilon) )
			ErrorMessage("Two parabolas intersect at infinity!");
	}
	else if ( (A0 == nullptr) && (A1 != nullptr) )
	{
		x = xmin;
		y = 0.5*(A1->S->x-x)*(A1->S->x-x)/(A1->S->y-yd) + 0.5*(A1->S->y+yd);
	}
	else if ( (A0 != nullptr) && (A1 == nullptr) )
	{
		x = xmax;
		y = 0.5*(A0->S->x-x)*(A0->S->x-x)/(A0->S->y-yd) + 0.5*(A0->S->y+yd);
	}
}

//_______________________________________________________________________________
// If the specified parabola is associated with a circle event, 
// this function places this circle event in the queue. 
//_______________________________________________________________________________
void CheckCircleEvent ( Arc **A, 
						const double yd, 
						const double xmin, 
						const double xmax, 
						const double ymin, 
						const double ymax, 
						const double CircleEventAccuracy )
{
	// If *A already has a circle event, then delete it from the doubly linked list of events first
	if ((*A)->CircleEvent == true)
	{
		if ((*A)->CircleEventPointer->PreviousEvent != nullptr)
			(*A)->CircleEventPointer->PreviousEvent->NextEvent = (*A)->CircleEventPointer->NextEvent;
		
		if ((*A)->CircleEventPointer->NextEvent != nullptr)
			(*A)->CircleEventPointer->NextEvent->PreviousEvent = (*A)->CircleEventPointer->PreviousEvent;
		
		delete (*A)->CircleEventPointer;
		
		(*A)->CircleEventPointer = nullptr;
		(*A)->CircleEvent = false;
	}
	
	// Fresh check for circle event
	double x0, y0, x1, y1, x2, y2, x, y, denominator;
	
	if ((*A)->PreviousArc == nullptr)
	{
		x1 = (*A)->S->x;
		y1 = (*A)->S->y;
		
		x2 = (*A)->NextArc->S->x;
		y2 = (*A)->NextArc->S->y;
		
		if ((y1-y2) > 1.0E-16)
		{
			x = xmin;
			y = -(x2-x1)*(x - 0.5*(x1+x2))/(y2-y1) + 0.5*(y1+y2);
		}	
		else
			return;
	}
	
	if ((*A)->NextArc == nullptr)
	{
		x0 = (*A)->PreviousArc->S->x;
		y0 = (*A)->PreviousArc->S->y;
		
		x1 = (*A)->S->x;
		y1 = (*A)->S->y;
		
		if ((y1-y0) > 1.0E-16)
		{
			x = xmax;
			y = -(x1-x0)*(x - 0.5*(x0+x1))/(y1-y0) + 0.5*(y0+y1);
		}
		else
			return;
	}
	
	if (((*A)->NextArc != nullptr) && ((*A)->PreviousArc != nullptr))
	{
		if ((*A)->NextArc->S == (*A)->PreviousArc->S)
			return;
		
		x0 = (*A)->PreviousArc->S->x;
		y0 = (*A)->PreviousArc->S->y;
		
		x1 = (*A)->S->x;
		y1 = (*A)->S->y;
		
		x2 = (*A)->NextArc->S->x;
		y2 = (*A)->NextArc->S->y;
		
		denominator = 2.0*((y1-y0)*(x2-x1) - (y2-y1)*(x1-x0));
		
		if (Absolute(denominator) > 1.0E-16)
		{
			x = ( (y1-y0)*(x2*x2-x1*x1) - (y2-y1)*(x1*x1-x0*x0) + (y1-y0)*(y2-y1)*(y2-y0) )/denominator;
			
			if (Absolute(y1-y0) > 1.0E-16)
				y = (x1-x0)*(0.5*(x0+x1) - x)/(y1-y0) + 0.5*(y0+y1);
			else
				y = (x2-x1)*(0.5*(x1+x2) - x)/(y2-y1) + 0.5*(y1+y2);
		}
		else
			return;
	}
	
	// Check if y is a NaN or, Inf
	if ((y != y) || std::isinf(y)) 
		return;
	
	y -= SquareRootSquaredSum((x-x1),(y-y1));
	
	if ( (y > yd) || (y < (ymin - 2.0*(ymax-ymin))) )
		return;
	
	// Verify this circle event
	double xl, yl, xr, yr, xm, ym;
	
	IntersectionParabola((*A)->PreviousArc,*A,y,xl,yl,xmin,xmax,ymin,ymax);
	IntersectionParabola(*A,(*A)->NextArc,y,xr,yr,xmin,xmax,ymin,ymax);
	IntersectionParabola((*A)->PreviousArc,(*A)->NextArc,y,xm,ym,xmin,xmax,ymin,ymax);
	
	if ( (Absolute(xl-xm) > CircleEventAccuracy) || 
		 (Absolute(yl-ym) > CircleEventAccuracy) || 
		 (Absolute(xr-xm) > CircleEventAccuracy) || 
		 (Absolute(yr-ym) > CircleEventAccuracy) )
		return;
	
	// It has a circle event at y. Insert it in the queue.
	if ( (TemporaryEvent = new Event) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	TemporaryEvent->y = y;
	TemporaryEvent->IsCircleEvent = true;
	
	TemporaryEvent->S = nullptr;
	TemporaryEvent->A = *A;
	
	bool loop = true;
	
	EventTail = EventHead;
	
	while (loop)
	{
		if (EventTail->y >= y)
		{
			TemporaryEvent->PreviousEvent = EventTail;
			EventTail = EventTail->NextEvent;
			
			if (EventTail == nullptr)
			{
				TemporaryEvent->PreviousEvent->NextEvent = TemporaryEvent;
				TemporaryEvent->NextEvent = nullptr;
				
				loop = false;
			}
		}
		else
		{
			TemporaryEvent->NextEvent = EventTail;
			
			EventTail->PreviousEvent = TemporaryEvent;
			TemporaryEvent->PreviousEvent->NextEvent = TemporaryEvent;
			
			loop = false;
		}
	}
	
	(*A)->CircleEvent = true;
	(*A)->CircleEventPointer = TemporaryEvent;
}

//_______________________________________________________________________________
// This fuction performs all the necessary actions for a site event
//_______________________________________________________________________________
void SiteEvent ( const double xmin, 
				 const double xmax, 
				 const double ymin, 
				 const double ymax, 
				 const double CircleEventAccuracy )
{
	// Step 1 : Find the arc just above the site
	// Step 2 : Insert a new arc splitting that arc in the beachline
	// Step 3 : Test two neighboring arcs of inserted arc for circle events
	// Step 4 : Delete the current event and go for the next event
	
	// Step 1 : Find the arc just above the site
	ArcTail = ArcHead;
	
	bool loop = true;
	
	while (loop)
	{
		double xl, yl, xr, yr;
		
		IntersectionParabola(ArcTail->PreviousArc,ArcTail,EventHead->y,xl,yl,xmin,xmax,ymin,ymax);
		IntersectionParabola(ArcTail,ArcTail->NextArc,EventHead->y,xr,yr,xmin,xmax,ymin,ymax);
		
		if ( (xl < EventHead->S->x) && (EventHead->S->x <= xr) )
			loop = false;
		else
			ArcTail = ArcTail->NextArc;
		
		if (ArcTail == nullptr)
			ErrorMessage("Beachline has a problem!");
	}
	
	// Step 2 : Insert a new arc splitting that arc in the beachline
	Arc *Arc0, *Arc1;
	
	if ( (Arc0 = new Arc) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	if ( (Arc1 = new Arc) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	Arc0->NextArc = Arc1;
	Arc0->PreviousArc = ArcTail;
	
	Arc1->NextArc = ArcTail->NextArc;
	Arc1->PreviousArc = Arc0;
	
	if (ArcTail->NextArc != nullptr)
		ArcTail->NextArc->PreviousArc = Arc1;
	
	ArcTail->NextArc = Arc0;
	
	Arc0->CircleEvent = false;
	Arc0->CircleEventPointer = nullptr;
	
	Arc1->CircleEvent = false;
	Arc1->CircleEventPointer = nullptr;
	
	Arc0->S = EventHead->S;
	Arc1->S = ArcTail->S;
	
	double x, y;
	
	IntersectionParabola(ArcTail,Arc0,EventHead->y,x,y,xmin,xmax,ymin,ymax);
	
	Arc0->xl = x;
	Arc0->yl = y;
	Arc0->xr = x;
	Arc0->yr = y;
	
	Arc1->xl = x;
	Arc1->yl = y;
	Arc1->xr = ArcTail->xr;
	Arc1->yr = ArcTail->yr;
	
	ArcTail->xr = x;
	ArcTail->yr = y;
	
	// Step 3 : Test two neighboring arcs of inserted arc for circle events
	CheckCircleEvent(&ArcTail,EventHead->y,xmin,xmax,ymin,ymax,CircleEventAccuracy);
	CheckCircleEvent(&Arc1,EventHead->y,xmin,xmax,ymin,ymax,CircleEventAccuracy);
	
	// Step 4 : Delete the current event and go for the next event
	ylast = EventHead->y;
	
	if (EventHead->NextEvent == nullptr)
	{
		delete EventHead;
		
		EventHead = nullptr;
	}
	else
	{
		EventHead = EventHead->NextEvent;
		
		delete EventHead->PreviousEvent;
	}
	
	EventNumber++;
}

//_______________________________________________________________________________
// This fuction performs all the necessary actions for a circle event associated 
// with an arc A. Assume that Al = A->PreviousArc and Ar = A->NextArc.
//_______________________________________________________________________________
void CircleEvent ( const double xmin, 
				   const double xmax, 
				   const double ymin, 
				   const double ymax, 
				   const double CircleEventAccuracy )
{
	// Step 1 : Add edges to the sites of Al, A and Ar. Introduce new edge to Al and Ar
	// Step 2 : Make Al and Ar neighbor to each other and refresh circle event, if any, associated to Al and Ar. Delete A.
	// Step 3 : Delete the current event and go for the next event
	
	// Step 1 : Add edges to the sites of Al, A and Ar. Introduce the new edge to Al and Ar
	TemporaryArc = EventHead->A;
	
	double x, y;
	
	IntersectionParabola(TemporaryArc->PreviousArc,TemporaryArc->NextArc,EventHead->y,x,y,xmin,xmax,ymin,ymax);
	
	Edge *TemporaryEdge;
	
	if ( (TemporaryEdge = new Edge) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	TemporaryEdge->X[0][0] = TemporaryArc->xl;
	TemporaryEdge->X[0][1] = TemporaryArc->yl;
	TemporaryEdge->X[1][0] = x;
	TemporaryEdge->X[1][1] = y;
	
	TemporaryEdge->NextEdge = nullptr;
	TemporaryEdge->PreviousEdge = TemporaryArc->S->EdgeLeaf;
	
	if (TemporaryArc->S->EdgeHead == nullptr)
		TemporaryArc->S->EdgeHead = TemporaryEdge;
	else
		TemporaryArc->S->EdgeLeaf->NextEdge = TemporaryEdge;
	
	TemporaryArc->S->EdgeLeaf = TemporaryEdge;
	
	if ( (TemporaryEdge = new Edge) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	TemporaryEdge->X[0][0] = TemporaryArc->xr;
	TemporaryEdge->X[0][1] = TemporaryArc->yr;
	TemporaryEdge->X[1][0] = x;
	TemporaryEdge->X[1][1] = y;
	
	TemporaryEdge->NextEdge = nullptr;
	TemporaryEdge->PreviousEdge = TemporaryArc->S->EdgeLeaf;
	
	TemporaryArc->S->EdgeLeaf->NextEdge = TemporaryEdge;
	TemporaryArc->S->EdgeLeaf = TemporaryEdge;
	
	if (TemporaryArc->PreviousArc != nullptr)
	{
		if ( (TemporaryEdge = new Edge) == nullptr )
			ErrorMessage("Error: Memory cannot be allocated!");
		
		TemporaryEdge->X[0][0] = TemporaryArc->PreviousArc->xr;
		TemporaryEdge->X[0][1] = TemporaryArc->PreviousArc->yr;
		TemporaryEdge->X[1][0] = x;
		TemporaryEdge->X[1][1] = y;
		
		TemporaryEdge->NextEdge = nullptr;
		TemporaryEdge->PreviousEdge = TemporaryArc->PreviousArc->S->EdgeLeaf;
		
		if (TemporaryArc->PreviousArc->S->EdgeHead == nullptr)
			TemporaryArc->PreviousArc->S->EdgeHead = TemporaryEdge;
		else
			TemporaryArc->PreviousArc->S->EdgeLeaf->NextEdge = TemporaryEdge;
		
		TemporaryArc->PreviousArc->S->EdgeLeaf = TemporaryEdge;
	}
	
	if (TemporaryArc->NextArc != nullptr)
	{
		if ( (TemporaryEdge = new Edge) == nullptr )
			ErrorMessage("Error: Memory cannot be allocated!");
		
		TemporaryEdge->X[0][0] = TemporaryArc->NextArc->xl;
		TemporaryEdge->X[0][1] = TemporaryArc->NextArc->yl;
		TemporaryEdge->X[1][0] = x;
		TemporaryEdge->X[1][1] = y;
		
		TemporaryEdge->NextEdge = nullptr;
		TemporaryEdge->PreviousEdge = TemporaryArc->NextArc->S->EdgeLeaf;
		
		if (TemporaryArc->NextArc->S->EdgeHead == nullptr)
			TemporaryArc->NextArc->S->EdgeHead = TemporaryEdge;
		else
			TemporaryArc->NextArc->S->EdgeLeaf->NextEdge = TemporaryEdge;
		
		TemporaryArc->NextArc->S->EdgeLeaf = TemporaryEdge;
	}
	
	if (TemporaryArc->PreviousArc != nullptr)
	{
		TemporaryArc->PreviousArc->xr = x;
		TemporaryArc->PreviousArc->yr = y;
	}
	
	if (TemporaryArc->NextArc != nullptr)
	{
		TemporaryArc->NextArc->xl = x;
		TemporaryArc->NextArc->yl = y;
	}
	
	// Step 2 : Make Al and Ar neighbor to each other and refresh circle event, if any, associated to Al and Ar. Delete A.
	if (TemporaryArc->PreviousArc != nullptr)
		TemporaryArc->PreviousArc->NextArc = TemporaryArc->NextArc;
	
	if (TemporaryArc->NextArc != nullptr)
		TemporaryArc->NextArc->PreviousArc = TemporaryArc->PreviousArc;
	
	if (TemporaryArc->PreviousArc != nullptr)
		CheckCircleEvent(&(TemporaryArc->PreviousArc),EventHead->y,xmin,xmax,ymin,ymax,CircleEventAccuracy);
	
	if (TemporaryArc->NextArc != nullptr)
		CheckCircleEvent(&(TemporaryArc->NextArc),EventHead->y,xmin,xmax,ymin,ymax,CircleEventAccuracy);
	
	if (TemporaryArc == ArcHead)
		ArcHead = TemporaryArc->NextArc;
	
	delete TemporaryArc;
	
	// Step 3 : Delete the current event and go for the next event
	ylast = EventHead->y;
	
	if (EventHead->NextEvent == nullptr)
	{
		delete EventHead;
		EventHead = nullptr;
	}
	else
	{
		EventHead = EventHead->NextEvent;
		delete EventHead->PreviousEvent;
	}
	
	EventNumber++;
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void GenerateVoronoiDiagram ( double **&Vertices, 
					  	      const double xmin, 
					  	      const double xmax, 
					  	      const double ymin, 
					  	      const double ymax, 
					          const int N, 
						      const double CircleEventAccuracy )
{
	char *s = new char[200];
	
	double marginx = 0.005*(xmax-xmin), marginy = 0.005*(ymax-ymin);
	
	// Prepare initial event queue (remove repeated sites)
	Event *TemporaryEvent;
	Site *TemporarySite;
	
	for (int i = 0; i < N; i++)
	{
		bool same = false;
		
		for (int j = (i+1); j < N; j++)
		{
			if ( (Absolute(Vertices[i][0]-Vertices[j][0]) < marginx) && (Absolute(Vertices[i][1]-Vertices[j][1]) < marginy) )
			{
				same = true;
				break;
			}
		}
		
		if (same == false)
		{
			if ( (TemporarySite = new Site) == nullptr )
				ErrorMessage("Error: Memory cannot be allocated!");
			
			TemporarySite->x = Vertices[i][0];
			TemporarySite->y = Vertices[i][1];
			
			TemporarySite->EdgeHead = nullptr;
			TemporarySite->EdgeLeaf = nullptr;
			
			TemporarySite->NextSite = nullptr;
			
			if (SiteHead == nullptr)
			{
				TemporarySite->PreviousSite = nullptr;
				SiteHead = TemporarySite;
			}
			else
			{
				TemporarySite->PreviousSite = SiteTail;
				SiteTail->NextSite = TemporarySite;
			}
			
			SiteTail = TemporarySite;
			
			if ( (TemporaryEvent = new Event) == nullptr )
				ErrorMessage("Error: Memory cannot be allocated!");
			
			TemporaryEvent->y = Vertices[i][1];
			TemporaryEvent->IsCircleEvent = false;
			
			TemporaryEvent->S = SiteTail;
			TemporaryEvent->A = nullptr;
			
			TemporaryEvent->PreviousEvent = nullptr;
			TemporaryEvent->NextEvent = nullptr;
			
			if (EventHead == nullptr)
			{
				EventHead = TemporaryEvent;
			}
			else
			{
				if (Vertices[i][1] >= EventHead->y)
				{
					TemporaryEvent->NextEvent = EventHead;
					EventHead = TemporaryEvent;
				}
				else
				{
					EventTail = EventHead;
					
					bool search = true;
					
					while (search)
					{
						if (EventTail->y > Vertices[i][1])
						{
							if (EventTail->NextEvent == nullptr)
							{
								// Add this event as the last event
								EventTail->NextEvent = TemporaryEvent;
								
								search = false;
							}
							else
							{
								if ((EventTail->NextEvent)->y > Vertices[i][1])
									EventTail = EventTail->NextEvent;
								else
								{
									// Add this as the next event of 'EventTail'
									TemporaryEvent->NextEvent = EventTail->NextEvent;
									EventTail->NextEvent = TemporaryEvent;
									
									search = false;
								}
							}
						}
					}
				}
			}
		}
	}
	
	EventTail = EventHead;
	
	bool eventlist = true;
	
	while (eventlist)
	{
		if (EventTail->NextEvent != nullptr)
		{
			std::cout << EventTail->y << std::endl;
			
			EventTail->NextEvent->PreviousEvent = EventTail;;
			EventTail = EventTail->NextEvent;
		}
		else
			eventlist = false;
	}
	
	std::cout << EventTail->y << std::endl;
	
	// Prepare initial beach line
	Site *FalseSite;
	
	if ( (FalseSite = new Site) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	FalseSite->x = 0.5*(xmin + xmax);
	FalseSite->y = 2.0*ymax - ymin;
	
	FalseSite->EdgeHead = nullptr;
	FalseSite->EdgeLeaf = nullptr;
	
	if ( (ArcHead = new Arc) == nullptr )
		ErrorMessage("Error: Memory cannot be allocated!");
	
	ArcHead->S = FalseSite;
	ArcHead->PreviousArc = nullptr;
	ArcHead->NextArc = nullptr;
	ArcHead->CircleEvent = false;
	ArcHead->CircleEventPointer = nullptr;
	
	ArcHead->xl = xmin;
	ArcHead->yl = 0.5*(ArcHead->S->x-xmin)*(ArcHead->S->x-xmin)/(ArcHead->S->y-ymax) + 0.5*(ArcHead->S->y+ymax);
	ArcHead->xr = xmax;
	ArcHead->yr = 0.5*(ArcHead->S->x-xmax)*(ArcHead->S->x-xmax)/(ArcHead->S->y-ymax) + 0.5*(ArcHead->S->y+ymax);
	
	// Execute events
	bool flag = true;
	
	while (flag)
	{
		std::cout << "y = " << EventHead->y << ", Event number = " << EventNumber << ", Event type : " << (EventHead->IsCircleEvent ? "Circle event" : "Site event") << std::endl;
		
		if (!EventHead->IsCircleEvent)
			SiteEvent(xmin,xmax,ymin,ymax,CircleEventAccuracy);
		else
			CircleEvent(xmin,xmax,ymin,ymax,CircleEventAccuracy);
		
		if (EventHead == nullptr)
			flag = false;
		
		// Collect all the final edges
		if (flag == false)
		{
			ArcTail = ArcHead;
			
			flag = true;
			
			while (flag)
			{
				if (ArcTail != nullptr)
				{
					Edge *TemporaryEdge;
					
					if ( (TemporaryEdge = new Edge) == nullptr )
						ErrorMessage("Error: Memory cannot be allocated!");
					
					double x, y;
					
					IntersectionParabola(ArcTail->PreviousArc,ArcTail,ylast,x,y,xmin,xmax,ymin,ymax);
					
					TemporaryEdge->X[0][0] = ArcTail->xl;
					TemporaryEdge->X[0][1] = ArcTail->yl;
					TemporaryEdge->X[1][0] = x;
					TemporaryEdge->X[1][1] = y;
					
					TemporaryEdge->NextEdge = nullptr;
					TemporaryEdge->PreviousEdge = ArcTail->S->EdgeLeaf;
					
					ArcTail->S->EdgeLeaf->NextEdge = TemporaryEdge;
					ArcTail->S->EdgeLeaf = TemporaryEdge;
					
					if ( (TemporaryEdge = new Edge) == nullptr )
						ErrorMessage("Error: Memory cannot be allocated!");
					
					IntersectionParabola(ArcTail,ArcTail->NextArc,ylast,x,y,xmin,xmax,ymin,ymax);
					
					TemporaryEdge->X[0][0] = ArcTail->xr;
					TemporaryEdge->X[0][1] = ArcTail->yr;
					TemporaryEdge->X[1][0] = x;
					TemporaryEdge->X[1][1] = y;
					
					TemporaryEdge->NextEdge = nullptr;
					TemporaryEdge->PreviousEdge = ArcTail->S->EdgeLeaf;
					
					ArcTail->S->EdgeLeaf->NextEdge = TemporaryEdge;
					ArcTail->S->EdgeLeaf = TemporaryEdge;
					
					ArcTail = ArcTail->NextArc;
				}
				else
					flag = false;
			}
		}
	}
	
	// Deallocate 'FalseSite'
	if (FalseSite->EdgeHead != nullptr)
	{
		flag = true;
		
		while (flag)
		{
			FalseSite->EdgeLeaf = FalseSite->EdgeHead->NextEdge;
			
			delete FalseSite->EdgeHead;
			
			if (FalseSite->EdgeLeaf == nullptr)
				flag = false;
			else
				FalseSite->EdgeHead = FalseSite->EdgeLeaf;
		}
	}
	
	delete FalseSite;
	
	// Detect any edge outside the domain
	Edge *TemporaryEdge;
	
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		bool loop = true;
		
		TemporaryEdge = TemporarySite->EdgeHead;
		
		while (loop)
		{
			if (TriangleArea(TemporarySite->x,TemporarySite->y,TemporaryEdge->X[0][0],TemporaryEdge->X[0][1],TemporaryEdge->X[1][0],TemporaryEdge->X[1][1]) < 0.0)
			{
				Swap(TemporaryEdge->X[0][0],TemporaryEdge->X[1][0]);
				Swap(TemporaryEdge->X[0][1],TemporaryEdge->X[1][1]);
			}
			
			if ( (TemporaryEdge->X[0][1] > ymax) && (TemporaryEdge->X[1][1] > ymax) )
			{
				// Remove that edge element
				if (TemporaryEdge->PreviousEdge == nullptr)
				{
					TemporarySite->EdgeHead = TemporaryEdge->NextEdge;
					TemporarySite->EdgeHead->PreviousEdge = nullptr;
				}
				else
					TemporaryEdge->PreviousEdge->NextEdge = TemporaryEdge->NextEdge;
				
				if (TemporaryEdge->NextEdge == nullptr)
				{
					TemporarySite->EdgeLeaf = TemporaryEdge->PreviousEdge;
					TemporaryEdge->PreviousEdge->NextEdge = nullptr;
				}
				else
					TemporaryEdge->NextEdge->PreviousEdge = TemporaryEdge->PreviousEdge;
				
				delete TemporaryEdge;
			}
			
			if ( ( (TemporaryEdge->X[0][1] > ymax) && (TemporaryEdge->X[1][1] < ymax) ) || ( (TemporaryEdge->X[0][1] < ymax) && (TemporaryEdge->X[1][1] > ymax) ) )
			{
				// Truncate this edge
				double x0 = TemporaryEdge->X[0][0];
				double y0 = TemporaryEdge->X[0][1];
				
				double x1 = TemporaryEdge->X[1][0];
				double y1 = TemporaryEdge->X[1][1];
				
				double x2 = x0 + (x1-x0)*(ymax-y0)/(y1-y0);
				
				if ( (TemporaryEdge->X[0][1] > ymax) && (TemporaryEdge->X[1][1] < ymax) )
				{
					TemporaryEdge->X[0][0] = x2;
					TemporaryEdge->X[0][1] = ymax;
				}
				
				if ( (TemporaryEdge->X[0][1] < ymax) && (TemporaryEdge->X[1][1] > ymax) )
				{
					TemporaryEdge->X[1][0] = x2;
					TemporaryEdge->X[1][1] = ymax;
				}
			}
			
			if ( (TemporaryEdge->X[0][1] < ymin) && (TemporaryEdge->X[1][1] < ymin) )
			{
				// Remove that edge element
				if (TemporaryEdge->PreviousEdge == nullptr)
				{
					TemporarySite->EdgeHead = TemporaryEdge->NextEdge;
					TemporarySite->EdgeHead->PreviousEdge = nullptr;
				}
				else
					TemporaryEdge->PreviousEdge->NextEdge = TemporaryEdge->NextEdge;
				
				if (TemporaryEdge->NextEdge == nullptr)
				{
					TemporarySite->EdgeLeaf = TemporaryEdge->PreviousEdge;
					TemporaryEdge->PreviousEdge->NextEdge = nullptr;
				}
				else
					TemporaryEdge->NextEdge->PreviousEdge = TemporaryEdge->PreviousEdge;
				
				delete TemporaryEdge;
			}
			
			if ( ( (TemporaryEdge->X[0][1] > ymin) && (TemporaryEdge->X[1][1] < ymin) ) || ( (TemporaryEdge->X[0][1] < ymin) && (TemporaryEdge->X[1][1] > ymin) ) )
			{
				// Truncate this edge
				double x0 = TemporaryEdge->X[0][0];
				double y0 = TemporaryEdge->X[0][1];
				
				double x1 = TemporaryEdge->X[1][0];
				double y1 = TemporaryEdge->X[1][1];
				
				double x2 = x0 + (x1-x0)*(ymin-y0)/(y1-y0);
				
				if ( (TemporaryEdge->X[0][1] < ymin) && (TemporaryEdge->X[1][1] > ymin) )
				{
					TemporaryEdge->X[0][0] = x2;
					TemporaryEdge->X[0][1] = ymin;
				}
				
				if ( (TemporaryEdge->X[0][1] > ymin) && (TemporaryEdge->X[1][1] < ymin) )
				{
					TemporaryEdge->X[1][0] = x2;
					TemporaryEdge->X[1][1] = ymin;
				}
			}
			
			if (TemporaryEdge->NextEdge == nullptr)
				loop = false;
			else
				TemporaryEdge = TemporaryEdge->NextEdge;
		}
		
		if (TemporarySite->NextSite == nullptr)
			flag = false;
		else
			TemporarySite = TemporarySite->NextSite;
	}
	
	// Close all the polygons
	double X[2], epsilon = 1.0E-10;
	
	int points;
	
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		bool loop = true, upper;
		
		points = 0;
		
		TemporaryEdge = TemporarySite->EdgeHead;
		
		while (loop)
		{
			if ( ( (Absolute(TemporaryEdge->X[0][1] - ymax) < epsilon) && (Absolute(TemporaryEdge->X[1][1] - ymax) > epsilon) ) || ( (Absolute(TemporaryEdge->X[0][1] - ymax) > epsilon) && (Absolute(TemporaryEdge->X[1][1] - ymax) < epsilon) ) )
			{
				if (Absolute(TemporaryEdge->X[0][1] - ymax) < epsilon)
					X[points] = TemporaryEdge->X[0][0];
				
				if (Absolute(TemporaryEdge->X[1][1] - ymax) < epsilon)
					X[points] = TemporaryEdge->X[1][0];
				
				upper = true;
				
				points++;
				
				if (points == 2)
					break;
			}
			
			if ( ( (Absolute(TemporaryEdge->X[0][1] - ymin) < epsilon) && (Absolute(TemporaryEdge->X[1][1] - ymin) > epsilon) ) || ( (Absolute(TemporaryEdge->X[0][1] - ymin) > epsilon) && (Absolute(TemporaryEdge->X[1][1] - ymin) < epsilon) ) )
			{
				if (Absolute(TemporaryEdge->X[0][1] - ymin) < epsilon)
					X[points] = TemporaryEdge->X[0][0];
				
				if (Absolute(TemporaryEdge->X[1][1] - ymin) < epsilon)
					X[points] = TemporaryEdge->X[1][0];
				
				upper = false;
				
				points++;
				
				if (points == 2)
					break;
			}
			
			TemporaryEdge = TemporaryEdge->NextEdge;
			
			if (TemporaryEdge == nullptr)
				loop = false;
		}
		
		if (loop)
		{
			// Add an edge
			Edge *LastEdge;
			
			if ( (LastEdge = new Edge) == nullptr )
				ErrorMessage("Error: Memory cannot be allocated!");
			
			LastEdge->X[0][0] = X[0];
			LastEdge->X[0][1] = (upper == true ? ymax : ymin);
			LastEdge->X[1][0] = X[1];
			LastEdge->X[1][1] = (upper == true ? ymax : ymin);
			
			if (TriangleArea(TemporarySite->x,TemporarySite->y,LastEdge->X[0][0],LastEdge->X[0][1],LastEdge->X[1][0],LastEdge->X[1][1]) < 0.0)
			{
				Swap(LastEdge->X[0][0],LastEdge->X[1][0]);
				Swap(LastEdge->X[0][1],LastEdge->X[1][1]);
			}
			
			TemporarySite->EdgeLeaf->NextEdge = LastEdge;
			LastEdge->PreviousEdge = TemporarySite->EdgeLeaf;
			LastEdge->NextEdge = nullptr;
			
			TemporarySite->EdgeLeaf = LastEdge;
		}
		
		TemporarySite = TemporarySite->NextSite;
		
		if (TemporarySite == nullptr)
			flag = false;
	}
	
	// Arrange all the edges (remove edges having zero length or duplicate edges, replace all small segments by a single line)
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		bool loop = true;
		
		TemporaryEdge = TemporarySite->EdgeHead;
		
		while (loop)
		{
			Edge *Edge0 = TemporaryEdge->NextEdge;
			
			bool loop0 = true;
			
			while (loop0)
			{
				if ( (Absolute(TemporaryEdge->X[1][0] - Edge0->X[0][0]) < epsilon) && (Absolute(TemporaryEdge->X[1][1] - Edge0->X[0][1]) < epsilon) )
				{
					if (Edge0 != TemporaryEdge->NextEdge)
					{
						Swap(TemporaryEdge->NextEdge->X[0][0],Edge0->X[0][0]);
						Swap(TemporaryEdge->NextEdge->X[0][1],Edge0->X[0][1]);
						Swap(TemporaryEdge->NextEdge->X[1][0],Edge0->X[1][0]);
						Swap(TemporaryEdge->NextEdge->X[1][1],Edge0->X[1][1]);						
					}
					
					break;
				}
				else
				{
					Edge0 = Edge0->NextEdge;
					
					if (Edge0 == nullptr)
						loop0 = false;
				}
			}
			
			if (loop0 == false)
				ErrorMessage("Next edge is not found!");
			
			TemporaryEdge = TemporaryEdge->NextEdge;
			
			if (TemporaryEdge->NextEdge == nullptr)
				loop = false;
		}
		
		TemporarySite = TemporarySite->NextSite;
		
		if (TemporarySite == nullptr)
			flag = false;
	}
	
	// Sort out zero-length edges
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		bool loop = true;
		
		TemporaryEdge = TemporarySite->EdgeHead;
		
		while (loop)
		{
			
			
			TemporaryEdge = TemporaryEdge->NextEdge;
			
			if (TemporaryEdge == nullptr)
				loop = false;
		}
		
		TemporarySite = TemporarySite->NextSite;
		
		if (TemporarySite == nullptr)
			flag = false;
	}
	
	// Sort out co-linear edges
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		bool loop = true;
		
		TemporaryEdge = TemporarySite->EdgeHead;
		
		while (loop)
		{
			
			
			TemporaryEdge = TemporaryEdge->NextEdge;
			
			if (TemporaryEdge == nullptr)
				loop = false;
		}
		
		TemporarySite = TemporarySite->NextSite;
		
		if (TemporarySite == nullptr)
			flag = false;
	}
	
	// Find out number of edges
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		bool loop = true;
		
		TemporaryEdge = TemporarySite->EdgeHead;
		
		points = 0;
		
		while (loop)
		{
			points++;
			
			TemporaryEdge = TemporaryEdge->NextEdge;
			
			if (TemporaryEdge == nullptr)
				loop = false;
		}
		
		TemporarySite->N = points;
		
		TemporarySite = TemporarySite->NextSite;
		
		if (TemporarySite == nullptr)
			flag = false;
	}
	
	// Remove edges and store polygons in a more compact form
	flag = true;
	
	TemporarySite = SiteHead;
	
	while (flag)
	{
		Allocate(TemporarySite->X,TemporarySite->N,2);
		
		int index = 0;
		
		bool loop = true;
		
		while (loop)
		{
			TemporarySite->X[index][0] = TemporarySite->EdgeHead->X[0][0];
			TemporarySite->X[index][1] = TemporarySite->EdgeHead->X[0][1];
			
			index++;
			
			TemporaryEdge = TemporarySite->EdgeHead->NextEdge;
			
			if (TemporaryEdge == nullptr)
				loop = false;
			else
			{
				delete TemporarySite->EdgeHead;
				
				TemporarySite->EdgeHead = TemporaryEdge;
			}
		}
		
		TemporarySite->EdgeLeaf = nullptr;
		
		TemporarySite = TemporarySite->NextSite;
		
		if (TemporarySite == nullptr)
			flag = false;
	}
	
	// Write tecplot file
	int Ne = 0, Nv = 0;
	
	flag = true;
	
	SiteTail = SiteHead;
	
	while (flag)
	{
		TemporarySite = SiteTail->NextSite;
		
		Ne += SiteTail->N;
		Nv += SiteTail->N+1;
		
		if (TemporarySite == nullptr)
			flag = false;
		else
			SiteTail = TemporarySite;
	}
	
	sprintf(s,"Output/Voronoi.tec");
	
	std::ofstream TecplotWrite(s, std::ios::out);
	TecplotWrite.flags( std::ios::dec | std::ios::scientific );
	TecplotWrite.precision(4);
	
	if ( !TecplotWrite )
		ErrorMessage("Output file couldnot be opened!");
	
	TecplotWrite << "TITLE = \"Voronoi diagram\"\nVariables = \"X\",\"Y\",\"u\"" << std::endl;
	TecplotWrite << "Zone N = " << Nv << ", E = " << Ne << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
	
	flag = true;
	
	SiteTail = SiteHead;
	
	while (flag)
	{
		double u = RandomNumber();
		
		TecplotWrite << SiteTail->x << "\t" << SiteTail->y << "\t" << u << std::endl;
		
		for (int i = 0; i < SiteTail->N; i++)
			TecplotWrite << SiteTail->X[i][0] << "\t" << SiteTail->X[i][1] << "\t" << u << std::endl;
		
		SiteTail = SiteTail->NextSite;
		
		if (SiteTail == nullptr)
			flag = false;
	}
	
	flag = true;
	
	int vertexnumber = 0;
	
	SiteTail = SiteHead;
	
	while (flag)
	{
		for (int i = 0; i < SiteTail->N; i++)
			TecplotWrite << (vertexnumber + 1) << "\t" << (vertexnumber + i + 2) << "\t" << (vertexnumber + (i == (SiteTail->N-1) ? 2 : (i + 3))) << std::endl;
		
		vertexnumber += (SiteTail->N+1);
		
		SiteTail = SiteTail->NextSite;
		
		if (SiteTail == nullptr)
			flag = false;
	}
	
	TecplotWrite.close();
	
	sprintf(s,"Output/Sites.tec");
	
	std::ofstream TecplotSitesWrite(s, std::ios::out);
	TecplotSitesWrite.flags( std::ios::dec | std::ios::scientific );
	TecplotSitesWrite.precision(4);
	
	if ( !TecplotSitesWrite )
		ErrorMessage("Output file couldnot be opened!");
	
	flag = true;
	
	SiteTail = SiteHead;
	
	int Nsites = 0;
	
	while (flag)
	{
		TecplotSitesWrite << SiteTail->x << "\t" << SiteTail->y << std::endl;
		
		SiteTail = SiteTail->NextSite;
		
		if (SiteTail == nullptr)
			flag = false;
		
		Nsites++;
	}
	
	TecplotSitesWrite.close();
	
	std::cout << "Number of sites = " << Nsites << std::endl;
	
	// Deallocate 'Site'
	flag = true;
	
	while (flag)
	{
		SiteTail = SiteHead->NextSite;
		
		Deallocate(SiteHead->X,SiteHead->N,2);
		
		delete SiteHead;
		
		if (SiteTail == nullptr)
			flag = false;
		else
			SiteHead = SiteTail;
	}
	
	// Deallocate 'Arc'
	flag = true;
	
	while (flag)
	{
		ArcTail = ArcHead->NextArc;
		
		delete ArcHead;
		
		if (ArcTail == nullptr)
			flag = false;
		else
			ArcHead = ArcTail;
	}
	
	delete [] s;
}
