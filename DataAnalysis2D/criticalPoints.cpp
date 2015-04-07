#include"criticalPoints.h"
#include "evaluator.h"
using namespace  std;
extern VTK vtk;
vector<CriticalPoint> criticalPoints, minPoints, saddles,maxPoints;

vector<DataIndex> vecIndices;

#ifndef DATA
#define DATA(x,y) data[x+dimX*y];
#define INDEX index.x+dimX*index.y
#endif

int ConstCriticalPoints()
{
	CriticalPoint cp;
	Point2D p;

#define  constcp toothcp

}