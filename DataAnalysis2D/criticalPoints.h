#include "vtk.h"
#include "geo.h"
class CriticalPoint
{
public:
	enum Type
	{
		Min,
		Saddle,
		Max
	};
	double x, y;
	Type type;
	int property;
	double val;
	Vector2D axis;
	CriticalPoint()
	{
		property = 0;
	}
	

};

int ConstCriticalPoints();
