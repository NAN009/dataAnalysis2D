#include "vtk.h"
#include "geo.h"

typedef struct Point2D
{
	double x, y;
	double val[125];
	int ix[5], iy[5];
	int tetra;
	double  bx, by, bw;
	Point2D(){}
	~Point2D()
	Point2D(double px, double py)
	{
		Set(px,py);
	}
	Point2D(const Point2D& p)
	{
		x = p.x;
		y = p.y;
		tetra = p.tetra;
		bx = p.bx; by = p.by; bw = p.bw;
		for (int i = 0; i < 125; ++i)
			val[i] = p.val[i];
	}
	void Set(const double px,const double py);
	void Update(const double px, const double py);

}Point2D;