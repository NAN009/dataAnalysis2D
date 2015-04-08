
#include <vector>
#define NORM2(x,y) (x*x+y*y)

typedef struct Vector2D
{
	union 
	{
		double x,r;
	};
	union 
	{
		double y, g;
	};

	Vector2D(){}
	Vector2D(double _x, double _y)
	{
		x = _x;
		y = _y;
	}
	Vector2D operator *(const Vector2D &a) const { return Vector2D(x*a.x, y*a.y); }

}Vector2D;