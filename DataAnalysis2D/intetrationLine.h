#include "criticalPoints.h"


class MorseSmale2D
{
public:
	CriticalPoint point;
	std::vector<Vector2D> line[2];
	int EndIdx[2];
	MorseSmale2D();
	~MorseSmale2D();

	void ConnectSaddleToMax();

};
void connectAllSaddleToMax();