#include "evaluator.h"
#include "intetrationLine.h"
#include <iostream>
 using namespace  std;
 extern VTK vtk;
 extern vector<CriticalPoint> criticalPoints, minPoints, saddles, maxPoints;
 vector<MorseSmale2D> SaddleToMax, SaddleToMin;

 MorseSmale2D::MorseSmale2D()
 {
	 line[0].clear();
	 line[1].clear();
 }
 MorseSmale2D::~MorseSmale2D()
 {
	 line[0].clear();
	 line[1].clear();
 }

#define  MAXSTEP 0.1
#define EPS 1e-14

 void MorseSmale2D::ConnectSaddleToMax()
 {
	 Point2D p, testp;
	 Vector2D dir, pos, dstep;
	 double p_val, test_val, step, l;
	 double dist;
	 unsigned int path, i, j, num, dimX = vtk.dimX, dimY = vtk.dimY;
	 for (path = 0; path < 2;++path)
	 {
		 line[path].clear();
		 pos.x = point.x; pos.y = point.y;
		 line[path].push_back(pos);
		 num = 0;
		 p.Set(point.x, point.y);
		 p.Update(point.x, point.y);
		 dir = point.axis;
		 if (path>0)
		 {
			 dir.x = -dir.x;
			 dir.y = -dir.y;
		 }
		 while (1)
		 {
			 pos.x = p.x+dir.x*MAXSTEP;
			 pos.y = p.y=dir.y*MAXSTEP;
			 if (pos.x<0 || pos.x>dimX - 1 ||
				 pos.y<0 || pos.y>dimY - 1)
			 {
				 EndIdx[path] = -1;
				 break;
			 }

			 p_val = EvaluateVal(p);
			 while (step>=EPS)
			 {
				 pos.x = p.x + dir.x*step;
				 pos.y = p.x + dir.y*step;
				 testp.Update(pos.x,pos.y);
				 test_val = EvaluateVal(testp);
				 step *= 0.618;
				 if (test_val > p_val)
				 {
					 pos.x = p.x + dir.x*step;
					 pos.y = p.x + dir.y*step;
					 p.Update(pos.x, pos.y);
					 if (p_val > EvaluateVal(p))
						 p.Update(testp.x,testp.y);
					 dir.x = EvaluateGradX(p);
					 dir.y = EvaluateGradY(p);
					 l = NORM2(dir.x, dir.y);
					 if (l > 0)
					 {
						 l = 1 / sqrt(l);
						 dir.x *= l;
						 dir.y *= l;
					 }
					 break;
				 }
			 }
			 //reach a max point
			 if (step < EPS)
			 {
				 dist = MAXSTEP;
				 for (i = 0; i < maxPoints.size(); ++i)
				 {
					 dstep.x = maxPoints[i].x - p.x;
					 dstep.y = maxPoints[i].x - p.y;
					 if (NORM2(dstep.x, dstep.y) < dist)
					 {
						 j = i;
						 dist = NORM2(dstep.x, dstep.y);
					 }
				 }
				 //connect to the max
				 if (dist < MAXSTEP)
				 {
					 pos.x = maxPoints[j].x;
					 pos.y = maxPoints[j].y;
					 line[path].push_back(pos);
					 EndIdx[path] = j;
				 }
				 else
				 {
					 CriticalPoint new_max;
					 new_max.x = pos.x; new_max.y = pos.y;
					 new_max.type = CriticalPoint::Max;
					 maxPoints.push_back(new_max);
					 criticalPoints.push_back(new_max);
					 line[path].push_back(pos);
					 EndIdx[path] = maxPoints.size() - 1;
					 printf(" %.16f,%.16f\n",pos.x,pos.y);
				 }
				 break;
			 }
			 else
			 {
				 dstep.x = (line[path])[num].x - pos.x;
				 dstep.y = (line[path])[num].y - pos.y;
				 if (NORM2(dstep.x, dstep.y) > MAXSTEP*MAXSTEP)
				 {
					 line[path].push_back(pos);
					 ++num;
				 }
			 }			 
		 }
	 }
 }
 void connectAllSaddleToMax()
 {
	 unsigned int i;
	 cout << "Connect Saddle to Max Point" << endl;
	 SaddleToMax.clear();
	 SaddleToMax.resize(saddles.size());
	 for (i = 0; i < saddles.size();++i)
	 {
		 SaddleToMax[i].point = saddles[i];
		 if (saddles[i].type == CriticalPoint::Saddle)
			 SaddleToMax[i].ConnectSaddleToMax();
		 printf("%.1f\r", i * 100 / saddles.size());
	 }
 }