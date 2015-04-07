#include "evaluator.h"
#include "integralline.h"

using namespace std;

extern VTK vtk;
extern vector<CriticalPoint> criticalpoints, minpoints, saddles1, saddles2, maxpoints;
vector<MorseSmale1D> Saddle2MaxMS1, Saddle1MinMS1, Saddle1MaxMS1;

bool isBackground(Vector3D pos) {
	DataIndex idx;
	signed char bgflag;

	idx.x = int(pos.x + 0.5); idx.y = int(pos.y + 0.5); idx.z = int(pos.z + 0.5);
	bgflag = vtk.bgflag[idx.x + vtk.dimX*idx.y + vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
		return true;
	/*

	idx.x=pos.x; idx.y=pos.y; idx.z=pos.z;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x+1; idx.y=pos.y; idx.z=pos.z;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x; idx.y=pos.y+1; idx.z=pos.z;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x+1; idx.y=pos.y+1; idx.z=pos.z;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x; idx.y=pos.y; idx.z=pos.z+1;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x+1; idx.y=pos.y; idx.z=pos.z+1;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x; idx.y=pos.y+1; idx.z=pos.z+1;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;
	idx.x=pos.x+1; idx.y=pos.y+1; idx.z=pos.z+1;
	bgflag = vtk.bgflag[idx.x+vtk.dimX*idx.y+vtk.dimY*idx.z];
	if (bgflag == FLAGBG || bgflag == FLAGBORDER)
	return true;*/
	return false;
}

MorseSmale1D::MorseSmale1D() {
	line[0].clear();
	line[1].clear();
}

MorseSmale1D::~MorseSmale1D() {
	line[0].clear();
	line[1].clear();
}

#define MAXSTEP 0.1
#define EPS 1e-14
void MorseSmale1D::Connect2SaddleToMax() {
	Point3D p, testp;
	Vector3D dir, pos, dstep;
	double pval, testval, step, l;
	double dist;
	unsigned int path, i, j, num, dimX = vtk.dimX, dimY = vtk.dimY, dimZ = vtk.dimZ;

	for (path = 0; path<2; ++path) {
		line[path].clear();
		pos.x = point.x; pos.y = point.y; pos.z = point.z;
		line[path].push_back(pos);
		num = 0;
		p.Set(point.x, point.y, point.z);
		testp.Set(point.x, point.y, point.z);
		dir = point.axis;
		if (path>0) {
			dir.x = -dir.x; dir.y = -dir.y; dir.z = -dir.z;
		}
		while (1) {
			pos.x = p.x + dir.x*MAXSTEP;
			pos.y = p.y + dir.y*MAXSTEP;
			pos.z = p.z + dir.z*MAXSTEP;
			if (pos.x < 0 || pos.x > dimX - 1 ||
				pos.y < 0 || pos.y > dimY - 1 ||
				pos.z < 0 || pos.z > dimZ - 1)
			{
				EndIdx[path] = -1;
				break;
			}
			pval = EvaluateVal(p);
			step = MAXSTEP;
			while (step >= EPS) {
				pos.x = p.x + dir.x*step;
				pos.y = p.y + dir.y*step;
				pos.z = p.z + dir.z*step;
				testp.Update(pos.x, pos.y, pos.z);
				testval = EvaluateVal(testp);
				step *= 0.618;
				if (testval > pval) {
					pos.x = p.x + dir.x*step;
					pos.y = p.y + dir.y*step;
					pos.z = p.z + dir.z*step;
					p.Update(pos.x, pos.y, pos.z);
					if (pval >= EvaluateVal(p)) {
						p.Update(testp.x, testp.y, testp.z);
					}
					dir = EvaluateGradient(p);
					l = NORM2(dir.x, dir.y, dir.z);
					if (l>0) {
						l = 1.0 / sqrt(l);
						dir.x *= l; dir.y *= l; dir.z *= l;
					}
					break;
				}
			}
			// Reach a max point
			if (step < EPS) {
				dist = MAXSTEP;
				for (i = 0; i<maxpoints.size(); ++i) {
					dstep.x = maxpoints[i].x - p.x;
					dstep.y = maxpoints[i].y - p.y;
					dstep.z = maxpoints[i].z - p.z;
					if (NORM2(dstep.x, dstep.y, dstep.z)<dist) {
						j = i;
						dist = NORM2(dstep.x, dstep.y, dstep.z);
					}
				}
				// connect to the max
				if (dist<MAXSTEP) {
					pos.x = maxpoints[j].x;
					pos.y = maxpoints[j].y;
					pos.z = maxpoints[j].z;
					line[path].push_back(pos);
					EndIdx[path] = j;
				}
				else {
					CriticalPoint newmax;
					newmax.x = pos.x; newmax.y = pos.y; newmax.z = pos.z;
					newmax.type = CriticalPoint::Max;
					maxpoints.push_back(newmax);
					criticalpoints.push_back(newmax);
					line[path].push_back(pos);
					EndIdx[path] = maxpoints.size() - 1;
					printf(" %.16f,%.16f,%.16f\n", pos.x, pos.y, pos.z);
				}
				break;
			}
			else {
				dstep.x = (line[path])[num].x - pos.x;
				dstep.y = (line[path])[num].y - pos.y;
				dstep.z = (line[path])[num].z - pos.z;
				if (NORM2(dstep.x, dstep.y, dstep.z)>MAXSTEP*MAXSTEP) {
					line[path].push_back(pos);
					++num;
				}
			}
		}
	}
}

void MorseSmale1D::Connect1SaddleToMax() {
	Point3D p, testp;
	Vector3D dir, pos, dstep;
	Matrix3x3 hesse, V;
	double pval, testval, step, l, diag[3];
	double dist;
	unsigned int path, i, j, num, dimX = vtk.dimX, dimY = vtk.dimY, dimZ = vtk.dimZ;

	for (path = 0; path<2; ++path) {
		line[path].clear();
		pos.x = point.x; pos.y = point.y; pos.z = point.z;
		line[path].push_back(pos);
		num = 0;
		p.Set(point.x, point.y, point.z);
		testp.Set(point.x, point.y, point.z);
		hesse = EvaluateHesse(p);
		hesse.SVD(V, diag);
		i = 0;
		if (diag[1]>diag[0])
			i = 1;
		if (diag[2]>diag[i])
			i = 2;
		dir.x = V.val[0][i];
		dir.y = V.val[1][i];
		dir.z = V.val[2][i];
		if (path>0) {
			dir.x = -dir.x; dir.y = -dir.y; dir.z = -dir.z;
		}
		while (1) {
			pos.x = p.x + dir.x*MAXSTEP;
			pos.y = p.y + dir.y*MAXSTEP;
			pos.z = p.z + dir.z*MAXSTEP;
			if (pos.x < 0 || pos.x > dimX - 1 ||
				pos.y < 0 || pos.y > dimY - 1 ||
				pos.z < 0 || pos.z > dimZ - 1)
			{
				EndIdx[path] = -1;
				break;
			}
			pval = EvaluateVal(p);
			step = MAXSTEP;
			while (step >= EPS) {
				pos.x = p.x + dir.x*step;
				pos.y = p.y + dir.y*step;
				pos.z = p.z + dir.z*step;
				testp.Update(pos.x, pos.y, pos.z);
				testval = EvaluateVal(testp);
				step *= 0.618;
				if (testval > pval) {
					pos.x = p.x + dir.x*step;
					pos.y = p.y + dir.y*step;
					pos.z = p.z + dir.z*step;
					p.Update(pos.x, pos.y, pos.z);
					if (pval >= EvaluateVal(p)) {
						p.Update(testp.x, testp.y, testp.z);
					}
					dir = EvaluateGradient(p);
					l = NORM2(dir.x, dir.y, dir.z);
					if (l>0) {
						l = 1.0 / sqrt(l);
						dir.x *= l; dir.y *= l; dir.z *= l;
					}
					break;
				}
			}
			// Reach a max point
			if (step < EPS) {
				dist = MAXSTEP;
				for (i = 0; i<maxpoints.size(); ++i) {
					dstep.x = maxpoints[i].x - p.x;
					dstep.y = maxpoints[i].y - p.y;
					dstep.z = maxpoints[i].z - p.z;
					if (NORM2(dstep.x, dstep.y, dstep.z)<dist) {
						j = i;
						dist = NORM2(dstep.x, dstep.y, dstep.z);
					}
				}
				// connect to the max
				if (dist<MAXSTEP) {
					pos.x = maxpoints[j].x;
					pos.y = maxpoints[j].y;
					pos.z = maxpoints[j].z;
					line[path].push_back(pos);
					EndIdx[path] = j;
				}
				else {
					CriticalPoint newmax;
					newmax.x = pos.x; newmax.y = pos.y; newmax.z = pos.z;
					newmax.type = CriticalPoint::Max;
					maxpoints.push_back(newmax);
					criticalpoints.push_back(newmax);
					line[path].push_back(pos);
					EndIdx[path] = maxpoints.size() - 1;
					printf(" %.16f,%.16f,%.16f\n", pos.x, pos.y, pos.z);
				}
				break;
			}
			else {
				dstep.x = (line[path])[num].x - pos.x;
				dstep.y = (line[path])[num].y - pos.y;
				dstep.z = (line[path])[num].z - pos.z;
				if (NORM2(dstep.x, dstep.y, dstep.z)>MAXSTEP*MAXSTEP) {
					line[path].push_back(pos);
					++num;
				}
			}
		}
	}
}

void MorseSmale1D::Connect1SaddleToMin() {
	Point3D p, testp;
	Vector3D dir, pos, dstep;
	double pval, testval, step, l;
	double dist;
	unsigned int path, i, j, num, dimX = vtk.dimX, dimY = vtk.dimY, dimZ = vtk.dimZ;

	for (path = 0; path<2; ++path) {
		line[path].clear();
		pos.x = point.x; pos.y = point.y; pos.z = point.z;
		line[path].push_back(pos);
		num = 0;
		p.Set(point.x, point.y, point.z);
		testp.Set(point.x, point.y, point.z);
		dir = point.axis;
		if (path>0) {
			dir.x = -dir.x; dir.y = -dir.y; dir.z = -dir.z;
		}
		while (1) {
			pos.x = p.x + dir.x*MAXSTEP;
			pos.y = p.y + dir.y*MAXSTEP;
			pos.z = p.z + dir.z*MAXSTEP;
			if (pos.x < 0 || pos.x > dimX - 1 ||
				pos.y < 0 || pos.y > dimY - 1 ||
				pos.z < 0 || pos.z > dimZ - 1)
			{
				line[path].clear();
				EndIdx[path] = -1;
				break;
			}
			pval = EvaluateVal(p);
			step = MAXSTEP;
			while (step >= EPS) {
				pos.x = p.x + dir.x*step;
				pos.y = p.y + dir.y*step;
				pos.z = p.z + dir.z*step;
				testp.Update(pos.x, pos.y, pos.z);
				testval = EvaluateVal(testp);
				step *= 0.618;
				if (testval < pval) {
					pos.x = p.x + dir.x*step;
					pos.y = p.y + dir.y*step;
					pos.z = p.z + dir.z*step;
					p.Update(pos.x, pos.y, pos.z);
					if (pval <= EvaluateVal(p)) {
						p.Update(testp.x, testp.y, testp.z);
					}
					dir = EvaluateGradient(p);
					l = NORM2(dir.x, dir.y, dir.z);
					if (l>0) {
						l = -1.0 / sqrt(l);
						dir.x *= l; dir.y *= l; dir.z *= l;
					}
					break;
				}
			}
			// Reach a max point
			if (step < EPS) {
				dist = MAXSTEP;
				for (i = 0; i<minpoints.size(); ++i) {
					dstep.x = minpoints[i].x - p.x;
					dstep.y = minpoints[i].y - p.y;
					dstep.z = minpoints[i].z - p.z;
					if (NORM2(dstep.x, dstep.y, dstep.z)<dist) {
						j = i;
						dist = NORM2(dstep.x, dstep.y, dstep.z);
					}
				}
				// connect to the min
				if (dist<MAXSTEP) {
					pos.x = minpoints[j].x;
					pos.y = minpoints[j].y;
					pos.z = minpoints[j].z;
					line[path].push_back(pos);
					EndIdx[path] = j;
				}
				else if (isBackground(pos)) {
					line[path].clear();
					EndIdx[path] = -1;
					//line[path].push_back(Vector3D(0,0,0));
				}
				else {
					//line[path].clear();
					//EndIdx[path] = -1;

					CriticalPoint newmin;
					newmin.x = pos.x; newmin.y = pos.y; newmin.z = pos.z;
					newmin.type = CriticalPoint::Max;
					minpoints.push_back(newmin);
					criticalpoints.push_back(newmin);
					line[path].push_back(pos);
					line[path].push_back(Vector3D(0, 0, 0));
					//printf(" %.16f,%.16f,%.16f: %g\n",pos.x,pos.y,pos.z,pval);

				}
				break;
			}
			else {
				dstep.x = (line[path])[num].x - pos.x;
				dstep.y = (line[path])[num].y - pos.y;
				dstep.z = (line[path])[num].z - pos.z;
				if (NORM2(dstep.x, dstep.y, dstep.z)>MAXSTEP*MAXSTEP) {
					line[path].push_back(pos);
					++num;
				}
			}
		}
	}
}
#undef EPS
#undef MAXSTEP

void ConnectAll2SaddleToMax() {
	unsigned int i;

	printf("Connect 2-Saddles to Max Points\n");

	Saddle2MaxMS1.clear();
	Saddle2MaxMS1.resize(saddles2.size());
	for (i = 0; i<saddles2.size(); ++i) {
		Saddle2MaxMS1[i].point = saddles2[i];
		//		printf("%d: %f %f %f\n",i,saddles2[i].x,saddles2[i].y,saddles2[i].z);
		if (saddles2[i].type == CriticalPoint::Saddle2) {
			Saddle2MaxMS1[i].Connect2SaddleToMax();
		}
		printf("%.1f\r", i*100.0 / saddles2.size());
	}
}

void ConnectAll1SaddleToMin() {
	unsigned int i;

	printf("Connect 1-Saddles to Min Points\n");

	Saddle1MinMS1.clear();
	Saddle1MinMS1.resize(saddles1.size());
	for (i = 0; i<saddles1.size(); ++i) {
		Saddle1MinMS1[i].point = saddles1[i];
		//printf("%f %f %f\n",saddles1[i].x,saddles1[i].y,saddles1[i].z);
		// 		if (saddles1[i].x<12.55 || saddles1[i].x>12.56 || saddles1[i].y<26.11 || saddles1[i].y>26.12 || saddles1[i].z<17.98 || saddles1[i].z>17.99)
		// 			continue;
		if (saddles1[i].type == CriticalPoint::Saddle1) {
			Saddle1MinMS1[i].Connect1SaddleToMin();
		}
		if (Saddle1MinMS1[i].line[0].size() + Saddle1MinMS1[i].line[1].size() == 0) {
			Saddle1MinMS1.erase(Saddle1MinMS1.begin() + i);
			saddles1.erase(saddles1.begin() + i);
			--i;
		}
		printf("%.1f\r", i*100.0 / saddles1.size());
	}
}

void ConnectAll1SaddleToMax() {
	unsigned int i;

	printf("Connect 1-Saddles to Max Points\n");

	Saddle1MaxMS1.clear();
	Saddle1MaxMS1.resize(saddles1.size());
	for (i = 0; i<saddles1.size(); ++i) {
		Saddle1MaxMS1[i].point = saddles1[i];
		//		printf("%d: %f %f %f\n",i,saddles2[i].x,saddles2[i].y,saddles2[i].z);
		if (saddles1[i].type == CriticalPoint::Saddle1) {
			Saddle1MaxMS1[i].Connect1SaddleToMax();
		}
		printf("%.1f\r", i*100.0 / saddles1.size());
	}
}