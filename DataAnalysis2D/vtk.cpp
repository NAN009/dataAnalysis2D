#include <math.h>
#include <iostream>
#include <fstream>
#include "vtk.h"

using namespace  std;

#define DATA(x,y) data[x+dimX*y];
#define INDEX index.x+dimX*index.y

#ifndef EPS
#define EPS 1e-8
#endif

int VTK::Load(string filename, double threshold)
{
	filebuf fb;
	double val;
	if (fb.open(filename, ios::in))
	{
		iostream is(&fb);
		string info;
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		is >> info >> dimX >> dimY;
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		
		unsigned int elems = dimY*dimX;
		data = new double[elems];
		min_val = 1e10;
		max_val = -1e10;
		for (unsigned int i = 0; i < elems;++i)
		{
			is >> val;
			if (val < threshold)
				val = threshold;
			data[i] = (int)((val - threshold));
			if (data[i] < 0)
				data[i] = 0;
			if (data[i] < min_val)
				min_val=data[i];
			if (data[i]>max_val)
				max_val = data[i];
		}
		cout << "Read vtk file(" << filename << ") compleled" << endl;
	}
	else return -1;
	return 0;
}