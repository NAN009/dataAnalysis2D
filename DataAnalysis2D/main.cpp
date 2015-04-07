#include "intetrationLine.h"
#include "vtk.h"
#include <string>
using namespace  std;

extern VTK vtk;
int main()
{
	string path;
	path = "G:/visbyGuWeiSi/visbyZengZuo/Hydrogen.vtk.vtk";

	vtk.Load(path,0);
	printf("%d %d\n", (int)vtk.min_val, (int)vtk.max_val);

	int k;
	k=
}