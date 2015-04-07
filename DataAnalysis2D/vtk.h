#include <string>
#include <list>
#include <vector>

typedef struct DataIndex
{
	int x, y;
	DataIndex(int xx = -1, int yy = -1)
	{
		x = xx;
		y = yy;
	}
}DataIndex;
#define  FLAG_BG 1;
#define  FLAG_BORDER -1;

class VTK
{
public:
	int dimX, dimY;
	std::string file_path;
	double *data;
	char *bg_flag; //whether the data point is the background
	double min_val, max_val;

	VTK()
	{
		dimX = dimY = 0;
		file_path = "";
		data = nullptr;
		bg_flag = nullptr;
	}
	~VTK()
	{
		if (data)
			delete[] data;
		if (bg_flag)
			delete[] bg_flag;
	}
	int Load(std::string filename, double threshold = 0);

};