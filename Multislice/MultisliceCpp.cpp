// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "XA_ini.h"
#include "XA_data.h"
#include "XAHWave.h"
#include <stdio.h>

using namespace xar;

int main(void)
{
	double x(5.0);
	const char* s = "C:\\Users\\tgureyev\\Downloads\\aaa.grd";
	double y = x * x;
	printf(s);

	XArray2D<float> xa(5, 5);
	Wavehead2D xh(1.e-4, 0, 1, 0, 1);
	xa.SetHeadPtr(xh.Clone());
	XArData::WriteFileGRD<float>(xa, s, eGRDBIN);

	return 0;
}
