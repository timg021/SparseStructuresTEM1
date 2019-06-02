// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "MultisliceCpp.h"
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

	int nx(512), ny(512);
	vector<double> vHead(5);
	double energ = 1;
	vHead[0] = 12.398E-4 / energ; //wl
	vHead[1] = -1; //xlo
	vHead[2] = 1; //xhi
	vHead[3] = -1; //ylo
	vHead[4] = 1; //yhi

	vector<double> R_in(10), x_in(10), y_in(10), z_in(10);
	R_in[0] = 0.3; R_in[1] = 0.1; R_in[2] = 0.1; R_in[3] = 0.1; R_in[4] = 0.1;
	R_in[5] = 0.1; R_in[6] = 0.1; R_in[7] = 0.1; R_in[8] = 0.1; R_in[9] = 0.1;
	x_in[0] = 0.5; x_in[1] = -0.8; x_in[2] = -0.6; x_in[3] = -0.4; x_in[4] = -0.2;
	x_in[5] = 0; x_in[6] = 0.2; x_in[7] = 0.4; x_in[8] = 0.6; x_in[9] = 0.8;
	y_in[0] = -0.5; y_in[1] = 0.5; y_in[2] = 0.5; y_in[3] = 0.5; y_in[4] = 0.5;
	y_in[5] = 0.5; y_in[6] = 0.5; y_in[7] = 0.5; y_in[8] = 0.5; y_in[9] = 0.5;
	z_in[0] = 0; z_in[1] = 0; z_in[2] = 0.2; z_in[3] = 0.4; z_in[4] = 0.6;
	z_in[5] = 0.8; z_in[6] = 0.6; z_in[7] = 0.4; z_in[8] = 0.2; z_in[9] = 0;

	vector<dcomplex> nc(10);
	vector<double> beta_in(10), delta_in(10);
	beta_in[0] = 3.0E-05; beta_in[1] = 4.0E-05; beta_in[2] = 5.0E-05; beta_in[3] = 4.0E-05; beta_in[4] = 3.0E-05;
	beta_in[5] = 4.0E-05; beta_in[6] = 5.0E-05; beta_in[7] = 4.0E-05; beta_in[8] = 3.0E-05; beta_in[9] = 4.0E-05;
	delta_in[0] = 3.0E-06; delta_in[1] = 4.0E-06; delta_in[2] = 5.0E-06; delta_in[3] = 4.0E-06; delta_in[4] = 3.0E-06;
	delta_in[5] = 4.0E-06; delta_in[6] = 5.0E-06; delta_in[7] = 4.0E-06; delta_in[8] = 3.0E-06; delta_in[9] = 4.0E-06;
	for (size_t i = 0; i < 10; i++)
		nc[i] = dcomplex(1.0 - delta_in[i], beta_in[i]);

	double zlo(-1), zhi(1);
	int nslices = 10;
	
	try
	{
		MultisliceSphereNF(nx, ny, vHead, nc, R_in, x_in, y_in, z_in, zlo, zhi, nslices, s);
	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	return 0;
}
