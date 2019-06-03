// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "ProjSpheres.h"
#include "XA_ini.h"
#include "XA_data.h"
#include "XAHWave.h"
#include <stdio.h>

using namespace xar;

int main(void)
{
	string outfilename("C:\\Users\\tgureyev\\Downloads\\aaa.grd");
	const size_t nslices = 2;
	const size_t nangles = 360;
	const size_t nx(512), ny(512);
	vector<double> vHead(5);
	const double energ(32.0);
	const double D(1.e+4);
	vHead[0] = 12.398E-4 / energ; //wl
	vHead[1] = -D; //xlo
	vHead[2] = D; //xhi
	vHead[3] = -D; //ylo
	vHead[4] = D; //yhi
	const double zlo(-D), zhi(D);

	const size_t nSpheres(10);

	XArray1D<double> R_in(nSpheres), x_in(nSpheres), y_in(nSpheres), z_in(nSpheres), x_in_i(nSpheres), z_in_i(nSpheres);
	R_in[0] = 0.3; R_in[1] = 0.1; R_in[2] = 0.1; R_in[3] = 0.1; R_in[4] = 0.1;
	R_in[5] = 0.1; R_in[6] = 0.1; R_in[7] = 0.1; R_in[8] = 0.1; R_in[9] = 0.1;
	x_in[0] = 0.5; x_in[1] = -0.8; x_in[2] = -0.6; x_in[3] = -0.4; x_in[4] = -0.2;
	x_in[5] = 0; x_in[6] = 0.2; x_in[7] = 0.4; x_in[8] = 0.6; x_in[9] = 0.8;
	y_in[0] = -0.5; y_in[1] = 0.5; y_in[2] = 0.5; y_in[3] = 0.5; y_in[4] = 0.5;
	y_in[5] = 0.5; y_in[6] = 0.5; y_in[7] = 0.5; y_in[8] = 0.5; y_in[9] = 0.5;
	z_in[0] = 0; z_in[1] = 0; z_in[2] = 0.2; z_in[3] = 0.4; z_in[4] = 0.6;
	z_in[5] = 0.8; z_in[6] = 0.6; z_in[7] = 0.4; z_in[8] = 0.2; z_in[9] = 0;
	R_in *= D; x_in *= D; y_in *= D;  z_in *= D;

	vector<dcomplex> nc(nSpheres);
	vector<double> beta_in(nSpheres), delta_in(nSpheres);
	beta_in[0] = 3.0E-11; beta_in[1] = 4.0E-11; beta_in[2] = 5.0E-11; beta_in[3] = 4.0E-11; beta_in[4] = 3.0E-11;
	beta_in[5] = 4.0E-11; beta_in[6] = 5.0E-11; beta_in[7] = 4.0E-11; beta_in[8] = 3.0E-11; beta_in[9] = 4.0E-11;
	delta_in[0] = 3.0E-07; delta_in[1] = 4.0E-07; delta_in[2] = 5.0E-07; delta_in[3] = 4.0E-07; delta_in[4] = 3.0E-07;
	delta_in[5] = 4.0E-07; delta_in[6] = 5.0E-07; delta_in[7] = 4.0E-07; delta_in[8] = 3.0E-07; delta_in[9] = 4.0E-07;
	for (size_t i = 0; i < nSpheres; i++)
		nc[i] = dcomplex(1.0 - delta_in[i], beta_in[i]);
	
	size_t i_dot = outfilename.rfind('.');
	size_t nfield_length = (nangles == 1) ? 1 : 1 + size_t(log10(double(nangles - 1))); //maximum number of digits in the output file name
	char ndig[8];
	sprintf(ndig, "%zd", nfield_length); //convert the calculated maximum number of digits into a string, e.g. 3 into "3"
	string myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded numbers into file names - see usage below

	double angle;
	char buffer[128];
	string outfilename_i;

	try
	{
		for (size_t i = 0; i < nangles; i++)
		{
			outfilename_i = outfilename;
			sprintf(buffer, myformat.data(), i);
			outfilename_i.insert(i_dot, buffer);
			angle = PI * i / nangles;
			for (size_t k = 0; k < nSpheres; k++)
			{
				x_in_i[k] = x_in[k] * cos(angle) + z_in[k] * sin(angle); //x - position of the centre of the internal sphere at the i - th rotation step
				z_in_i[k] = -x_in[k] * sin(angle) + z_in[k] * cos(angle); //z - position of the centre of the internal sphere at the i - th rotation step
			}
			MultisliceSphereNF(nx, ny, vHead, nc, R_in, x_in_i, y_in, z_in_i, zlo, zhi, nslices, outfilename_i.data());
		}
	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	return 0;
}
