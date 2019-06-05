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
	string outfilename("C:\\Users\\tgureyev\\Downloads\\Cpp\\aaa.grd");
	const size_t nslices = 2;
	const size_t nangles = 1;
	const size_t nx(512), ny(512);
	vector<double> vHead(5); // analogue of Wavehead2D
	const double energ(1.0); // E in keV
	const double LengthScale(1.0); //scaling factor for length-type parameters
	vHead[0] = 12.398E-4 / energ; //wl in microns
	vHead[1] = -LengthScale; //xlo
	vHead[2] = LengthScale; //xhi
	vHead[3] = -LengthScale; //ylo
	vHead[4] = LengthScale; //yhi
	const double zlo(-LengthScale), zhi(LengthScale);

	const size_t nSpheres(10);
	XArray1D<double> R_in(nSpheres), x_in(nSpheres), y_in(nSpheres), z_in(nSpheres), x_in_i(nSpheres), z_in_i(nSpheres);

	R_in[0] = 0.3 * LengthScale;
	for (size_t i = 1; i < nSpheres; i++) R_in[i] = 0.1 * LengthScale;
	
	x_in[0] = 0;
	for (size_t i = 1; i < nSpheres; i++) x_in[i] = pow(-1, i) * double(i % 3 + 1) * 0.2 * LengthScale;

	y_in[0] = -0.5 * LengthScale;
	for (size_t i = 1; i < nSpheres; i++) y_in[i] = 0.5 * LengthScale;

	z_in[0] = 0;
	for (size_t i = 1; i < nSpheres; i++) z_in[i] = pow(-1, int(i / 5)) * double(i % 3 + 1) * 0.2 * LengthScale;

	vector<dcomplex> nc(nSpheres);
	vector<double> beta_in(nSpheres), delta_in(nSpheres);
	double BetaOrder = 1.E-5, DeltaOrder = 1.e-6; //scaling factors for beta and delta
	for (size_t i = 0; i < nSpheres; i++) beta_in[i] = double(i % 5 + 1) * BetaOrder;
	for (size_t i = 0; i < nSpheres; i++) delta_in[i] = double(i % 5 + 1) * DeltaOrder;
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
