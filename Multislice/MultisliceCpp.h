#pragma once
#include "XArray2D.h"
#include "XArray3D.h"
#include "IXAHWave.h"
#include "XA_head2.h"
#include "XA_data.h"
#include "XA_fft2.h"

namespace xar
{
	void MultisliceSphereNF(int nx, int ny, vector<double> vHead, vector<dcomplex> nc, vector<double> R, vector<double> xr, vector<double> yr,
		vector<double> zr, double zlo, double zhi, int nslices, const char* outfilename);
}
