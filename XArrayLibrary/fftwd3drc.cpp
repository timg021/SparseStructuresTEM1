// Implementation of the fftwd3drc class

#include "fftwd3drc.h"

void Fftwd3drc::Cleanup()
{
	if (aplan != 0) { fftw_destroy_plan(aplan); aplan = 0; }
	if (bplan != 0) { fftw_destroy_plan(bplan); bplan = 0; }
	if (pin != 0) { fftw_free(pin); pin = 0; }
	if (pout != 0) { fftw_free(pout); pout = 0; }
	nz = ny = nx = 0;
}


void Fftwd3drc::GetRealXArray3D(xar::XArray3D<double>& aaa)
{
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the output real array are different from the internal one in Fftwd3drc");
	int m = 0;
	double fnorm = 1.0 / double(GetNr());
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				aaa[k][j][i] = pin[m++] * fnorm;
}


void Fftwd3drc::SetRealXArray3D(xar::XArray3D<double> aaa)
{
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the input real array are different from the internal one in Fftwd3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				pin[m++] = aaa[k][j][i];
}


void Fftwd3drc::GetComplexXArray3D(xar::XArray3D<xar::dcomplex>& aaa)
{
	int nx2 = GetNx2();
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx2)
		throw std::runtime_error("dimensions of the output complex array are different from the internal one in Fftwd3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{ 
				aaa[k][j][i] = xar::dcomplex(pout[m][0], pout[m][1]); 
				m++; 
			}
}


void Fftwd3drc::SetComplexXArray3D(xar::XArray3D<xar::dcomplex> aaa)
{
	int nx2 = GetNx2();
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx2)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwd3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				pout[m][0] = aaa[k][j][i].real();
				pout[m][1] = aaa[k][j][i].imag();
				m++;
			}
}


void Fftwd3drc::PrintRealArray(const char* message)
{
	double* pin = GetReal();
	printf(message);
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				printf("\npin[%zd,%zd,%zd] = %g", k, j, i, pin[m++]);
}


void Fftwd3drc::PrintComplexArray(const char* message)
{
	printf(message);
	int m = 0;
	int nx2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				printf("\npout[%zd,%zd,%zd] = (%g, %g)", k, j, i, pout[m][0], pout[m][1]); 
				m++; 
			}
}

