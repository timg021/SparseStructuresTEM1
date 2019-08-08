// Implementation of the fftwf3drc class

#include "fftwf3drc.h"

void Fftwf3drc::Cleanup()
{
	if (aplan != 0) { fftwf_destroy_plan(aplan); aplan = 0; }
	if (bplan != 0) { fftwf_destroy_plan(bplan); bplan = 0; }
	if (pin != 0) { fftwf_free(pin); pin = 0; }
	if (pout != 0) { fftwf_free(pout); pout = 0; }
	nx = ny = nz = 0;
}


void Fftwf3drc::GetRealXArray3D(xar::XArray3D<float>& aaa)
{
	if ((int)aaa.GetDim1() != nx || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nz)
		throw std::runtime_error("dimensions of the output real array are different from the internal one in Fftwf3drc");
	int m = 0;
	float fnorm = 1.0f / float(GetNr());
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				aaa[i][j][k] = pin[m++] * fnorm;
}


void Fftwf3drc::SetRealXArray3D(xar::XArray3D<float> aaa)
{
	if ((int)aaa.GetDim1() != nx || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nz)
		throw std::runtime_error("dimensions of the input real array are different from the internal one in Fftwf3drc");
	int m = 0;
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				pin[m++] = aaa[i][j][k];
}


void Fftwf3drc::GetComplexXArray3D(xar::XArray3D<xar::fcomplex>& aaa)
{
	int nz2 = GetNz2();
	if ((int)aaa.GetDim1() != nx || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nz2)
		throw std::runtime_error("dimensions of the output complex array are different from the internal one in Fftwf3drc");
	int m = 0;
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz2; k++)
			{ 
				aaa[i][j][k] = xar::fcomplex(pout[m][0], pout[m][1]); 
				m++; 
			}
}


void Fftwf3drc::SetComplexXArray3D(xar::XArray3D<xar::fcomplex> aaa)
{
	int nz2 = GetNz2();
	if ((int)aaa.GetDim1() != nx || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nz2)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwf3drc");
	int m = 0;
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz2; k++)
			{
				pout[m][0] = aaa[i][j][k].real();
				pout[m][1] = aaa[i][j][k].imag();
				m++;
			}
}


void Fftwf3drc::PrintRealArray(const char* message)
{
	float* pin = GetReal();
	printf(message);
	int m = 0;
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				printf("\npin[%zd,%zd,%zd] = %g", i, j, k, pin[m++]);
}


void Fftwf3drc::PrintComplexArray(const char* message)
{
	printf(message);
	int m = 0;
	int nz2 = GetNz2();
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz2; k++)
			{ 
				printf("\npout[%zd,%zd,%zd] = (%g, %g)", i, j, k, pout[m][0], pout[m][1]); 
				m++; 
			}
}

