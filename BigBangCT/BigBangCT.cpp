// BigBangCT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <complex.h>
#include <fftw3.h>
#include <XArray3D.h>
#include <fftwf3drc.h>

int main()
{
	try
	{
		printf("\nStarting BigBangCT program ...");
		index_t nx(4), ny(4), nz(4);
		index_t nz2 = nz / 2 + 1;
		xar::XArray3D<float> aaa(nx, ny, nz, 0.0);
		xar::XArray3D<xar::fcomplex> ccc(nx, ny, nz2);

		//allocate space and create FFTW plans
		Fftwf3drc fftf((int)nx, (int)ny, (int)nz);

		// first array to transform
		aaa[0][0][0] = 1.0f; // delta-function

		// FFT of 1st array
		fftf.SetRealXArray3D(aaa);
		fftf.PrintRealArray("\nBefore FFT:");
		fftf.ForwardFFT();
		fftf.PrintComplexArray("\nAfter FFT:");

		// store away the result of the 1st FFT
		fftf.GetComplexXArray3D(ccc);

		// second array to transform
		for (index_t i = 0; i < nx; i++)
			for (index_t j = 0; j < ny; j++)
				for (index_t k = 0; k < nz; k++)
					aaa[i][j][k] = float(i + j + k);

		// FFT of the 2nd array
		fftf.SetRealXArray3D(aaa);
		fftf.PrintRealArray("\nBefore FFT:");
		fftf.ForwardFFT();
		fftf.PrintComplexArray("\nAfter FFT:");

		/// multiply FFTs of 2 arrays
		float ftemp;
		fftwf_complex* pout = fftf.GetComplex();
		int m = 0;
		for (index_t i = 0; i < nx; i++)
			for (index_t j = 0; j < ny; j++)
				for (index_t k = 0; k < nz2; k++)
				{
					ftemp = pout[m][0] * ccc[i][j][k].real() - pout[m][1] * ccc[i][j][k].imag();
					pout[m][1] = pout[m][1] * ccc[i][j][k].real() + pout[m][0] * ccc[i][j][k].imag();
					pout[m][0] = ftemp;
					m++;
				}
		fftf.PrintComplexArray("\nAfter multiplication:");

		// inverse FFT of the product
		fftf.InverseFFT();
				
		// get the result
		fftf.GetRealXArray3D(aaa);
		printf("\nAfter inverse FFT:");
		for (index_t i = 0; i < nx; i++)
			for (index_t j = 0; j < ny; j++)
				for (index_t k = 0; k < nz; k++)
					printf("\naaa[%zd,%zd,%zd] = %g", i, j, k, aaa[i][j][k]);
	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	printf("\nPress any key to exit..."); getchar();
	return 0;

}

