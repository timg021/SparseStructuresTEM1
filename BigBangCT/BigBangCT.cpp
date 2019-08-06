// BigBangCT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <complex.h>
#include <fftw3.h>
#include <XArray3D.h>

int main()
{
	index_t nx(4), ny(4), nz(4);
    std::cout << "Hello World!\n";
	xar::XArray3D<float> aaa(nx, ny, nz, 1.0);

	printf("\nBefore FFT:");
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				printf("\ni = %zd, j = %zd, k = %zd, arr[i,j,k] = %g", i, j, k, aaa[i][j][k]);
	fftwf_plan fp = fftwf_plan_dft_r2c_3d(2, 2, 2, &aaa[0][0][0], (fftwf_complex*)& aaa[0][0][0], FFTW_ESTIMATE);
	fftwf_execute(fp);
	printf("\nAfter FFT:");
		for (index_t i = 0; i < nx; i++)
			for (index_t j = 0; j < ny; j++)
				for (index_t k = 0; k < nz; k++)
					printf("\ni = %zd, j = %zd, k = %zd, arr[i,j,k] = %g", i, j, k, aaa[i][j][k]);
}

