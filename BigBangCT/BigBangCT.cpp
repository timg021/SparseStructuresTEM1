// BigBangCT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <complex.h>
#include <fftw3.h>
#include <XArray3D.h>

int main()
{
	index_t nx(4), ny(4), nz(4);
	index_t nz1 = nz / 2 + 1, nz2 = nz1 * 2;
	index_t nr = nx * ny * nz, nr2 = nx * ny * nz2, nc = nx * ny * nz1;
    std::cout << "Hello World!\n";
	xar::XArray3D<float> aaa(nx, ny, nz2, 0.0);
	xar::XArray3D<float> bbb(nx, ny, nz2, 0.0);

	aaa[0][0][0] = 1.0f;
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
			{
				bbb[i][j][k] = float(i + j + k);
			}
	
	float* inaaa = (float*)fftwf_malloc(sizeof(float) * nr);
	fftwf_complex* outaaa = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nc);

	// FFT of aaa
	index_t m(0);
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				inaaa[m++] = aaa[i][j][k];

	printf("\nBefore FFT:");
	for (index_t m = 0; m < nr; m++)
		printf("\n inaaa[%zd] = %g", m, inaaa[m]);

	fftwf_plan fpa = fftwf_plan_dft_r2c_3d((int)nx, (int)ny, (int)nz, inaaa, outaaa, FFTW_ESTIMATE);
	fftwf_execute(fpa);

	printf("\nAfter FFT:");
	for (index_t m = 0; m < nc; m++)
		printf("\n outaaa[%zd] = (%g, %g)", m, outaaa[m][0], outaaa[m][1]);

	index_t i(0);
	float* paaa = &aaa[0][0][0];
	for (index_t m = 0; m < nc; m++)
	{
		*paaa++ = outaaa[m][0]; *paaa++ = outaaa[m][1];
	}

	// FFT of bbb
	m = 0;
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				inaaa[m++] = bbb[i][j][k];

	printf("\nBefore FFT:");
	for (index_t m = 0; m < nr; m++)
		printf("\n inaaa[%zd] = %g", m, inaaa[m]);

	fftwf_execute(fpa);

	i = 0;
	paaa = &bbb[0][0][0];
	for (index_t m = 0; m < nc; m++)
	{
		*paaa++ = outaaa[m][0]; *paaa++ = outaaa[m][1];
	}

	printf("\nAfter FFT:");
	for (auto it = bbb.begin(); it != bbb.end(); it++)
		printf("\n bbb = %g", *it);

	/// multiply FFTs of 2 arrays
	paaa = &aaa[0][0][0];
	float* pbbb = &bbb[0][0][0];
	index_t m2, m21;
	for (index_t m = 0; m < nc; m++)
	{
		m2 = m * 2; m21 = m2 + 1;
		outaaa[m][0] = paaa[m2] * pbbb[m2] - paaa[m21] * pbbb[m21];
		outaaa[m][1] = paaa[m2] * pbbb[m21] + paaa[m21] * pbbb[m2];
	}

	printf("\nAfter multiplication:");
	for (index_t m = 0; m < nc; m++)
		printf("\n outaaa[%zd] = (%g, %g)", m, outaaa[m][0], outaaa[m][1]);

	// inverse FFT of the product
	fftwf_plan fpa1 = fftwf_plan_dft_c2r_3d((int)nx, (int)ny, (int)nz, outaaa, inaaa, FFTW_ESTIMATE);
	fftwf_execute(fpa1);

	m = 0;
	float fnorm = 1.0f / float(nr);
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				aaa[i][j][k] = inaaa[m++] * fnorm;

	printf("\nAfter inverse FFT:");
	for (index_t i = 0; i < nx; i++)
		for (index_t j = 0; j < ny; j++)
			for (index_t k = 0; k < nz; k++)
				printf("\ni = %zd, j = %zd, k = %zd, aaa[i,j,k] = %g", i, j, k, aaa[i][j][k]);
	
	fftwf_destroy_plan(fpa1);
	fftwf_destroy_plan(fpa);
	fftwf_free(outaaa);
	fftwf_free(inaaa);
}

