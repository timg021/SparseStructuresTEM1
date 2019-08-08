// BigBangCT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <complex.h>
#include <chrono>
#include <fftw3.h>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "fftwf3drc.h"

using namespace xar;

vector<string> FileNames(size_t nangles, size_t ndefocus, string filenamebase);

int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting BigBangCT program ...");
		index_t nx(256), ny(256), nz(256);
		index_t nz2 = nz / 2 + 1;
		double zDefocusRange = 10.0;
		XArray3D<float> aaa(nx, ny, nz, 0.0);
		XArray3D<xar::fcomplex> ccc(nx, ny, nz2);

		//allocate space and create FFTW plans
		Fftwf3drc fftf((int)nx, (int)ny, (int)nz);

		// first array to transform
		//aaa[0][0][0] = 1.0f; // delta-function
		printf("\nReading 1st set of input files ...");
		IXAHWave2D* ph2new = CreateWavehead2D();
		XArray2D<float> inten;
		inten.SetHeadPtr(ph2new);
		size_t nangles = 1;
		size_t ndefocus = 256;
		string filenamebase("C:\\Users\\TimGu\\Downloads\\TempData\\aaa.grd");
		vector<string> infiles = FileNames(nangles, ndefocus, filenamebase);
		size_t kk = 0;
		for (size_t i = 0; i < nangles; i++)
		{
			for (size_t j = 0; j < ndefocus; j++)
			{
				XArData::ReadFileGRD(inten, infiles[i * ndefocus + j].c_str(), 0.025);
				for (size_t ii = 0; ii < inten.GetDim1(); ii++)
					for (size_t jj = 0; jj < inten.GetDim2(); jj++)
						aaa[ii][jj][kk] = inten[ii][jj];
				kk++;
			}
		}

		double xstep = GetXStep(inten);
		double ystep = GetYStep(inten);
		double zstep = zDefocusRange / (kk - 1);

		// FFT of 1st array
		printf("\nFFT of the 1st 3D set ...");
		fftf.SetRealXArray3D(aaa);
		//fftf.PrintRealArray("\nBefore FFT:");
		fftf.ForwardFFT();
		//fftf.PrintComplexArray("\nAfter FFT:");

		// store away the result of the 1st FFT
		fftf.GetComplexXArray3D(ccc);

		// second array to transform
		printf("\nReading 2nd set of input files ...");
		//for (index_t i = 0; i < nx; i++)
		//	for (index_t j = 0; j < ny; j++)
		//		for (index_t k = 0; k < nz; k++)
		//			aaa[i][j][k] = float(i + j + k);
		filenamebase = "C:\\Users\\TimGu\\Downloads\\TempData\\bbb.grd";
		infiles = FileNames(nangles, ndefocus, filenamebase);
		kk = 0;
		for (size_t i = 0; i < nangles; i++)
		{
			for (size_t j = 0; j < ndefocus; j++)
			{
				XArData::ReadFileGRD(inten, infiles[i * ndefocus + j].c_str(), 0.025);
				for (size_t ii = 0; ii < inten.GetDim1(); ii++)
					for (size_t jj = 0; jj < inten.GetDim2(); jj++)
						aaa[ii][jj][kk] = inten[ii][jj];
				kk++;
			}
		}
		
		// FFT of the 2nd array
		printf("\nFFT of the 2nd 3D set ...");
		fftf.SetRealXArray3D(aaa);
		//fftf.PrintRealArray("\nBefore FFT:");
		fftf.ForwardFFT();
		//fftf.PrintComplexArray("\nAfter FFT:");

		/// multiply FFTs of 2 arrays
		printf("\nMultiplying the two 3D FFTs ...");
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
		//fftf.PrintComplexArray("\nAfter multiplication:");

		// inverse FFT of the product
		printf("\nInverse FFT of the product ...");
		fftf.InverseFFT();
				
		// get the result
		printf("\nAnalysing the result ...");
		fftf.GetRealXArray3D(aaa);
		//printf("\nAfter inverse FFT:");
		//for (index_t i = 0; i < nx; i++)
			//for (index_t j = 0; j < ny; j++)
				//for (index_t k = 0; k < nz; k++)
					//printf("\naaa[%zd,%zd,%zd] = %g", i, j, k, aaa[i][j][k]);

		size_t dindex = size_t(aaa.Norm(eNormIndexOfMax));
		size_t imax = dindex / (ny * nz);
		size_t jmax = (dindex - imax * ny * nz) / nz;
		size_t kmax = dindex - imax * ny * nz - jmax * nz;

		printf("\nPoint of maximum correlation is (%zd, %zd, %zd).", imax, jmax, kmax);
		printf("\nPoint of maximum correlation is (%g, %g, %g).", imax * xstep, jmax * ystep, kmax * zstep);
	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\nMain program finished. Execution time = %I64d s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	printf("\nPress any key to exit..."); getchar();
	return 0;

}


vector<string> FileNames(size_t nangles, size_t ndefocus, string filenamebase)
{
	if (ndefocus < 1 || nangles < 1)
		throw std::runtime_error("bad number of angles and/or defocus distances in FileNames()");

	char buffer[128];
	string strAngle, outfilename_i, outfilename_j;

	vector<string> vstrfileout(ndefocus * nangles); // vector of full output filenames

	// create formatting string to add properly formatted indexes at the end of the output file names
	size_t i_dot = filenamebase.rfind('.'), nfieldA_length, nfieldB_length;
	char ndig[8];
	string myformat("");
	if (ndefocus > 1)
	{
		nfieldA_length = 1 + size_t(log10(double(ndefocus - 1))); //maximum number of digits corresponding to defocuses in the output file name
		sprintf(ndig, "%zd", nfieldA_length); //convert the calculated maximum number of digits corresponding to defocuses into a string, e.g. 5 into "5"
		myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded defocus indexes into file names
	}
	if (nangles > 1)
	{
		nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the output file name
		sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
		myformat += "_%0" + string(ndig) + "d"; //construct format string for inserting two 0-padded angle indexes into file names - see usage below
	}

	for (size_t i = 0; i < nangles; i++)
	{
		for (size_t j = 0; j < ndefocus; j++)
		{
			outfilename_j = filenamebase;
			if (ndefocus == 1 && nangles > 1) sprintf(buffer, myformat.data(), i);
			else if (ndefocus > 1 && nangles == 1) sprintf(buffer, myformat.data(), j);
			else sprintf(buffer, myformat.data(), j, i);
			outfilename_j.insert(i_dot, buffer);
			vstrfileout[i * ndefocus + j] = outfilename_j;
		}
	}

	return vstrfileout;
}


