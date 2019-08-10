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

#define TEST_RUN 0

using namespace xar;

vector<string> FileNames(index_t nangles, index_t ndefocus, string filenamebase);

int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting BigBangCT program ...");
		index_t nx = 256, ny = 256, nz = 256;
		index_t nz2 = nz / 2 + 1;
		index_t nangles = 1; // !!! nangles values other than 1 are currently not fully supported in the code below
		index_t ndefocus = 256;
		double zDefocusRange = 10.0;
		double zstep = zDefocusRange / double(ndefocus - 1);
		double xstep = 1.0, ystep = 1.0; // default values - may be overwritten below by data read from input files
		double wl = 0.025; // wavelength in input file units (usually, Angstroms)
		string filenamebaseIn1("C:\\Users\\TimGu\\Downloads\\TempData\\aaa.grd");
		string filenamebaseIn2("C:\\Users\\TimGu\\Downloads\\TempData\\bbb.grd");
		string filenamebaseOut("C:\\Users\\TimGu\\Downloads\\TempData\\ccc.grd");

		XArray3D<float> aaa(nx, ny, nz, 0.0);
		XArray3D<xar::fcomplex> ccc(nx, ny, nz2);

		//allocate space and create FFTW plans
		Fftwf3drc fftf((int)nx, (int)ny, (int)nz);

		// first array to transform
#if TEST_RUN
		aaa[0][0][0] = 1.0f; // delta-function
#else
		printf("\nReading 1st set of input files %s ...", filenamebaseIn1.c_str());
		IXAHWave2D* ph2new = CreateWavehead2D();
		XArray2D<float> inten;
		inten.SetHeadPtr(ph2new);
		vector<string> infiles = FileNames(nangles, ndefocus, filenamebaseIn1);
		for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
		{
			for (index_t kk = 0; kk < ndefocus; kk++)
			{
				XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
				for (index_t ii = 0; ii < inten.GetDim1(); ii++)
					for (index_t jj = 0; jj < inten.GetDim2(); jj++)
						aaa[ii][jj][kk] = inten[ii][jj];
			}
		}
		aaa -= 1.0f; // can take log() instead
		xstep = GetXStep(inten);
		ystep = GetYStep(inten);
#endif

		// FFT of 1st array
		printf("\nFFT of the 1st 3D set ...");
		fftf.SetRealXArray3D(aaa);
#if TEST_RUN		
		fftf.PrintRealArray("\nBefore FFT:");
#endif
		fftf.ForwardFFT();
#if TEST_RUN		
		fftf.PrintComplexArray("\nAfter FFT:");
#endif

		// store away the result of the 1st FFT
		fftf.GetComplexXArray3D(ccc);

		// second array to transform
#if TEST_RUN
		aaa[0][0][0] = 1.0; // delta-function
#else
		printf("\nReading 2nd set of input files %s ...", filenamebaseIn2.c_str());
		infiles = FileNames(nangles, ndefocus, filenamebaseIn2);
		for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
		{
			for (index_t kk = 0; kk < ndefocus; kk++)
			{
				XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
				for (index_t ii = 0; ii < inten.GetDim1(); ii++)
					for (index_t jj = 0; jj < inten.GetDim2(); jj++)
						aaa[ii][jj][kk] = inten[ii][jj];
			}
		}
		aaa -= 1.0f; // can also take log instead
#endif
		
		// FFT of the 2nd array
		printf("\nFFT of the 2nd 3D set ...");
		fftf.SetRealXArray3D(aaa);
#if TEST_RUN		
		fftf.PrintRealArray("\nBefore FFT:");
#endif
		fftf.ForwardFFT();
#if TEST_RUN		
		fftf.PrintComplexArray("\nAfter FFT:");
#endif

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
#if TEST_RUN		
		fftf.PrintComplexArray("\nAfter multiplication:");
#endif

		// inverse FFT of the product
		printf("\nInverse FFT of the product ...");
		fftf.InverseFFT();
				
		// get the result
		fftf.GetRealXArray3D(aaa);

#if TEST_RUN		
		printf("\nAfter inverse FFT:");
		for (index_t i = 0; i < nx; i++)
			for (index_t j = 0; j < ny; j++)
				for (index_t k = 0; k < nz; k++)
					printf("\naaa[%zd,%zd,%zd] = %g", i, j, k, aaa[i][j][k]);
#else
		printf("\nWriting the output files %s ...", filenamebaseOut.c_str());
		infiles = FileNames(nangles, ndefocus, filenamebaseOut);
		for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
		{
			for (index_t kk = 0; kk < ndefocus; kk++)
			{
				for (index_t ii = 0; ii < inten.GetDim1(); ii++)
					for (index_t jj = 0; jj < inten.GetDim2(); jj++)
						inten[ii][jj] = aaa[ii][jj][kk];
				XArData::WriteFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), xar::eGRDBIN);
			}
		}
#endif

		index_t dindex = index_t(aaa.Norm(eNormIndexOfMax));
		index_t imax = dindex / (ny * nz);
		index_t jmax = (dindex - imax * ny * nz) / nz;
		index_t kmax = dindex - imax * ny * nz - jmax * nz;

		printf("\nIndexes of the point of maximum correlation are (%zd, %zd, %zd).", imax, jmax, kmax);
		printf("\nCoordinates of the point of maximum correlation is (%g, %g, %g).", imax * xstep, jmax * ystep, kmax * zstep);


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


vector<string> FileNames(index_t nangles, index_t ndefocus, string filenamebase)
{
	if (ndefocus < 1 || nangles < 1)
		throw std::runtime_error("bad number of angles and/or defocus distances in FileNames()");

	char buffer[128];
	string strAngle, outfilename_i, outfilename_j;

	vector<string> vstrfileout(ndefocus * nangles); // vector of full output filenames

	// create formatting string to add properly formatted indexes at the end of the output file names
	index_t i_dot = filenamebase.rfind('.'), nfieldA_length, nfieldB_length;
	char ndig[8];
	string myformat("");
	if (ndefocus > 1)
	{
		nfieldA_length = 1 + index_t(log10(double(ndefocus - 1))); //maximum number of digits corresponding to defocuses in the output file name
		sprintf(ndig, "%zd", nfieldA_length); //convert the calculated maximum number of digits corresponding to defocuses into a string, e.g. 5 into "5"
		myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded defocus indexes into file names
	}
	if (nangles > 1)
	{
		nfieldB_length = 1 + index_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the output file name
		sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
		myformat += "_%0" + string(ndig) + "d"; //construct format string for inserting two 0-padded angle indexes into file names - see usage below
	}

	for (index_t i = 0; i < nangles; i++)
	{
		for (index_t j = 0; j < ndefocus; j++)
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


