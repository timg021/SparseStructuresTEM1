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

//!!! NOTE that fftw3 and XArray3D have the same notation for the order of the array dimensions. Both use C-style array structure, i.e. the last index changes the fastest, 
//i.e. in XArray3D(dim1, dim2, dim3) the fastest changing dimension is dim3, and in fftw_r2c_3d(n0, n1, n2) the fastest dimension is n2.
//Note also that these dimensions are associated with the following physical coordinates by default: dim1(n0) <-> nz, dim2(n1) <-> ny, dim3(n2) <-> nx.

#define TEST_RUN 0

using namespace xar;

vector<string> FileNames(index_t nangles, index_t ndefocus, string filenamebase);
double FindMax(XArray3D<float>& aaa, index_t karad, index_t jarad, index_t iarad, index_t& kmax, index_t& jmax, index_t& imax);
int mod(int n, index_t m) { return (n - int(index_t(floor(double(n) / m)) * m)); }

int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting BigBangCT program ...");
		double zmin = -10.0, zmax = 0.0, zstep = 0.0392157;
		index_t ndefocus = 1 + size_t((zmax - zmin) / zstep + 0.5); // number of defocus planes to propagate to at each rotation angle	
		index_t nz = ndefocus, ny = 256, nx = 256;
		index_t nx2 = nx / 2 + 1;
		index_t nangles = 1; // !!! nangles values other than 1 are currently not fully supported in the code below
		index_t natom = 4; // how many atoms to locate
		double atomsize = 2.0; // atom diameter in physical units
		double xmin = 0.0, ymin = 0.0;  // default values - may be overwritten below by data read from input files
		double xstep = 1.0, ystep = 1.0; // default values - may be overwritten below by data read from input files
		double xpos2 = 1.0, ypos2 = 1.0, zpos2 = 1.0; // position of the template atom in the 2nd input 3D array - may be overwritten below by data read from input files
		double wl = 0.025; // wavelength in input file units (usually, Angstroms)
		string filenamebaseIn1("C:\\Users\\TimGu\\Downloads\\TempData\\aaa.grd");
		string filenamebaseIn2("C:\\Users\\TimGu\\Downloads\\TempData\\bbb.grd");
		string filenamebaseOut("C:\\Users\\TimGu\\Downloads\\TempData\\ccc.grd");

		XArray3D<float> aaa(nz, ny, nx);
		XArray3D<xar::fcomplex> ccc(nz, ny, nx2);

		//allocate space and create FFTW plans
		Fftwf3drc fftf((int)nz, (int)ny, (int)nx);

		// first array to transform
		aaa.Fill(0);
#if TEST_RUN
		aaa[3][5][7] = 55.0f; // delta-function
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
				if (kk == 0) ny = inten.GetDim1(); 
				else if (inten.GetDim1() != ny) throw std::runtime_error("different ny dimension in input file");
				if (kk == 0) nx = inten.GetDim2(); 
				else if (inten.GetDim2() != nx) throw std::runtime_error("different nx dimension in input file"); 
				for (index_t jj = 0; jj < ny; jj++)
					for (index_t ii = 0; ii < nx; ii++)
						aaa[kk][jj][ii] = inten[jj][ii] - 1.0f; // can take log() instead;
			}
		}
		xmin = GetXlo(inten);
		ymin = GetYlo(inten);
		xstep = GetXStep(inten);
		ystep = GetYStep(inten);
		printf("\nDimensions of input images (nx,ny,nz) = (%zd, %zd, %zd); minima = (%g, %g, %g); steps = (%g, %g, %g).", nx, ny, nz, xmin, ymin, zmin, xstep, ystep, zstep);
#endif

		// FFT of 1st array
		printf("\nFFT of the 1st 3D set ...");
		fftf.SetRealXArray3D(aaa);
#if TEST_RUN		
		//fftf.PrintRealArray("\nBefore 1st FFT:");
#endif
		fftf.ForwardFFT();
#if TEST_RUN		
		//fftf.PrintComplexArray("\nAfter 1st FFT:");
#endif

		// store away the result of the 1st FFT
		fftf.GetComplexXArray3D(ccc);

		// second array to transform
		aaa.Fill(0);
#if TEST_RUN
		aaa[2][2][2] = 1.0; // delta-function
#else
		printf("\nReading 2nd set of input files %s ...", filenamebaseIn2.c_str());
		infiles = FileNames(nangles, ndefocus, filenamebaseIn2);
		for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
		{
			for (index_t kk = 0; kk < ndefocus; kk++)
			{
				XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
				if (inten.GetDim1() != ny) throw std::runtime_error("different ny dimension in input file");
				if (inten.GetDim2() != nx) throw std::runtime_error("different nx dimension in input file");
				for (index_t jj = 0; jj < ny; jj++)
					for (index_t ii = 0; ii < nx; ii++)
						aaa[kk][jj][ii] = inten[jj][ii] - 1.0f; // can take log() instead;;
			}
		}
		// find the centre of gravity
		double integ = 0.0, xpos = 0.0, ypos = 0.0, zpos = 0.0, dtemp;
		for (index_t kk = 0; kk < ndefocus; kk++)
			for (index_t jj = 0; jj < ny; jj++)
				for (index_t ii = 0; ii < nx; ii++)
				{
					dtemp = abs(aaa[kk][jj][ii]);
					integ += dtemp;
					xpos += dtemp * ii;
					ypos += dtemp * jj;
					zpos += dtemp * kk;
				}
		xpos /= integ; ypos /= integ; zpos /= integ;
		index_t ipos = index_t(xpos + 0.5), jpos = index_t(ypos + 0.5), kpos = index_t(zpos + 0.5);
		xpos2 = xmin + xstep * ipos; ypos2 = ymin + ystep * jpos; zpos2 = zmin + zstep * kpos;
		printf("\nCentre of gravity position of the 2nd array in pixels = (%zd, %zd, %zd), and in physical units = (%g, %g, %g).", ipos, jpos, kpos, xpos2, ypos2, zpos2);
#endif
		
		// FFT of the 2nd array
		printf("\nFFT of the 2nd 3D set ...");
		fftf.SetRealXArray3D(aaa);
#if TEST_RUN		
		//fftf.PrintRealArray("\nBefore 2nd FFT:");
#endif
		fftf.ForwardFFT();
#if TEST_RUN		
		//fftf.PrintComplexArray("\nAfter 2nd FFT:");
#endif

		/// multiply FFTs of 2 arrays, taking the conjugate of the second one
		printf("\nMultiplying FFT of the first by the conjugate of the FFT of the second ...");
		float ftemp;
		fftwf_complex* pout = fftf.GetComplex();
		int m = 0;
		for (index_t k = 0; k < nz; k++)
			for (index_t j = 0; j < ny; j++)
				for (index_t i = 0; i < nx2; i++)
				{
					ftemp = pout[m][0] * ccc[k][j][i].real() + pout[m][1] * ccc[k][j][i].imag();
					pout[m][1] = -pout[m][1] * ccc[k][j][i].real() + pout[m][0] * ccc[k][j][i].imag();
					pout[m][0] = ftemp;
					m++;
				}
#if TEST_RUN		
		//fftf.PrintComplexArray("\nAfter multiplication:");
#endif

		// inverse FFT of the product
		printf("\nInverse FFT of the product ...");
		fftf.InverseFFT();
				
		// get the result
		fftf.GetRealXArray3D(aaa);

#if TEST_RUN		
		//printf("\nAfter inverse FFT:");
		//for (index_t k = 0; k < nz; k++)
		//	for (index_t j = 0; j < ny; j++)
		//		for (index_t i = 0; i < nx; i++)
		//			printf("\naaa[%zd,%zd,%zd] = %g", k, j, i, aaa[k][j][i]);
#else
		printf("\nWriting the output files %s ...", filenamebaseOut.c_str());
		infiles = FileNames(nangles, ndefocus, filenamebaseOut);
		for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
		{
			for (index_t kk = 0; kk < ndefocus; kk++)
			{
				for (index_t jj = 0; jj < ny; jj++)
					for (index_t ii = 0; ii < nx; ii++)
						inten[jj][ii] = aaa[kk][jj][ii];
				XArData::WriteFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), xar::eGRDBIN);
			}
		}

		// find the maximums
		index_t kmax = 0, jmax = 0, imax = 0;
		index_t karad = index_t(atomsize / zstep / 2.0 + 0.5), jarad = index_t(atomsize / ystep / 2.0 + 0.5), iarad = index_t(atomsize / xstep / 2.0 + 0.5);
		for (index_t nn = 0; nn < natom; nn++)
		{
			double amax = FindMax(aaa, karad, jarad, iarad, kmax, jmax, imax);
			printf("\n\nDetected atom number %zd:", nn);
			printf("\nOptimal shift (i,j,k) of the 2nd array to the 1st one in pixels = (%zd, %zd, %zd).", imax, jmax, kmax);
			printf("\nOptimal shift (x,y,z) of the 2nd array to the 1st one in physics units = (%g, %g, %g).", imax * zstep, jmax * ystep, kmax * xstep);
			printf("\nMaximum correlation = %g.", amax);

			double x1 = xstep * nx, xmax = xpos2 + imax * xstep;
			if (xmax > x1) xmax -= x1;
			double y1 = ystep * ny, ymax = ypos2 + jmax * ystep;
			if (ymax > y1) ymax -= y1;
			double z1 = zstep * (nz - 1), zmax = zpos2 + kmax * zstep;
			if (zmax > z1) zmax -= z1;
			printf("\nAbsolute position (x,y,z) of the detected feature in 1st array in physics units = (%g, %g, %g).", xmax, ymax, zmax);
		}

#endif
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
// Creates a sequence of file names properly indexed by rotation angles and defocus distances (using the same algorithm as in MultisliceK.cpp)
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


double FindMax(XArray3D<float>& aaa, index_t karad, index_t jarad, index_t iarad, index_t& kmax, index_t& jmax, index_t& imax)
// Finds the value and position of the maximum in a 3D array, and fills a 3D vicinity of that point with zeros (to enable the search of subsequent maximums)
{
	kmax = jmax = imax = 0;
	double amax = aaa[kmax][jmax][imax];
	index_t nz = aaa.GetDim1(), ny = aaa.GetDim2(), nx = aaa.GetDim3();
	for (index_t kk = 0; kk < nz; kk++)
		for (index_t jj = 0; jj < ny; jj++)
			for (index_t ii = 0; ii < nx; ii++)
				if (aaa[kk][jj][ii] > amax)
				{
					amax = aaa[kk][jj][ii];
					kmax = kk; jmax = jj; imax = ii;
				}
	
	int kk1, jj1, ii1;
	for (int kk = int(kmax)- int(karad); kk < kmax + karad; kk++)
	{
		kk1 = mod(kk, nz); assert(kk1 >= 0 && kk1 < nz);
		for (int jj = int(jmax) - int(jarad); jj < jmax + jarad; jj++)
		{
			jj1 = mod(jj, ny); assert(jj1 >= 0 && jj1 < ny);
			for (int ii = int(imax) - int(iarad); ii < imax + iarad; ii++)
			{
				ii1 = mod(ii, nx); assert(ii1 >= 0 && ii1 < nx);
				aaa[kk1][jj1][ii1] = 0.0f;
			}
		}
	}

	return amax;
}