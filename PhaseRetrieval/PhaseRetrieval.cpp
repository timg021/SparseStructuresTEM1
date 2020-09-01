// PhaseRetrieval.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// The IWFR algorithm is implemented here according to the description in the paper
// L. Allen, "Electron Microscope Cs Correction Using Iterative Wave - Function Reconstruction", Microscopy and Analysis 20(4):15-17 (UK), 2006

#include <complex.h>
#include <chrono>
#include <omp.h>

#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"

using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting IWFR PhaseRetrieval program ...");
		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024];
		FILE* ff0 = fopen("PhaseRetrieval.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PhaseRetrieval.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 1. Defocus_distances_in_Angstroms
		std::vector<size_t> vwhite(0); // vector of white spaces separating different defocus distances (there should be exactly one white space before each defocus distance and no spaces at the end)
		for (size_t i = 1; i < strlen(cline); i++) if (isspace(cline[i])) vwhite.push_back(i); // count the number of different defocus distances in the input parameter file
		size_t ndefocus = vwhite.size(); // number of detected defocus distances
		if (!ndefocus) throw std::exception("Error reading defocus distances from input parameter file.");
		std::vector<char[1024]> vstr_defocus(ndefocus); // vector of strings into which defocus distances will be read
		std::vector<double> defocus(ndefocus); // vector of defocus distances (double precision numbers)
		vwhite.push_back(strlen(cline)); // add one more entry corresponding to the end of the parameter string
		printf("\nDefocus planes positions (%zd in total): ", ndefocus);
		for (size_t j = 0; j < ndefocus; j++)
		{
			for (size_t i = vwhite[j]; i < vwhite[j + 1]; i++)
				vstr_defocus[j][i - vwhite[j]] = cline[i];
			vstr_defocus[j][vwhite[j + 1] - vwhite[j]] = '\0'; // string terminator
			defocus[j] = atof(vstr_defocus[j]);
			printf("%g ", defocus[j]);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2. Input_filename_base_of_defocus_series_of_the_sample_in_GRD_format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		printf("\nDefocus series file name base = %s", filenamebaseIn.c_str());
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3. Wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in input file units (usually, Angstroms). Unfortunately, it is not saved in the GRD files
		printf("\nWavelength = %g (A)", wl);
		if (wl < 0 || wl > 1)
			throw std::exception("Error: wavelength value appears to be wrong.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam); 
		printf("\nObjective aperture = %g (mrad)", aobj);
		if (aobj < 0 || aobj > 1000)
			throw std::exception("Error: objective aperture value appears to be wrong.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. Spherical aberrations Cs3 and Cs5 in mm
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading spherical aberrations from input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> Angstroms
		Cs5 *= 1.e+7; // mm --> Angstroms

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6. Maximal number of IWFR iterations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading maximal number of iterations from input parameter file.");
		int kmax = atoi(cparam);
		printf("\nMaximal number of iterations = %d", kmax);
		if (kmax < 1)
			throw std::exception("Error: the maximal number of iterations should be >= 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Minimal phase reconstruction error
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading minimal phase reconstruction error from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nMinimal phase reconstruciton error = %g", epsilon);
		if (epsilon < 0)
			throw std::exception("Error: minimal phase reconstruction error must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Output file name in GRC format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name from input parameter file.");
		string filenameOut = cparam;
		printf("\nOutput file name = %s", filenameOut.c_str());
		if (GetFileExtension(filenameOut) != string(".GRC"))
			throw std::exception("Error: output file extension must be grc or GRC.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		// define main work objects
		double dtemp; // auxilliary variable
		double ssej(0.0), ssejm1(0.0); // current and previous average reconstruction errors
		vector<double> verr(ndefocus); // reconstruction errors in individual defocus planes
		vector<XArray2D<double>> vint0(ndefocus), vint(ndefocus); // initial and iterated defocused intensities
		vector<XArray2D<dcomplex>> vcamp(ndefocus); // defocused complex amplitudes
		XArray2D<dcomplex> campOut; // complex amplitude in the plane z=0
		vector<double> vint0_L1(ndefocus); // L1 norms of the initial defocused intensities

		// read input GRD files and create the initial complex amplitudes
		vector<string> vinfilenames;
		FileNames(1, ndefocus, filenamebaseIn, vinfilenames); // create vector of input filenames
		for (index_t n = 0; n < ndefocus; n++)
		{
			XArData::ReadFileGRD(vint0[n], vinfilenames[n].c_str(), wl);
			vint0_L1[n] = vint0[n].Norm(eNormL1);
			if (vint0_L1[n] == 0) throw std::exception("Error: input intensity file is empty");
			vint0[n] ^= 0.5; // intensity --> real amplitude
		}

		// start point of IWFR iterations
		for (int k = 0; k < kmax; k++)
		{
			// apply initial or newly reconstructed phases (the second case is equal to restoring the original moduli)
			for (index_t n = 0; n < ndefocus; n++)
				if (k == 0)	MakeComplex(vint0[n], 0.0, vcamp[n], true);
				else
					for (index_t j = 0; j < vcamp[n].GetDim1(); j++)
						for (index_t i = 0; i < vcamp[n].GetDim2(); i++)
						{
							dtemp = abs(vcamp[n][j][i]);
							if (dtemp) vcamp[n][j][i] *= vint0[n][j][i] / dtemp;
							else MakeComplex(vint0[n], 0.0, vcamp[n], true);
						}

			// backpropagate each defocused amplitude to the plane z = 0
			#pragma omp parallel for
			for (int n = 0; n < ndefocus; n++)
			{
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(-defocus[n], false, k2maxo, Cs3, Cs5); // propagate to z = 0
			}

			// average backpropagated complex amplitudes
			campOut = vcamp[0];
			for (index_t n = 1; n < ndefocus; n++) campOut += vcamp[n];
			campOut /= double(ndefocus);

			// forward propagate the averaged complex amplitude to the defocus planes
			#pragma omp parallel for
			for (int n = 0; n < ndefocus; n++)
			{
				vcamp[n] = campOut;
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(defocus[n], false, k2maxo, Cs3, Cs5); // propagate to z = z[n]
				Abs(vcamp[n], vint[n]);
				vint[n] -= vint0[n];
				verr[n] = pow(vint[n].Norm(eNormL2), 2.0) / vint0_L1[n];
			}

			// calculate the current reconstruction error and
			// if the difference with the previous error is smaller than the defined minimum
			// or, if the error started to increase, interrupt the iterations
			ssej = 0.0;
			for (index_t n = 0; n < ndefocus; n++) ssej += verr[n];
			ssej /= double(ndefocus);
			printf("\nIteration number %d; total error SSE_aver = %g\n", k, ssej);
			for (index_t n = 0; n < ndefocus; n++) printf("SSE(%zd) = %g ", n, verr[n]);
			if (k && (ssejm1 - ssej) < epsilon) break;
			else ssejm1 = ssej;
		} 
		// finish point of IWFR iterations

		// write the reconstructed complex amplitude in the plane z=0 into a GRC file
		XArData::WriteFileGRC(campOut, filenameOut.c_str(), eGRCBIN);
	}
	catch (std::exception& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}

	std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\n\nMain program finished. Execution time = %I64d s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	printf("\nPress any key to exit..."); getchar();
	return 0;

}