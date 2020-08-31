// PhaseRetrieval.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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
		printf("\nStarting PhaseRetrieval program ...");
		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024];
		FILE* ff0 = fopen("PhaseRetrieval.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PhaseRetrieval.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // Defocus_distances_in_Angstroms
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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // Input_filename_base_of_defocus_series_of_the_sample_in_GRD_format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		printf("\nDefocus series file name base = %s", filenamebaseIn.c_str());
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength parameter from input parameter file.");
		double wl = atof(cparam); // wavelength in input file units (usually, Angstroms). Unfortunately, it is not saved in the GRD files
		printf("\nWavelength = %g (A)", wl);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // Output_file_name_in_GRD_format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name from input parameter file.");
		string filenameOut = cparam;
		printf("\nOutput file name = %s", filenameOut.c_str());

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("Error: the number of parallel threads in input parameter file should be >= 1.");

		fclose(ff0); // close input parameter file

		// minimal reconstruction error for interrupting the reconstruction cycle
		double epsilon(0.01);

		// microscope aberrations - for future extension
		double Caobj(0.0), k2maxo = pow(Caobj * 0.001f / wl, 2.0); // objective aperture, mrad --> rad
		double Cs3(0.0), C3 = double(Cs3 * 1.e+7); // C3 aberration, mm --> Angstroms
		double Cs5(0.0), C5 = double(Cs5 * 1.e+7); // C5 aberration, mm --> Angstroms

		//************************************ end reading input parameters from file

		// read input GRD files and create the initial complex amplitudes
		vector<string> vinfilenames;
		FileNames(1, ndefocus, filenamebaseIn, vinfilenames); // create vector of input filenames
		vector<XArray2D<double>> vint0(ndefocus), vint(ndefocus);
		vector<XArray2D<dcomplex>> vcamp(ndefocus);
		vector<double> vint0_norms(ndefocus);
		for (index_t n = 0; n < ndefocus; n++)
		{
			XArData::ReadFileGRD(vint0[n], vinfilenames[n].c_str(), wl);
			vint0_norms[n] = vint0[n].Norm(eNormL1);
			if (vint0_norms[n] == 0) throw std::exception("Error: input intensity file is empty");
			vint0[n] ^= 0.5; // intensity into real amplitude
			MakeComplex(vint0[n], 0.0, vcamp[n], true);
		}

		// starting point of the IWFR iterations
		double ssej(1.0), ssejm1; // current and previous phase reconstruction error, error tolerance
		vector<double> verr(ndefocus); // reconstruction errors in individual defocus planes
		do
		{
			// remember the reconstruction error from the previous cycle
			ssejm1 = ssej;

			// backpropagate each defocus amplitude to the plane z = 0
			// OMP parallelization here
			for (index_t n = 0; n < ndefocus; n++)
			{
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(-defocus[n], false, k2maxo, C3, C5); // propagate to z = 0
			}
			// have to make sure that all OMP threads have finished

			// average backpropagated complex amplitudes
			for (index_t n = 1; n < ndefocus; n++) vcamp[0] += vcamp[n];
			vcamp[0] /= double(ndefocus);

			// forward propagate the averaged complex amplitude to the defocus planes
			// OMP parallelization here
			for (index_t n = 0; n < ndefocus; n++)
			{
				if (n) vcamp[n] = vcamp[0];
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(defocus[n], false, k2maxo, C3, C5); // propagate to z = z[n]
				Abs(vcamp[n], vint[n]);
				vint[n] -= vint0[n];
				verr[n] = pow(vint[n].Norm(eNormL2), 2.0) / vint0_norms[n];
			}
			// have to make sure that all OMP threads have finished

			// calculate the current reconstruction error
			ssej = 0.0;
			for (index_t n = 0; n < ndefocus; n++) ssej += verr[n];
			ssej /= double(ndefocus);

		} while ((ssejm1 - ssej) > epsilon);

		// write the locations of the detected atoms into an XYZ file compatible with Vesta XYZ format
		FILE* ff = fopen(filenameOut.c_str(), "wt");
		printf("\n\nWriting output file %s in Vesta XYZ format ...\n", filenameOut.c_str());
		fclose(ff);

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