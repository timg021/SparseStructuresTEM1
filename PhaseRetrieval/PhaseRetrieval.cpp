// PhaseRetrieval.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <complex.h>
#include <chrono>
#include <omp.h>

#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"


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
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];
		FILE* ff0 = fopen("PhaseRetrieval.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PhaseRetrieval.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // Defocus_distances_in_Angstroms
		std::vector<size_t> vwhite(0); // vector of white spaces separating different defocus distances (there should be exactly one white space before each defocus distance and no spaces at the end)
		for (size_t i = 1; i < strlen(cline); i++) if (isspace(cline[i])) vwhite.push_back(i); // count the number of different defocus distances in the input parameter file
		int ndefocus = vwhite.size(); // number of detected defocus distances
		if (!ndefocus) throw std::exception("Error reading defocus distances from input parameter file.");
		std::vector<char[1024]> vstr_defocus(ndefocus); // vector of strings into which defocus distances will be read
		std::vector<double> defocus(ndefocus); // vector of defocus distances (double precision numbers)
		vwhite.push_back(strlen(cline)); // add one more entry corresponding to the end of the parameter string
		printf("\nDefocus planes positions (%d in total): ", ndefocus);
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
		string filenamebaseIn1 = cparam;
		printf("\nDefocus series file name base = %s", filenamebaseIn1.c_str());
		
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
			throw std::exception("The number of parallel threads in input parameter file should be >= 1.");

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		// calculate some useful parameters
		if (ndefocus < 1) throw std::exception("Error: the number of defocus planes is less than 1.");
		index_t nz = ndefocus;
		printf("\nNumber of defocus planes = %zd.", nz);
		index_t ny, nx, nx2; // x and y sizes of the arrays
		index_t nangles = 1; // !!! nangles values other than 1 are currently not supported in the code below
		double xmin = 0.0, ymin = 0.0;  // default values - may be overwritten below by data read from input files
		double xstep = 1.0, ystep = 1.0; // default values - may be overwritten below by data read from input files
		

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