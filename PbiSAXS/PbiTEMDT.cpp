//PbiTEMDT.cpp this file contains code for pre-processing in-line phase-contrast CT projections for the 1st order TIE-based TEM DT reconstruction

#include <chrono>
#include <omp.h>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_tie.h"

using namespace xar;

#define PBITEMDT
#ifdef PBITEMDT

int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting PbiTEMDT program ...");

		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024], cparam5[1024];
		printf("\nReading PbiTEMDT.txt input parameter file ...");
		FILE* ff0 = fopen("PbiTEMDT.txt", "rt");
			if (!ff0) throw std::exception("Error: cannot open parameter file PbiTEMDT.txt.");

			fgets(cline, 1024, ff0); // 1st line - comment

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for input files
			if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) != 2) throw std::exception("Error reading input file name base from input parameter file.");
			string strfilepathin = cparam;
			printf("\nInput file name base = %s", strfilepathin.c_str());

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for output files
			if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
			string strfilepathout = cparam;
			printf("\nOutput file name base = %s", strfilepathout.c_str());

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of CT projection angles
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading the number of CT projection angles parameter from input parameter file.");
			int nangles = atoi(cparam); // needs to be 'int' for OpenMP 
			printf("\nTotal number of CT projection angles over 360 degrees = %d", nangles);
			int nangles2 = nangles / 2; 
			if (nangles <= 0 || nangles != nangles2 * 2) throw std::exception("Total number of CT projection angles over 360 degrees must be even.");
			double angle_step = 360.0 / nangles, angle; // the CT scan range must be 360 degrees!!!

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // parameters for trimming or padding the input image (to bring the dims to powers of 2)
			if (sscanf(cline, "%s %s %s %s %s %s", ctitle, cparam1, cparam2, cparam3, cparam4, cparam5) != 6) throw std::exception("Error reading image trim/pad parameters from input parameter file.");
			int iXLeft = atoi(cparam1);
			int iXRight = atoi(cparam2);
			int iYTop = atoi(cparam3);
			int iYBottom = atoi(cparam4);
			float fPadVal = (float)atof(cparam5);
			printf("\nImage trim(-)/pad(+) parameters: iXLeft = %d, iXRight = %d, iYTop = %d, iYBottom = %d, PadValue = %g", iXLeft, iXRight, iYTop, iYBottom, fPadVal);
			if (iXLeft > 0 && (iXRight < 0 || iYTop < 0 || iYBottom < 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
			else if (iXLeft < 0 && (iXRight > 0 || iYTop > 0 || iYBottom > 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
			else if (iXRight > 0 &&( iYTop < 0 || iYBottom < 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
			else if (iXRight < 0 && (iYTop > 0 || iYBottom > 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
			else if (iYTop > 0 && iYBottom < 0) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
			else if (iYTop < 0 && iYBottom > 0) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // wavelength in angstroms
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength parameter from input parameter file.");
			double wl = atof(cparam);
			printf("\nWavelength (A) = %g", wl);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // defocus distance in angstroms
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus distance parameter from input parameter file.");
			double defocus = atof(cparam);
			printf("\nDefocus distance (A) = %g", defocus);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // delta/beta
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading delta/beta parameter from input parameter file.");
			double delta2beta = atof(cparam);
			printf("\nDelta / beta = %g", delta2beta);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of worker threads
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of worker threads from input parameter file.");
			int nThreads = atoi(cparam);
			printf("\nNumber of worker threads = %d", nThreads);

			fclose(ff0);

		// create formatting string to add properly formatted indexes at the end of the output file names
		size_t i_dot = strfilepathin.rfind('.'), o_dot = strfilepathout.rfind('.'), nfieldB_length;
		char ndig[8];
		string myformat(""), myformat1("");
		if (nangles > 1)
		{
			nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the input file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
			nfieldB_length = 1 + size_t(log10(double(nangles2 - 1))); //maximum number of digits corresponding to angles in the output file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat1 += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
		}

		omp_set_num_threads(nThreads);
		#pragma omp parallel default(none)
		#pragma omp for schedule(dynamic) nowait
		for (int i = 0; i < nangles2; i++)
		{
			try // this is needed in OMP mode in order to catch exceptions within the worker threads
			{
				int j = i + nangles2;
				angle = angle_step * double(i);
				printf("\nAngle = %f", angle);

				// generate input and output file names
				char buffer[128], buffer2[128], bufferout[128];
				string infile_i, infile_i2; // input files with an in-line projection images
				string outfile_j; // output file with pre-processed image
				infile_i = strfilepathin;
				infile_i2 = strfilepathin;
				outfile_j = strfilepathout;
				sprintf(buffer, myformat.data(), i);
				sprintf(buffer2, myformat.data(), j);
				sprintf(bufferout, myformat1.data(), i);
				infile_i.insert(i_dot, buffer);
				infile_i2.insert(i_dot, buffer2);
				outfile_j.insert(o_dot, bufferout);

				// read the in-line projection image from input file
				XArray2D<float> xaint; // input image
				printf("\nReading input file = %s ...", infile_i.c_str());
				XArData::ReadFileGRD(xaint, infile_i.c_str(), wl);
				XArray2DMove<float> xamoveint(xaint);
				if (iXRight > 0 || iXLeft > 0 || iYTop > 0 || iYBottom > 0) xamoveint.Pad(index_t(iYTop), index_t(iYBottom), index_t(iXLeft), index_t(iXRight), fPadVal);
				else if (iXRight < 0 || iXLeft < 0 || iYTop < 0 || iYBottom < 0) xamoveint.Trim(index_t(-iYTop), index_t(-iYBottom), index_t(-iXLeft), index_t(-iXRight));
				//double l2dim1 = log2(xaint.GetDim1()), l2dim2 = log2(xaint.GetDim2());
				//if (int(l2dim1) != l2dim1 || int(l2dim2) != l2dim2) throw std::exception("Dimensions of the input image after pad/trim are not integer powers of 2.");

				// read the symmetric in-line projection image from input file
				XArray2D<float> xaint2; // input image2
				printf("\nReading input file = %s ...", infile_i2.c_str());
				XArData::ReadFileGRD(xaint2, infile_i2.c_str(), wl);
				XArray2DMove<float> xamoveint2(xaint2);
				if (iXRight > 0 || iXLeft > 0 || iYTop > 0 || iYBottom > 0) xamoveint2.Pad(index_t(iYTop), index_t(iYBottom), index_t(iXLeft), index_t(iXRight), fPadVal);
				else if (iXRight < 0 || iXLeft < 0 || iYTop < 0 || iYBottom < 0) xamoveint.Trim(index_t(-iYTop), index_t(-iYBottom), index_t(-iXLeft), index_t(-iXRight));
				//double l2dim1 = log2(xaint.GetDim1()), l2dim2 = log2(xaint.GetDim2());
				//if (int(l2dim1) != l2dim1 || int(l2dim2) != l2dim2) throw std::exception("Dimensions of the input image after pad/trim are not integer powers of 2.");

				// symmetrize the contrast function
				printf("\nSymmetrising the contrast function ...");
				xamoveint2.FlipX();
				float* arrI = &xaint.front();
				float* arrI2 = &xaint2.front();
				for (index_t k = 0; k < xaint.size(); k++)
				{
					arrI[k] = 0.5f * (arrI[k] + arrI2[k] - 2.0f);
				}

				// calculate regularized inverse Laplacian
				printf("\nCaclulating regulaized inverse Laplacian ...");
				XA_2DTIE<float> xatie;
				xatie.DP(xaint, delta2beta, defocus);

				// trim back and write the result to output file
				if (iXRight > 0 || iXLeft > 0 || iYTop > 0 || iYBottom > 0) xamoveint.Trim(index_t(iYTop), index_t(iYBottom), index_t(iXLeft), index_t(iXRight));
				printf("\nWriting output file = %s ...", outfile_j.c_str());
				XArData::WriteFileGRD(xaint, outfile_j.c_str(), eGRDBIN); 
			}
			catch (std::exception & E)
			{
				printf("\n\n!!!Exception: %s\n", E.what());
				exit(1);
			}
		}
	}
	catch (std::exception & E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}

	std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\n\nMain program finished. Execution time = %I64d s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	//printf("\nPress any key to exit..."); getchar();

	return 0;
}

#endif //PBITEMDT