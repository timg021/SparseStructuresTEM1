//PbiSAXS.cpp this file contains code for extraction of SASX signal from in-line phase-contrast images

#include <chrono>
#include <omp.h>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_tie.h"
#include "XA_born.h"

using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting PbiSAXS program ...");

		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];
		printf("\nReading PbiSAXS.txt input parameter file ...");
		FILE* ff0 = fopen("PbiSAXS.txt", "rt");
			if (!ff0) throw std::exception("Error: cannot open parameter file PbiSAXS.txt.");

			fgets(cline, 1024, ff0); // 1st line - comment

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for input files
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading input file name base from input parameter file.");
			string strfilepathin = cparam;
			printf("\nInput file name base = %s", strfilepathin.c_str());

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for output files
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
			string strfilepathout = cparam;
			printf("\nInput file name base = %s", strfilepathout.c_str());

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of CT projection angles
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading the number of CT projection angles parameter from input parameter file.");
			int nangles = atoi(cparam); // needs to be 'int' for OpenMP 
			printf("\nNumber of CT projection angles = %d", nangles);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // CT angle range in degrees
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading CT angle range parameter from input parameter file.");
			double angle_range = atof(cparam);
			printf("\nCT angle range = %g", angle_range);
			double angle_step = angle_range / nangles, angle;

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // parameters for trimming the input image (e.g. to bring the dims to powers of 2)
			if (sscanf(cline, "%s %s %s %s %s", ctitle, cparam1, cparam2, cparam3, cparam4) != 5) throw std::exception("Error reading image trim parameters from input parameter file.");
			index_t iXLeft = (index_t)atoi(cparam1);
			index_t iXRight = (index_t)atoi(cparam2);
			index_t iYTop = (index_t)atoi(cparam3);
			index_t iYBottom = (index_t)atoi(cparam4);
			printf("\nImage trim parameters: iXLeft = %zd, iXRight = %zd, iYTop = %zd, iYBottom = %zd", iXLeft, iXRight, iYTop, iYBottom);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // wavelength in microns
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength parameter from input parameter file.");
			double wl = atof(cparam);
			printf("\nWavelength (microns) = %g", wl);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // defocus distance in microns
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus distance parameter from input parameter file.");
			double defocus = atof(cparam);
			printf("\nDefocus distance (microns) = %g", defocus);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // delta/beta
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading delta/beta parameter from input parameter file.");
			double delta2beta = atof(cparam);
			printf("\nDelta / beta = %g", delta2beta);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // regularization parameter for 1st Born
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading regularization parameter for 1st Born from input parameter file.");
			double alpha = atof(cparam);
			printf("\nRegularization parameter for 1st Born = %g", alpha);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of worker threads
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of worker threads from input parameter file.");
			int nThreads = atoi(cparam);
			printf("\nNumber of worker threads = %d", nThreads);

			fclose(ff0);

		string infile_i; // input file with an in-line projection image
		string outfile_i; // output file with SAXS image

		XArray2D<float> xaobjtie; // TIE-Hom retrieved intensity array
		XArray2D<float> xaint, xaint0; // auxillary real array
		XArray2D<fcomplex> xacamp; // auxillary complex array

		// create formatting string to add properly formatted indexes at the end of the output file names
		size_t i_dot = strfilepathin.rfind('.'), o_dot = strfilepathout.rfind('.'), nfieldB_length;
		char ndig[8];
		string myformat("");
		if (nangles > 1)
		{
			nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the output file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
		}

		omp_set_num_threads(nThreads);
		#pragma omp parallel default(none) private(xaobjtie, xaint, xacamp, infile_i, outfile_i)
		#pragma omp for schedule(dynamic) nowait
		for (int i = 0; i < nangles; i++)
		{
			try // this is needed in OMP mode in order to catch exceptions within the worker threads
			{
				angle = angle_step * double(i);
				printf("\nAngle = %f", angle);

				// generate input and output file names
				char buffer[128];
				infile_i = strfilepathin;
				outfile_i = strfilepathout;
				sprintf(buffer, myformat.data(), i);
				infile_i.insert(i_dot, buffer);
				outfile_i.insert(o_dot, buffer);

				// read the in-line projection image from input file
				printf("\nReading input file = %s ...", infile_i.c_str());
				XArData::ReadFileGRD(xaint, infile_i.c_str(), wl);
				XArray2DMove<float> xamoveint(xaint);
				xamoveint.Trim(iYTop, iYBottom, iXLeft, iXRight);
				xaint0 = xaint; // save the trimmed original image for later use

				// do TIE-Hom phase retrieval
				XA_2DTIE<float> xatie;
				xatie.DP(xaint, delta2beta, defocus);
				xaobjtie = xaint; // save the TIE-Hom retrieved intensity for later use

				// homogenise the object-plane complex amplitude and do forward propagation
				xacamp.Resize(xaint.GetDim1(), xaint.GetDim2());
				fcomplex* arrC = &xacamp.front();
				float* arrI = &xaint.front();
				float famp, fdelta2beta = float(delta2beta);
				for (index_t i = 0; i < xacamp.size(); i++)
				{
					famp = pow(arrI[i], 0.5f);
					arrC[i] = std::polar<float>(famp, fdelta2beta * log(famp));
				}
				xacamp.SetHeadPtr(xaint.GetHeadPtr() ? xaint.GetHeadPtr()->Clone() : 0);
				XArray2DFFT<float> xafft2(xacamp);
				xafft2.Fresnel(defocus); // propagate forward to the image plane

				// Do 1st Born on the difference between the original image and DP-repropagated image
				float* arrI0 = &xaint0.front();
				for (index_t i = 0; i < xacamp.size(); i++) arrI[i] = arrI0[i] - std::norm(arrC[i]) + 1.0f;
				float fIin = (float)xaint.Norm(xar::eNormAver);
				XA_2DBorn<float> xaborn;
				xaborn.BornSC(xaint, defocus, delta2beta, alpha);
				arrI = &xaint.front();
				float* arrItie = &xaobjtie.front();
				// we add 1 below in order to output not mu_Born, but 1 - mu_Born ~ exp(-mu_Born)
				for (index_t i = 0; i < xaint.size(); i++) arrI[i] = 1.0f + (arrI[i] / fIin - 1.0f) / (2.0f * arrItie[i]);

				// write the result to output file
				printf("\nWriting output file = %s ...", outfile_i.c_str());
				XArData::WriteFileGRD(xaint, outfile_i.c_str(), eGRDBIN);
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