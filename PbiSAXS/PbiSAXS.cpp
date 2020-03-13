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

		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024], cparam5[1024];
		printf("\nReading PbiSAXS.txt input parameter file ...");
		FILE* ff0 = fopen("PbiSAXS.txt", "rt");
			if (!ff0) throw std::exception("Error: cannot open parameter file PbiSAXS.txt.");

			fgets(cline, 1024, ff0); // 1st line - comment

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for input files
			if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) != 2) throw std::exception("Error reading input file name base from input parameter file.");
			string strfilepathin = cparam;
			printf("\nInput file name base = %s", strfilepathin.c_str());

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for output files
			if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
			string strfilepathout = cparam;
			printf("\nOutput file name base = %s", strfilepathout.c_str());

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of CT projection angles, stride
			if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading the number of CT projection angles & stride parameters from input parameter file.");
			int nangles = atoi(cparam); // needs to be 'int' for OpenMP 
			int istride = atoi(cparam1);
			printf("\nNumber of CT projection angles = %d", nangles);

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // CT angle range in degrees
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading CT angle range parameter from input parameter file.");
			double angle_range = atof(cparam);
			printf("\nCT angle range = %g", angle_range);
			double angle_step = angle_range / nangles, angle;

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

			fgets(cline, 1024, ff0); strtok(cline, "\n"); // calculate BornHom(0) or TIEHom+BornHom(1)
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading BornHom / Tie+BornHom / FresnelHom parameter from input parameter file.");
			int iMode = atoi(cparam);
			switch (iMode)
			{
				case 0:
					printf("\nThis program will calculate BornHom retrieval increment to TieHom.");
					break;
				case 1:
					printf("\nThis program will calculate Tie+BornHom retrieval.");
					break;
				case 2:
					printf("\nThis program will calculate FresnelHom retrieval.");
					break;
				default:
					throw std::exception("BornHom / Tie+BornHom / FresnelHom parameter values can only be 0, 1 or 2.");
			}
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of worker threads
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of worker threads from input parameter file.");
			int nThreads = atoi(cparam);
			printf("\nNumber of worker threads = %d", nThreads);

			fclose(ff0);

		string infile_i; // input file with an in-line projection image
		string outfile_j; // output file with SAXS image

		XArray2D<float> xaobjtie; // TIE-Hom retrieved intensity array
		XArray2D<float> xaint, xaint0; // auxillary real array
		XArray2D<fcomplex> xacamp; // auxillary complex array

		// create formatting string to add properly formatted indexes at the end of the output file names
		size_t i_dot = strfilepathin.rfind('.'), o_dot = strfilepathout.rfind('.'), nfieldB_length;
		char ndig[8];
		string myformat(""), myformat1("");
		if (nangles > 1)
		{
			nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the input file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
			nfieldB_length = 1 + size_t(log10(double(nangles / istride - 1))); //maximum number of digits corresponding to angles in the output file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat1 += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
		}

		omp_set_num_threads(nThreads);
		#pragma omp parallel default(none) private(xaobjtie, xaint, xaint0, xacamp, infile_i, outfile_j)
		#pragma omp for schedule(dynamic) nowait
		for (int i = 0; i < nangles; i += istride)
		{
			try // this is needed in OMP mode in order to catch exceptions within the worker threads
			{
				int j = i / istride;
				angle = angle_step * double(i);
				printf("\nAngle = %f", angle);

				// generate input and output file names
				char buffer[128], buffer1[128];
				infile_i = strfilepathin;
				outfile_j = strfilepathout;
				sprintf(buffer, myformat.data(), i);
				sprintf(buffer1, myformat1.data(), j);
				infile_i.insert(i_dot, buffer);
				outfile_j.insert(o_dot, buffer1);

				// read the in-line projection image from input file
				printf("\nReading input file = %s ...", infile_i.c_str());
				XArData::ReadFileGRD(xaint, infile_i.c_str(), wl);
				XArray2DMove<float> xamoveint(xaint);
				if (iXRight > 0 || iXLeft > 0 || iYTop > 0 || iYBottom > 0) xamoveint.Pad(index_t(iYTop), index_t(iYBottom), index_t(iXLeft), index_t(iXRight), fPadVal);
				else if (iXRight < 0 || iXLeft < 0 || iYTop < 0 || iYBottom < 0) xamoveint.Trim(index_t(-iYTop), index_t(-iYBottom), index_t(-iXLeft), index_t(-iXRight));
				double l2dim1 = log2(xaint.GetDim1()), l2dim2 = log2(xaint.GetDim2());
				if (int(l2dim1) != l2dim1 || int(l2dim2) != l2dim2) throw std::exception("Dimensions of the input image after pad/trim are not integer powers of 2.");

				xaint.Abs(); // @@@@@ temporary fix for bad input images with small negative values (but no zeros)
				if (xaint.Norm(xar::eNormMin) <= 0)
					throw std::exception("Input image intensity distribution contains some negative or zero values.");

				if (iMode == 0 || iMode == 1)
				{
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
					for (index_t k = 0; k < xacamp.size(); k++)
					{
						famp = pow(arrI[k], 0.5f);
						arrC[k] = std::polar<float>(famp, fdelta2beta * log(famp));
					}
					xacamp.SetHeadPtr(xaint.GetHeadPtr() ? xaint.GetHeadPtr()->Clone() : 0);
					XArray2DFFT<float> xafft2(xacamp);
					xafft2.Fresnel(defocus); // propagate forward to the image plane

					// Do 1st Born on the difference between the original image and DP-repropagated image
					float* arrI0 = &xaint0.front();
					for (index_t k = 0; k < xacamp.size(); k++) arrI[k] = arrI0[k] - std::norm(arrC[k]) + 1.0f;
					float fIin = (float)xaint.Norm(xar::eNormAver);
					XA_2DBorn<float> xaborn;
					xaborn.BornSC(xaint, defocus, delta2beta, alpha, false);
					arrI = &xaint.front();
					float* arrItie = &xaobjtie.front();
					//if (iMode == 0) // we output 1 - 2*Itie*Mu_Born - i.e. do nothing here
						//for (index_t i = 0; i < xaint.size(); i++) arrI[i] = 1.0f + (arrI[i] / fIin - 1.0f) / arrItie[i]; // to output just Mu_Born
					if (iMode == 1)	// we add TieHom to mu_Born here
						for (index_t k = 0; k < xaint.size(); k++) arrI[k] = arrItie[k] + arrI[k] - 1.0f;
				}

				if (iMode == 2)
				{
					XA_2DBorn<float> xaborn;
					xaborn.BornSC(xaint, defocus, delta2beta, alpha, true);
				}

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