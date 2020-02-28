//PbiSAXS.cpp this file contains code for extraction of SASX signal from in-line phase-contrast images

#include <chrono>
#include <omp.h>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_tie.h"
#include "unwrap_2d_ljmu.h"

using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting PbiSAXS program ...");

		int nThreads = 1;
		//string strfilepathin = "C:\\Users\\tgureyev\\OneDrive - The University of Melbourne\\SAXS\\LysLesPhaseCT200keV_900rot\\";
		//string strfilepathout = "C:\\Users\\tgureyev\\Downloads\\Temp\\LysLesPhaseCT200keV_900rotSAXS\\";
		//string strfilepathin = "C:\\Users\\gur017\\OneDrive - The University of Melbourne\\SAXS\\LysLesPhaseCT200keV_900rot\\";
		//string strfilepathout = "C:\\Users\\gur017\\OneDrive - The University of Melbourne\\SAXS\\LysLesPhaseCT200keV_900rot\\";
		string strfilepathin = "C:\\Users\\gur017\\Downloads\\Temp\\";
		string strfilepathout = "C:\\Users\\gur017\\Downloads\\Temp\\";
		//string strfilepathin = "C:\\Users\\tgureyev\\Downloads\\Temp\\"; 
		//string strfilepathout = "C:\\Users\\tgureyev\\Downloads\\Temp\\";
		string infile = "img1m.grd", infile_i; // input file with an in-line projection image
		string outfile = "saxsnew_" + infile, outfile_i; // output file with SAXS image

		XArray2D<float> xaampin; // input real amplitude array
		XArray2D<float> xaobjtie; // TIE-Hom retrieved intensity array
		XArray2D<float> xaint; // auxillary real array
		XArray2D<fcomplex> xacamp; // auxillary complex array

		index_t numIter = 10; // number of GS iterations to perform after TIE-Hom(DP) for each input image
		index_t iYLeft = 256, iYRight = 256, iXLeft = 256, iXRight = 256;
		fcomplex tMaskVal = fcomplex(1.0f, 0.0f);
		double wl = 0.0001; // 2.5e-6; // X-ray wavelength in microns
		double defocus = 1.e+6; // 0.015; // defocus distance in microns
		double delta2beta = 100; // 1; // delta/beta

		int nangles = 1; // needs to be 'int' for OpenMP
		double angle_range = 180.0;
		double angle_step = angle_range / nangles, angle;

		char buffer[128];
		// create formatting string to add properly formatted indexes at the end of the output file names
		size_t i_dot = infile.rfind('.'), o_dot = outfile.rfind('.'), nfieldB_length;
		char ndig[8];
		string myformat("");
		if (nangles > 1)
		{
			nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the output file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
		}

		omp_set_num_threads(nThreads);
		#pragma omp parallel default(none) private(xaampin, xaobjtie, xaint, xacamp, infile_i, outfile_i, buffer)
		#pragma omp for schedule(dynamic) nowait
		for (int i = 0; i < nangles; i++)
		{
			angle = angle_step * double(i);
			printf("\nAngle = %f", angle);

			// generate input and output file names
			infile_i = infile;
			outfile_i = outfile;
			sprintf(buffer, myformat.data(), i);
			infile_i.insert(i_dot, buffer);
			outfile_i.insert(o_dot, buffer);

			// read the in-line projection image from input file
			printf("\nReading input file = %s ...", infile_i.c_str());
			XArData::ReadFileGRD(xaint, (strfilepathin + infile_i).c_str(), wl);
			//xaint += 1.0; //@@@@@@ temporary code for bad input data
			xaampin = xaint; 
			xaampin ^= 0.5; // convert input image into real amplitude for future use

			// do TIE-Hom phase retrieval
			XA_2DTIE<float> xatie;
			xatie.DP(xaint, delta2beta, defocus);
			xaobjtie = xaint; // save the TIE-Hom retrieved intensity for later use
			//XArData::WriteFileGRD(xaobjtie, (strfilepathout + "TIE" + infile_i).c_str(), eGRDBIN);

			// Do Gerchberg-Saxton (maybe worth trying HIO later!!)
			xaint ^= 0.5; // convert TIE-Hom retrieved object-plane intensity into real amplitude
			MakeComplex(xaint, 0.0f, xacamp, true); // create complex amplitude with TIE-Hom retrieved real amplitude and zero phase
			XArray2DFFT<float> xafft2(xacamp);
			for (index_t kk = 0; kk < numIter; kk++)
			{
				printf("\nDoing GS, k = %zd ...", kk);
				if (kk == 0) 
					xatie.Homogenise(xacamp, delta2beta); // replace the phase with delta2beta*log(amplitude)
				else
				{
					xatie.Homogenise1(xacamp, delta2beta); // replace the modulus with exp(beta2delta*phase)
					//xatie.EnforceSupport(xacamp, iYLeft, iYRight, iXLeft, iXRight, tMaskVal);
				}
				//XArData::WriteFileGRC(xacamp, (strfilepathout + "campTIE.grc").c_str(), eGRCBIN);
				xafft2.Fresnel(defocus); // propagate forward to the image plane
				//XArData::WriteFileGRC(xacamp, (strfilepathout + "camp1TIE.grc").c_str(), eGRCBIN);
				xatie.ReplaceModulus(xacamp, xaampin); // replace the modulus with that of the input image
				xafft2.Fresnel(-defocus); // propagate back to the object plane
				//XArData::WriteFileGRC(xacamp, (strfilepathout + "camp0GS1.grc").c_str(), eGRCBIN);
			}
			//Abs2(xacamp, xaint);
			CArg(xacamp, xaint); // it is very important to extract the result from the phase, rather than from intensity
			// Unwrap the phase
			XArray2D<double> xaintd(xaint.GetDim1(), xaint.GetDim2()), xaintd1(xaint.GetDim1(), xaint.GetDim2());;
			for (index_t i = 0; i < xaint.GetDim1(); i++)
				for (index_t j = 0; j < xaint.GetDim2(); j++)
					xaintd[i][j] = (double)xaint[i][j];
			XArray2D<char> input_mask(xaint.GetDim1(), xaint.GetDim2(), 255);
			XArray2DMove<char> xamove(input_mask);
			xamove.Mask(iYLeft, iYRight, iXLeft, iXRight, 0);
			unwrap2D(&xaintd.front(), &xaintd1.front(), (unsigned char*)&input_mask.front(), (int)xaintd.GetDim2(), (int)xaintd.GetDim1(), 0, 0, 0, 0);
			for (index_t i = 0; i < xaint.GetDim1(); i++)
				for (index_t j = 0; j < xaint.GetDim2(); j++)
					xaint[i][j] = (float)xaintd1[i][j];

			xaint /= float(0.5 * delta2beta);
			xaint.Exp();
			XArData::WriteFileGRD(xaint, (strfilepathout + "GSnew_" + infile_i).c_str(), eGRDBIN);
			xaint -= xaobjtie; // GS minus TIE_Hom
			xatie.DP(xaint, delta2beta / 10.0, defocus); // mysterious addional processing that seems to bring the result much closer to the "true SAXS" signal

			// write the result to output file
			//xaint += tMaskVal.real(); //@@@@@@ temporary code for bad input data
			printf("\nWriting output file = %s ...", outfile_i.c_str());
			XArData::WriteFileGRD(xaint, (strfilepathout + outfile_i).c_str(), eGRDBIN);
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