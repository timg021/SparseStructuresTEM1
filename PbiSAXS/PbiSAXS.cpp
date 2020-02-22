//PbiSAXS.cpp this file contains code for extraction of SASX signal from in-line phase-contrast images

#include <chrono>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_tie.h"

using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting PbiSAXS program ...");

		string strfilepath = "C:\\Users\\gur017\\Downloads\\Temp\\";
		//string strfilepath = "C:\\Users\\tgureyev\\Downloads\\Temp\\"; // LysLesPhaseCT200keV_900rot\\";
		string infile = "img1.grd", infile_i; // input file with an in-line projection image
		string outfile = "saxs_" + infile, outfile_i; // output file with SAXS image

		index_t kkGS = 1; // number of GS refinement cycles to perform after TIE-Hom(DP) for each input image
		
		XArray2D<float> xaampin; // real amplitude of the input image
		XArray2D<float> xaobjtie; // TIE-Hom retrieved intensity array
		XArray2D<float> xaint; // auxillary real array
		XArray2D<fcomplex> xacamp; // auxillary complex array

		index_t nangles = 1; // 900;
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

		double wl = 0.0001; // 2.5e-6; // X-ray wavelength in microns
		double defocus = 1.e+6; // 0.015; // defocus distance in microns
		double delta2beta = 1000.0; // 0.1; // delta/beta

		for (index_t i = 0; i < nangles; i++)
		{
			angle = angle_step * double(i);
			printf("\nAngle = %f", angle);

			// generate input and output file names
			infile_i = infile;
			outfile_i = outfile;
			sprintf(buffer, myformat.data(), i);
			infile_i.insert(i_dot, buffer);
			outfile_i.insert(o_dot, buffer);
			printf("\nInput file = %s", infile_i.c_str());
			printf("\nOutput file = %s", outfile_i.c_str());

			// read the in-line projection image from input file
			XArData::ReadFileGRD(xaint, (strfilepath + infile_i).c_str(), wl);
			xaampin = xaint; 
			xaampin ^= 0.5; // save the input real amplitude for later use

			// do TIE-Hom phase retrieval
			XA_2DTIE<float> xatie;
			xatie.DP(xaint, delta2beta, defocus);
			xaobjtie = xaint; // save the TIE-Hom retrieved intensity for later use
			//XArData::WriteFileGRD(xaobjtie, (strfilepath + "TIE" + infile_i).c_str(), eGRDBIN);

			// Do Gerchberg-Saxton
			xaint ^= 0.5; // convert object-plane intensity into real amplitude
			MakeComplex(xaint, 0.0f, xacamp, true); // create complex amplitude with TIE-Hom retrieved real amplitude and zero phase
			XArray2DFFT<float> xafft2(xacamp);
			for (index_t kk = 0; kk < kkGS; kk++)
			{
				printf("\nDoing GS, k = %zd ...", kk);
				xatie.Homogenise(xacamp, delta2beta); // replace phase with the delta2beta*log(amplitude)
				xafft2.Fresnel(defocus); // propagate forward to the image plane
				xatie.ReplaceModulus(xacamp, xaampin); // replace the modulus with that of the input image
				xafft2.Fresnel(-defocus); // propagate back to the object plane
			}

			// write the result to output file
			Abs2(xacamp, xaint);
			XArData::WriteFileGRD(xaint, (strfilepath + "GS" + infile_i).c_str(), eGRDBIN);
			xaint -= xaobjtie;
			XArData::WriteFileGRD(xaint, (strfilepath + outfile_i).c_str(), eGRDBIN);
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