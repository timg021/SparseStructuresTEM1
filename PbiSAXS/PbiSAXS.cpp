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
		string infile = "LENA.grd"; // input file with an in-line projection image
		string outfile = "SAXS0" + infile; // output file with SAXS image
		XArray2D<float> xaimagein; // input in-line projection image array
		XArray2D<float> xaobjtie; // TIE-Hom retrieved intensity array
		XArray2D<float> xapha; // phase
		XArray2D<float> xaint; // intensity

		double wl = 1.e-4; // X-ray wavelength in microns
		double defocus = 1.e+6; // defocus distance in microns
		double delta2beta = 350; // delta/beta

		// read the in-line projection image from input file
		XArData::ReadFileGRD(xaimagein, (strfilepath + infile).c_str(), wl);
		xaint = xaimagein; // make a temporary copy of the input image

		// do TIE-Hom phase retrieval
		XA_2DTIE<float> xatie;
		xatie.DP(xaint, delta2beta, defocus);
		xaobjtie = xaint; // save the TIE-Hom retrieved intensity for later use
		//XArData::WriteFileGRD(xaobjtie, (strfilepath + "TIE" + infile).c_str(), eGRDBIN);

		// Fresnel-propagate the TIE-Hom retrieved object-plane complex amplitude forward to the image plane
		xapha = xaint;
		xapha.Log();
		xapha *= float(0.5 * delta2beta);
		xaint ^= 0.5; // convert intensity into real amplitude
		XArray2D<fcomplex> xacamp;
		MakeComplex(xaint, xapha, xacamp, true);
		XArray2DFFT<float> xafft2(xacamp);
		xafft2.Fresnel(defocus);
		Abs2(xacamp, xaint);

		// Subtract the retrieved-repropagated image from the original image, take the inverse Laplacian and divide by I0
		xaimagein -= xaint; // this is the SAXS component of the image
		//XArData::WriteFileGRD(xaimagein, (strfilepath + "SAXS1" + infile).c_str(), eGRDBIN);
		xatie.DP(xaimagein, delta2beta/100.0, defocus); // regularized inverse Laplacian
		//XArData::WriteFileGRD(xaimagein, (strfilepath + "SAXS0TIE" + infile).c_str(), eGRDBIN);
		xaimagein /= xaobjtie; // this is the SAXS component of the object

		// write the result to output file
		XArData::WriteFileGRD(xaimagein, (strfilepath + outfile).c_str(), eGRDBIN);
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