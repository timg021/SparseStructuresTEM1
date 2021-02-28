// DHT.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <complex.h>
#include <chrono>
#include <omp.h>

#include "IXAHWave.h"
#include "XAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "fftwd3drc.h"
#include "XA_spln2.h"
#include "XA_spln3.h"
#include "XA_iwfr.h"
#include "XA_TIE.h"
#include "XA_tiff.h"

using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting Differential Holographic Tomography program ...");
		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;

		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024];
		FILE* ff0 = fopen("PhaseRetrieval.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PhaseRetrieval.txt.");

		// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
		while (true)
		{
			fgets(cline, 1024, ff0);
			if (!(cline[0] == '/' && cline[1] == '/')) break;
		}

		strtok(cline, "\n"); // 1. Input file with rotation angles and defocus distances
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading file name with rotation angles and defocus distances from input parameter file.");
		printf("\nReading defocus parameter file %s ...", cparam);
		ReadDefocusParamsFile(cparam, v2angles, vvdefocus, false);
		index_t nangles = v2angles.size(); // number of rotation steps 
		vector<index_t> vndefocus(nangles); // vector of numbers of defocus planes at different illumination angles
		for (index_t i = 0; i < nangles; i++) vndefocus[i] = vvdefocus[i].size();

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2. Input filename base of defocus series of the sample in TIFF, GRD or GRC format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		bool bTIFFinput;
		if (GetFileExtension(filenamebaseIn) == string(".TIFF") || GetFileExtension(filenamebaseIn) == string(".TIF")) bTIFFinput = true;
		else if (GetFileExtension(filenamebaseIn) == string(".GRD") || GetFileExtension(filenamebaseIn) == string(".GRC")) bTIFFinput = false;
		else throw std::exception("Error: input filename extension must be TIF, GRD or GRC.");
		if (GetFileExtension(filenamebaseIn) == string(".GRC")) // check that there is only one input defocused complex amplitude file per each illumination angle
		{
			for (index_t i = 0; i < nangles; i++)
				if (vndefocus[i] != 1)
					throw std::exception("Error: only one input defocused complex amplitude file per each illumination angle is allowed.");
		}
		//printf("\nInput defocus series file name base = %s", filenamebaseIn.c_str());  - this has to be delayed until the 3D Laplacian filter mode is known
		// FileNames2(vndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames - this has to be delayed until the 3D Laplacian filter mode is known

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3. Xmin Xmax Ymin Ymax in Angstroms (for TIFF input files only)
		double xlo, xhi, xst; // physical x-bounds in angstroms
		double ylo, yhi, yst; // physical y-bounds in angstroms
		if (bTIFFinput)
		{
			if (sscanf(cline, "%s %s %s %s %s", ctitle, cparam, cparam1, cparam2, cparam3) != 5) 
				throw std::exception("Error reading Xmin, Xmax, Ymin or Ymax from input parameter file.");
			xlo = atof(cparam);
			xhi = atof(cparam1);
			ylo = atof(cparam2);
			yhi = atof(cparam3);
			printf("\nPhysical boundaries of input images: Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g (Angstroms)", xlo, xhi, ylo, yhi);
		}
				
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in Angstroms
		printf("\nWavelength = %g (A)", wl);
		if (wl < 0 || wl > 1)
			throw std::exception("Error: wavelength value appears to be wrong.");
		double EE; // incident electron energy in volts (recalculated from the wavelength below)
		constexpr double hp = 6.62607004e-34; // Planck's constant (m2 kg / s)
		constexpr double cc = 299792458; // speed of light (m / s)
		constexpr double ee = 1.602176634e-19; // electron charge (coulomb)
		constexpr double m0 = 9.1093837015e-31; // electron rest mass (kg)
		constexpr long double mc2 = m0 * cc * cc; // mc^2
		long double chl2 = long double(cc * cc * hp * hp) / long double(wl * wl * 1.e-20);
		long double abra = sqrt(mc2 * mc2 + chl2);
		EE = double((-mc2 + abra) / long double(ee));
		printf("\nIncident electron energy E = %g (keV)", EE / 1000.0);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam); 
		if (aobj != 0) printf("\nObjective aperture = %g (mrad)", aobj);
		else  printf("\nObjective aperture is infinite.");
		if (aobj < 0 || aobj > 1000)
			throw std::exception("Error: objective aperture value appears to be wrong.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6. Spherical aberrations Cs3 and Cs5 in mm
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading spherical aberrations from input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> Angstroms
		Cs5 *= 1.e+7; // mm --> Angstroms

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7.Phase retrieval method: 1 = IWFR, 2 = CTFL2, 3 = MinLogAmp, 4 = PhaseB7
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading phase retrieval method from input parameter file.");
		int nPhaseRetrieval = atoi(cparam);
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Save or not phase phase-retrieved defocused complex amplitudes in files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading save_or_not phase-retrieved defocused complex amplitudes switch from input parameter file.");
		int nSaveDefocCAmps = atoi(cparam);
		if (nSaveDefocCAmps == 1)
			printf("\nPhase-retrieved defocused complex amplitudes will be saved in GRC files");
		else if (nSaveDefocCAmps == 0) 
				printf("\nPhase-retrieved defocused complex amplitudes will not be saved in files");
			else
				throw std::exception("Error: save_or_not phase-retrieved defocused complex amplitudes switch must be 0 or 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9. Maximal number of IWFR iterations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading maximal number of iterations from input parameter file.");
		int kmax = atoi(cparam);
		printf("\nMaximal number of iterations = %d", kmax);
		if (kmax < 1)
			throw std::exception("Error: the maximal number of iterations should be >= 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10. Minimal phase reconstruction error(IWFR) or Tikhonov regularization parameter alpha(CTFL2)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading minimal phase reconstruction error/regulariztion parameter from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nMinimal phase reconstruciton error (IWFR) or Tikhonov regularization parameter (CTFL2) = %g", epsilon);
		if (epsilon < 0)
			throw std::exception("Error: minimal phase reconstruction error must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11. Output defocus distances min max and step in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading output defocus distances from input parameter file.");
		double zlo = atof(cparam); // minimum output defocus in Angstroms - !!! will be corrected with dzextra below
		double zhi = atof(cparam1); // maximum output defocus in Angstroms - !!! will be corrected with dzextra below 
		double zst = abs(atof(cparam2)); // output defocus step in Angstroms 
		if (zlo > zhi) std::swap(zlo, zhi);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //12. Extra defocus for 3D reconstruction in_Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading extra defocus for 3D reconstruction parameter from input parameter file.");
		double dzextra = atof(cparam);
		printf("\nExtra defocus for 3D reconstruction = %g (Angstroms)", dzextra);
		double zlodz(zlo), zhidz(zhi);
		zlodz += dzextra; zhidz += dzextra;

		printf("\nOutput defocus distances: min = %g, max = %g, step = %g (Angstroms)", zlo, zhi, zst);
		if (zst == 0)
			throw std::exception("Error: output defocus step cannot be zero.");
		int noutdefocus = int((zhidz - zlodz) / zst + 0.5); // number of defocus planes to propagate to
		if (noutdefocus <= 0)
			throw std::exception("Error: number of output defocus planes must be positive.");
		vector<double> voutdefocus(noutdefocus); // vector of output defocus distances
		printf("\nThere are %d output defocus plane positions", noutdefocus);
		for (index_t n = 0; n < noutdefocus; n++)
		{
			voutdefocus[n] = zlodz + zst * n;
			//printf("%g ", voutdefocus[n]);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13. Low threshold value in Volts and multiplicative scaling factor (dimensionless) for 3D potential
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading low threshold value and scale factor for 3D_potential from input parameter file.");
		double dLowThreshold = atof(cparam);
		double dScalingFactor = atof(cparam1);
		printf("\nLow threshold value for 3D_potential = %g (V)", dLowThreshold);
		printf("\nMultiplicative scaling factor for 3D potential = %g (dimensionless)", dScalingFactor);
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14. 3D Laplacian filter mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading 3D Laplacian filter mode from input parameter file.");
		int imode3Dfilter = atoi(cparam);
		switch (imode3Dfilter)
		{
		case 0:
			printf("\n3D Laplacian filter won't be applied.");
			break;
		case 1:
			printf("\n3D Laplacian filter will be applied.");
			break;
		case 2:
			printf("\nThe program will only apply the 3D Laplacian filter to the input files.");
			break;
		default:
			throw std::exception("Error: unknown value for 3D Laplacian filter mode in input parameter file.");
		}
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15. Regularization parameter for inverse 3D Laplacian
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading regularization parameter for inverse 3D Laplacian from input parameter file.");
		double alpha = atof(cparam);
		printf("\nRegularization parameter for inverse 3D Laplacian = %g", alpha);
		if (alpha < 0 && (imode3Dfilter == 1 || imode3Dfilter == 2))
			throw std::exception("Error: regularization parameter for 3D Laplacian filter must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 16. Multislice reprojection of the 3D electrostatic potential mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading reprojectin of 3D potential mode from input parameter file.");
		int imodeReproj = atoi(cparam);
		switch (imodeReproj)
		{
		case 0:
			printf("\nMultislice reprojection of the 3D electrostatic potential won't be applied.");
			break;
		case 1:
			printf("\nMultislice reprojection of the 3D electrostatic potential will be applied.");
			break;
		case 2:
			printf("\nThe program will only apply multislice reprojection of the 3D electrostatic potential imported from the input files.");
			break;
		default:
			throw std::exception("Error: unknown value for 3D Laplacian filter mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17. Slice thickness for multislice reprojection in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading slice thickness for multislice reprojection from input parameter file.");
		double sliceTh = atof(cparam);
		printf("\nSlice thickness for multislice reprojection = %g (Angstroms)", sliceTh);
		if (sliceTh < (4.0 * zst) && (imodeReproj == 1 || imodeReproj == 2))
			throw std::exception("Error: slice thickness for multislice reprojection must be 4 x z_step or larger.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 18. Output file name base in GRD or TIFF format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		bool bTIFFoutput;
		if (GetFileExtension(filenamebaseOut) == string(".TIFF") || GetFileExtension(filenamebaseOut) == string(".TIF")) bTIFFoutput = true;
		else if (GetFileExtension(filenamebaseOut) == string(".GRD")) bTIFFoutput = false;
		else throw std::exception("Error: output filename extension must be TIF ot GRD.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19. Verbose output during execution? Yes = 1, No = 0
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading verbose output parameter from input parameter file.");
		bool bVerboseOutput(true); // if this is TRUE, additional information is printed during execution
		(atoi(cparam) == 0 || atoi(cparam) == 1) ? bVerboseOutput = (bool)atoi(cparam) : throw std::exception("Error: verbose output parameter must be 0 or 1 in input parameter file.");
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 20. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		//************************************ create vectors of input and output file names
		
		vector<string> vinfilenamesTot, vinfilenames1Tot; // input filenames for the defocused intensities or complex amplitudes
		vector< vector<string> > vvinfilenames(nangles), vvinfilenames1(nangles); // same input filenames for the defocused series in the form of vector of vectors
		if (imode3Dfilter == 2 || imodeReproj == 2) // only read in pre-caculated 3D potential from "output" files and filter or reproject it
		{
			printf("\nInput file name base for pre-existing 3D potential = %s", filenamebaseOut.c_str());
			FileNames(1, noutdefocus, filenamebaseOut, vinfilenamesTot); // create 1D array of input filenames to read the input 2D slices of a previously reconstructed 3D object
		}
		else // read in the input files as directed
		{
			printf("\nInput defocus series file name base = %s", filenamebaseIn.c_str());
			FileNames2(vndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames
		}
		index_t ndefcurrent = 0;
		for (index_t na = 0; na < nangles; na++)
		{
			vvinfilenames[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
			for (index_t n = 0; n < vndefocus[na]; n++) vvinfilenames[na][n] = vinfilenamesTot[ndefcurrent++];
		}

		string filenamebaseOutNew("BAD_STRING"), filenamebaseOutCAmp("BAD_STRING"), filenamebaseOutDefocCAmp("BAD_STRING"); // don't use it, unless it is redefined later
		vector<string> voutfilenamesTot; // output filenames for the reconstructed 3D potential
		if (imode3Dfilter == 2) // output the Laplace-filtered 3D potential
		{
			filenamebaseOutNew = filenamebaseOut;
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "L");
			printf("\nOutput file name base for the filtered 3D potential = %s", filenamebaseOutNew.c_str());
			FileNames(1, noutdefocus, filenamebaseOutNew, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the filtered 3D object
		}
		else if (imodeReproj != 2) // in the case imodeReproj == 2, there is no point to save the same 3D potential as read from the input files
		{
			printf("\nOutput file name base for 3D potential = %s", filenamebaseOut.c_str());
			FileNames(1, noutdefocus, filenamebaseOut, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
		}

		vector<string> voutfilenamesTotCAmp; // output filenames the reprojected defocused complex amplitudes
		vector< vector<string> > vvoutfilenamesCAmp(nangles);  // same output filenames for the reprojected defocused complex amplitudes in the form of vector of vectors
		if (imodeReproj == 1 || imodeReproj == 2) // create filenames for saving the reprojected intensities, based on the names of input defocused intensities
		{
			filenamebaseOutCAmp = filenamebaseIn;
			filenamebaseOutCAmp.replace(filenamebaseOutCAmp.find_last_of("."), filenamebaseOutCAmp.length() - filenamebaseOutCAmp.find_last_of("."), "R.grc");
			printf("\nOutput file name base for reprojected complex amplitudes = %s", filenamebaseOutCAmp.c_str());
			FileNames2(vndefocus, filenamebaseOutCAmp, voutfilenamesTotCAmp); // create "2D array" of output filenames for reprojected complex amplitudes
			
			// Also copy the filenames of the input defocused intensities for the purpose of comparing with the reprojected intensities
			// note that in the case imode3Dfilter == 2 or imodeReproj == 2, these files are different from input images (which come from the potential)
			printf("\nDefocus series file name base for comparing with re-projected images = %s", filenamebaseIn.c_str());
			FileNames2(vndefocus, filenamebaseIn, vinfilenames1Tot); // create "total 2D array" of input filenames

			// Create the vector of vectors of output filenames for reprojected complex amplitudes at different defocus distances at each illumination angle
			// and the vector of vectors of filenames of the input defocused images for the purpose of comparing with the reprojected intensities
			index_t ndefcurrent = 0;
			for (index_t na = 0; na < nangles; na++)
			{
				vvoutfilenamesCAmp[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				vvinfilenames1[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				for (index_t n = 0; n < vndefocus[na]; n++)
				{
					vvoutfilenamesCAmp[na][n] = voutfilenamesTotCAmp[ndefcurrent];
					vvinfilenames1[na][n] = vinfilenames1Tot[ndefcurrent];
					ndefcurrent++;
				}
			}
		}

		vector<string> voutfilenamesTotDefocCAmp; // output filenames the reprojected defocused complex amplitudes
		if (nSaveDefocCAmps == 1) // create filenames for saving the phase-retrieved defocused complex amplitudes
		{
			filenamebaseOutDefocCAmp = filenamebaseIn;
			filenamebaseOutDefocCAmp.replace(filenamebaseOutDefocCAmp.find_last_of("."), filenamebaseOutDefocCAmp.length() - filenamebaseOutDefocCAmp.find_last_of("."), "D.grc");
			printf("\nOutput file name base for phase-retrieved defocused complex amplitudes = %s", filenamebaseOutDefocCAmp.c_str());
			FileNames(nangles, 1, filenamebaseOutDefocCAmp, voutfilenamesTotDefocCAmp);

		}

		//************************************ end creating vectors of input and output file names


		//*********************************** start main calculations
		bool bAbort(false);
		index_t nx, ny; // transverse array sizes in pixels

		std::unique_ptr<IXAHead> pHead(nullptr);
		XArray3D<double> K3out; // big 3D reconstructed array (needs to fit into RAM alongside with with everything else)

		if (imode3Dfilter != 2 && imodeReproj != 2) // do phase retrieval and backpropagation prior to 3D filtering or reprojection and output
		{
			printf("\n\nPerforming phase retrieval and backpropagation into the 3D volume ...");
			// start of cycle over illumination angles
			//index_t ndefcurrent(0);
			for (index_t na = 0; na < nangles; na++)
			{
				double angleY = v2angles[na].a * PI180;
				double angleX = v2angles[na].b * PI180;
				double cosangleY = cos(angleY);
				double sinangleY = sin(angleY);
				double cosangleX = cos(angleX);
				double sinangleX = sin(angleX);
				if (bVerboseOutput) printf("\n\n*** Illumination angle[%zd] = (%g, %g) (degrees)", na, angleY / PI180, angleX / PI180);
				else printf("\n*** Illumination angle[%zd] = (%g, %g) (degrees)", na, angleY / PI180, angleX / PI180);

				index_t ndefocus = vndefocus[na]; // number of defocus planes at the current illumination angle
				vector<Pair> vdefocus = vvdefocus[na]; // vector of input Z" rotation angles and defocus positions at the current illumination angle
				vector<string> vinfilenames = vvinfilenames[na]; // input filenames of defocused images at the current illumination angle
				vector<string> voutfilenames(noutdefocus);

				// define main work objects
				double zout(0.0); // z plane at which campOut will be defined
				XArray2D<dcomplex> campOut; // complex amplitude in the "middle" defocus plane

				if (GetFileExtension(vinfilenames[0]) == string(".GRC")) // phase retrieval is not required - work with the provided complex amplitude
				{
					if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[0].a, vdefocus[0].b);
					zout = vdefocus[0].b;
					if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[0].c_str());
					XArData::ReadFileGRC(campOut, vinfilenames[0].c_str()); //	read the only one (!) input GRC file
					pHead.reset(campOut.GetHeadPtr()->Clone());
					IXAHWave2D* ph2 = GetIXAHWave2D(campOut);
					ny = campOut.GetDim1();
					nx = campOut.GetDim2();
					xlo = ph2->GetXlo();
					xhi = ph2->GetXhi();
					xst = ph2->GetXStep(nx);
					ylo = ph2->GetYlo();
					yhi = ph2->GetYhi();
					yst = ph2->GetYStep(ny);
					if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the 3D reconstructed object is supposed to have cubic voxels");
					if (vdefocus[0].a != 0) // rotate around Z" if needed
					{
						// rotate input defocused complex amplitude around Z" back to zero angle
						XArray2D<double> ampRe0, ampIm0, ampRe, ampIm;
						Re(campOut, ampRe0); Im(campOut, ampIm0);
						XArray2DSpln<double> xaSplnRe(ampRe0), xaSplnIm(ampIm0);
						xaSplnRe.Rotate(ampRe, -vdefocus[0].a, 0.5 * (yhi + ylo), 0.5 * (xhi + xlo), 1.0f); // expecting uniform background with unit amplitude
						xaSplnIm.Rotate(ampIm, -vdefocus[0].a, 0.5 * (yhi + ylo), 0.5 * (xhi + xlo), 0.0f); // expecting uniform background with unit amplitude
						MakeComplex(ampRe, ampIm, campOut, false);
					}
				}
				else // do phase retrieval and find the complex wave amplitude in the plane z = zout
				{
					if (ndefocus == 1)
					{
						if (nPhaseRetrieval != 3 && nPhaseRetrieval != 4)
							throw std::exception("Error: unsuitable phase retrieval method in the input parameter file for the case of single defocus distance.");
					}
					else
					{
						if (nPhaseRetrieval != 1 && nPhaseRetrieval != 2)
							throw std::exception("Error: unsuitable phase retrieval method in the input parameter file for the case of multiple defocus distances.");
					}

					// read defocused images from files
					vector<XArray2D<double>> vint0(ndefocus); // initial defocused intensities
					for (int n = 0; n < ndefocus; n++)
					{
						if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[n].a, vdefocus[n].b);
						if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[n].c_str());
						if (bTIFFinput)
						{
							TIFFReadFile(vint0[n], vinfilenames[n].c_str()); // read input TIFF files
							if (n == 0)
							{
								pHead.reset(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
								ny = vint0[0].GetDim1();
								nx = vint0[0].GetDim2();
								xst = (nx > 1) ? (xhi - xlo) / double(nx) : 1.0;
								yst = (ny > 1) ? (yhi - ylo) / double(ny) : 1.0;
							}
							else
							{
								ny = vint0[n].GetDim1();
								nx = vint0[n].GetDim2();
								if (ny != vint0[0].GetDim1() || nx != vint0[0].GetDim2())
									throw std::exception("Error: input TIFF files have different dimensions");
							}
							vint0[n].SetHeadPtr(pHead->Clone());
						}
						else
						{
							XArData::ReadFileGRD(vint0[n], vinfilenames[n].c_str(), wl); //	read input GRD files
							if (n == 0) pHead.reset(vint0[0].GetHeadPtr()->Clone());
							IXAHWave2D* ph2 = GetIXAHWave2D(vint0[n]);
							ny = vint0[n].GetDim1();
							nx = vint0[n].GetDim2();
							xlo = ph2->GetXlo();
							xhi = ph2->GetXhi();
							xst = ph2->GetXStep(nx);
							ylo = ph2->GetYlo();
							yhi = ph2->GetYhi();
							yst = ph2->GetYStep(ny);
							if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the 3D reconstructed object is supposed to have cubic voxels");
						}
						// rotate input defocused image around Z" back to zero angle
						if (vdefocus[n].a != 0) // rotate around Z" if needed
						{
							XArray2D<double> vintTemp(vint0[n]);
							XArray2DSpln<double> xaSpln(vintTemp);
							xaSpln.Rotate(vint0[n], -vdefocus[n].a, 0.5 * (yhi + ylo), 0.5 * (xhi + xlo), 1.0f); // expecting uniform background with unit amplitude
						}
					}

					// do phase retrieval
					switch (nPhaseRetrieval)
					{
					case 1:
					{
						vector<double> vdefocusdist(ndefocus);
						for (int n = 0; n < ndefocus; n++) vdefocusdist[n] = vdefocus[n].b;
						if (bVerboseOutput) printf("\nPerforming IWFR phase retrieval ...");
						zout = 0.0;
						for (size_t n = 0; n < ndefocus; n++) zout += vdefocus[n].b;
						zout /= double(ndefocus);
						XA_IWFR<double> xa_iwfr;
						xa_iwfr.Iwfr(vint0, campOut, vdefocusdist, zout, k2maxo, Cs3, Cs5, kmax, epsilon, bVerboseOutput);
					}
					break;
					case 2:
					{
						vector<double> vdefocusdist(ndefocus);
						for (int n = 0; n < ndefocus; n++) vdefocusdist[n] = vdefocus[n].b;
						if (bVerboseOutput) printf("\nPerforming CTFL2 phase retrieval ...");
						zout = vdefocus[0].b; //CTFL2 always retrieves the phase in the first defocused plane
						XA_IWFR<double> xa_iwfr;
						XArray2D<double> fiOut, ampOut;
						ampOut = vint0[0]; // use this before the call to xa_iwfr.CTFL2 spoils all vint0 images
						ampOut ^= 0.5; // preserve the amplitude at the zeros distance for reuse after phase retrieval
						double alphaCTF2 = epsilon; // alphaCTF2 will be changed by the call to CTFL2()
						xa_iwfr.CTFL2(vint0, fiOut, vdefocusdist, k2maxo, Cs3, Cs5, alphaCTF2);
						if (bVerboseOutput) printf("\nMinimal [sum(CTF_n^2)]^2 in denominator = %g", alphaCTF2);
						MakeComplex(ampOut, fiOut, campOut, true);
					}
					break;
					case 3:
					{
						zout = vdefocus[0].b;
						XArray2D<double> fiOut, ampOut;
						ampOut = vint0[0]; ampOut ^= 0.5;
						XA_IWFR<double> xa_iwfr;
						if (bVerboseOutput) printf("\nPerforming -0.5*Log(I) phase retrieval ...");
						xa_iwfr.MinLogAmp(vint0[0], fiOut);
						MakeComplex(ampOut, fiOut, campOut, true);
					}
					break;
					case 4:
					{
						zout = vdefocus[0].b;
						XArray2D<double> fiOut, ampOut;
						ampOut = vint0[0]; ampOut ^= 0.5;
						XA_IWFR<double> xa_iwfr;
						if (bVerboseOutput) printf("\nPerforming -0.5*sqrt(Kmax^2-K^2) phase retrieval ...");
						xa_iwfr.PhaseB7(vint0[0], fiOut);
						MakeComplex(ampOut, fiOut, campOut, true);
					}
					break;
					default:
						throw std::exception("Error: unknown phase retrieval method in the input parameter file.");
					}

					// optionally save phase retrieved defocused complex amplitudes in GRC files
					if (nSaveDefocCAmps == 1)
						XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmp[na].c_str(), eGRCBIN);
				}

				// now start calculations of the output defocused images
				if (bVerboseOutput) printf("\nPropagating to output defocus planes ...");
				vector<XArray2D<dcomplex>> vcamp(noutdefocus); // defocused complex amplitudes
				#pragma omp parallel for
				for (int n = 0; n < noutdefocus; n++)
				{
					try
					{
						// propagate to the output defocused plane
						vcamp[n] = campOut;
						xar::XArray2DFFT<double> xafft(vcamp[n]);
						xafft.Fresnel(voutdefocus[n] - zout, false, k2maxo, Cs3, Cs5); // propagate to z_out[n]
					}
					catch (std::exception& E)
					{
						printf("\n\n!!!Exception: %s\n", E.what());
						bAbort = true;
					}
				}
				if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

				if (bVerboseOutput) printf("\nUpdating 3D reconstructed object ...");
				double xc = (xhi + xlo) / 2.0; // x-coordinate of the centre of rotation
				double yc = (yhi + ylo) / 2.0; // y-coordinate of the centre of rotation
				double zc = (zhidz + zlodz) / 2.0; // z-coordinate of the centre of rotation

				// calculate the coordinate illumination angle parameters
				double xxx;
				vector<double> x_sinangleY(nx), x_cosangleY(nx);
				for (index_t i = 0; i < nx; i++)
				{
					xxx = xlo + xst * i - xc;
					x_sinangleY[i] = xxx * sinangleY;
					x_cosangleY[i] = xxx * cosangleY;
				}
				double yyy;
				vector<double> y_sinangleX(ny), y_cosangleX(ny);
				for (index_t j = 0; j < ny; j++)
				{
					yyy = ylo + yst * j - yc;
					y_sinangleX[j] = yyy * sinangleX;
					y_cosangleX[j] = yyy * cosangleX;
				}
				double zzz;
				vector<double> z_sinangleX(noutdefocus), z_cosangleX(noutdefocus);
				for (index_t n = 0; n < noutdefocus; n++)
				{
					zzz = voutdefocus[n] - zc;
					z_sinangleX[n] = zzz * sinangleX;
					z_cosangleX[n] = zzz * cosangleX;
				}

				// allocate the large 3D output array
				if (na == 0)
				{
					K3out.Resize(noutdefocus, ny, nx, 0.0);
					K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
				}

				// rotate the defocus plane around the "vertical" y axis by -angle, instead of rotating the 3D object by the angle
				int ii, jj, nn;
				double zz;
				int nx2 = int(nx) - 2, ny2 = int(ny) - 2, noutdefocus2 = noutdefocus - 2;
				double zz_sinangleY, zz_cosangleY, dK, dx0, dx1, dy0, dy1, dz0, dz1, dz0K, dz1K;

				#pragma omp parallel for shared (K3out, vcamp) private(yyy, zz, zz_sinangleY, zz_cosangleY, xxx, zzz, dx1, dx0, dy1, dy0, dz1, dz0, ii, jj, nn, dK, dz0K, dz1K)
				for (int n = 0; n < noutdefocus; n++)
				{
					try
					{
						for (index_t j = 0; j < ny; j++)
						{
							// inverse rotation around X' axis
							yyy = yc + y_cosangleX[j] - z_sinangleX[n]; // y coordinate with respect to the rotated 3D sample
							dy1 = abs(yyy - ylo) / yst;
							jj = (int)dy1; if (jj < 0 || jj > ny2) continue;
							dy1 -= jj; dy0 = 1.0 - dy1;
							zz = y_sinangleX[j] + z_cosangleX[n];
							zz_sinangleY = xc - zz * sinangleY;
							zz_cosangleY = zc + zz * cosangleY;
								
							for (index_t i = 0; i < nx; i++)
							{
								// inverse rotation around Y axis
								xxx = x_cosangleY[i] + zz_sinangleY; // x coordinate with respect to the rotated 3D sample
								zzz = x_sinangleY[i] + zz_cosangleY; // z coordinate with respect to the rotated 3D sample
								dx1 = abs(xxx - xlo) / xst;
								ii = (int)dx1; if (ii < 0 || ii > nx2) continue;
								dx1 -= ii; dx0 = 1.0 - dx1;
								dz1 = abs(zzz - zlodz) / zst;
								nn = (int)dz1; if (nn < 0 || nn > noutdefocus2) continue;
								dz1 -= nn; dz0 = 1.0 - dz1;
								dK = 1.0 - std::norm(vcamp[n][j][i]);
								dz0K = dz0 * dK;
								dz1K = dz1 * dK;
								K3out[nn][jj][ii] += dx0 * dy0 * dz0K;
								K3out[nn][jj][ii + 1] += dx1 * dy0 * dz0K;
								K3out[nn + 1][jj][ii] += dx0 * dy0 * dz1K;
								K3out[nn + 1][jj][ii + 1] += dx1 * dy0 * dz1K;
								K3out[nn][jj + 1][ii] += dx0 * dy1 * dz0K;
								K3out[nn][jj + 1][ii + 1] += dx1 * dy1 * dz0K;
								K3out[nn + 1][jj + 1][ii] += dx0 * dy1 * dz1K;
								K3out[nn + 1][jj + 1][ii + 1] += dx1 * dy1 * dz1K;
							}
						}
					}
					catch (std::exception& E)
					{
						printf("\n\n!!!Exception: %s\n", E.what());
						bAbort = true;
					}
				}
				if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
			} // end of cycle over illumination angles

			// renormalization of the reconstructed 3D distribution
			double ww = 1.0; // "average width" of the potential distribution around a single atom in A, this may change in the future
			double drenorm = -EE / ww * xst; // @@@@@ in principle, it is supposed to be more like xst * yst instead of just xst, and also divided by dzextra
			drenorm *= dScalingFactor;
			//drenorm = -6000.0; // @@@@@ numerical tests seem to indicate that this is an "optimal" value for this factor
			//printf("\n\nRenormalization parameter for the reconstructed 3D distribution = %g", drenorm);
			drenorm /= double(nangles);
			K3out.ThresholdLow(dLowThreshold, 0.0, drenorm);

		} // end of case if imode3Dfilter != 2
		else // read in pre-existing 3D distribution of the electrostatic potential
		{
			printf("\n\nReading pre-existing 3D distribution of the electrostatic potential from files ...");
			// read input 2D slices
			XArray2D<double> ipIn;
			for (int n = 0; n < noutdefocus; n++)
			{
				if (bTIFFinput)
				{
					TIFFReadFile(ipIn, vinfilenamesTot[n].c_str()); // read input TIFF files
					if (n == 0)
					{
						pHead.reset(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
						ny = ipIn.GetDim1();
						nx = ipIn.GetDim2();
						xst = (nx > 1) ? (xhi - xlo) / double(nx) : 1.0;
						yst = (ny > 1) ? (yhi - ylo) / double(ny) : 1.0;
						if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the input 3D object is supposed to have cubic voxels");
						K3out.Resize(noutdefocus, ny, nx, 0.0);
						K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
					}
					else
					{
						if (ny != ipIn.GetDim1() || nx != ipIn.GetDim2())
							throw std::exception("Error: input TIFF files have different dimensions");
					}
				}
				else
				{

					XArData::ReadFileGRD(ipIn, vinfilenamesTot[n].c_str(), wl); //	read input GRD files
					if (n == 0)
					{
						IXAHWave2D* ph2 = GetIXAHWave2D(ipIn);
						ny = ipIn.GetDim1();
						nx = ipIn.GetDim2();
						xlo = ph2->GetXlo();
						xhi = ph2->GetXhi();
						xst = ph2->GetXStep(nx);
						ylo = ph2->GetYlo();
						yhi = ph2->GetYhi();
						yst = ph2->GetYStep(ny);
						if (xst != yst || xst != zst)	throw std::exception("Error: the input 3D object is supposed to have cubic voxels");
						pHead.reset(ipIn.GetHeadPtr()->Clone());
						K3out.Resize(noutdefocus, ny, nx, 0.0);
						K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
					}
				}
				for (index_t j = 0; j < ny; j++)
					for (index_t i = 0; i < nx; i++)
						K3out[n][j][i] = ipIn[j][i];
			}
		}
		 
		// regularized inverse 3D Laplacian filter
		if (imode3Dfilter == 1 || imode3Dfilter == 2)
		{
			printf("\n\nInverse Laplace filtering 3D reconstructed object ...");
			double dnorm1 = K3out.Norm(eNormL1);
				
			//allocate space for FFT transform and create FFTW plans
			Fftwd3drc fftf((int)noutdefocus, (int)ny, (int)nx);

			// FFT of test array
			fftf.SetRealXArray3D(K3out);
			fftf.ForwardFFT();

			/// multiply FFT of K3out arrays by the FFT version of regularized inverse 3D Laplacian
			double fact = 4.0 * PI * PI * wl * (abs(dzextra) * 2.0);
			double dksi2 = fact / ((xhi - xlo) * (xhi - xlo));
			double deta2 = fact / ((yhi - ylo) * (yhi - ylo));
			double dzeta2 = fact / ((zhi - zlo) * (zhi - zlo));
			double dtemp, dtemp1, dk2, djk2;

			fftw_complex* pout = fftf.GetComplex();
			int m = 0, nc2 = fftf.GetNx2();
			index_t k1, j1, nyd2 = ny / 2, nzd2 = noutdefocus / 2;
			for (index_t k = 0; k < noutdefocus; k++)
			{
				k <= nzd2 ? k1 = k : k1 = noutdefocus - k;
				dk2 = k1 * k1 * dzeta2 + alpha;
				for (index_t j = 0; j < ny; j++)
				{
					j <= nyd2 ? j1 = j : j1 = ny - j;
					djk2 = j1 * j1 * deta2 + dk2;
					for (index_t i = 0; i < nc2; i++)
					{
						dtemp = i * i * dksi2 + djk2;
						dtemp != 0 ? dtemp1 = alpha / dtemp : dtemp1 = 0.0; // protection against division by zero
						pout[m][0] *= dtemp1;
						pout[m][1] *= dtemp1;
						m++;
					}
				}
			}

			// inverse FFT of the product
			fftf.InverseFFT();
			fftf.GetRealXArray3D(K3out);

			// restore L1 norm
			K3out *= dnorm1 / K3out.Norm(eNormL1); 
		}

		// output the 3D array
		if (imodeReproj != 2)
		{
			printf("\n\nSaving the reconstructed 3D object into output files ...");
			#pragma omp parallel for shared(K3out)
			for (int n = 0; n < noutdefocus; n++)
			{
				try
				{
					XArray2D<double> ipOut(ny, nx);
					ipOut.SetHeadPtr(pHead ? pHead->Clone() : 0);

					for (index_t j = 0; j < ny; j++)
						for (index_t i = 0; i < nx; i++)
							ipOut[j][i] = K3out[n][j][i];

					if (bVerboseOutput) printf("\nOutput defocus slice no. = %d; output file = %s", n, voutfilenamesTot[n].c_str());
					if (bTIFFoutput) TIFFWriteFile(ipOut, voutfilenamesTot[n].c_str(), eTIFF32);
					else XArData::WriteFileGRD(ipOut, voutfilenamesTot[n].c_str(), eGRDBIN);
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
		}

		// multislice reprojection
		if (imodeReproj == 1 || imodeReproj == 2)
		{
			printf("\n\nPerforming forward re-propagation through the 3D distribution of the electrostatic potential ...");
			double yc = (ylo + yhi) / 2.0, xc = (xlo + xhi) / 2.0;
			K3out *= PI / (wl * EE) * zst; // this converts the potential V into -2 * PI / lambda  * delta * zst, with zst required in the subsequent integration inside Multislice_eV
			xar::XArray3DSpln<double> spln3d(K3out);

			vector<vector<double> > vverr(nangles); // reconstruction errors in individual defocus planes at different illumination angles
			vector<double> verraver(nangles, 0.0); // reconstruction errors at different illumination angles averaged over the corresponding defocus distances
			vector<vector<double> > vvnorm0(nangles); // C0 norms of reprojected intensities-1 in individual defocus planes at different illumination angles
			vector<double> vnormaver0(nangles, 0.0); // C0 norms of reprojected intensities-1 different illumination angles averaged over the corresponding defocus distances
			vector<vector<double> > vvnorm1(nangles); // C0 norms of original intensities-1 in individual defocus planes at different illumination angles
			vector<double> vnormaver1(nangles, 0.0); // C0 norms of original intensities-1 different illumination angles averaged over the corresponding defocus distances

			#pragma omp parallel for shared (K3out, spln3d, vndefocus, vvdefocus, v2angles, vvoutfilenamesCAmp)
			for (int na = 0; na < nangles; na++)
			{
				try
				{
					vector<Pair> vdefocus;
					vector<string> voutfilenamesCAmp, vinfilenames1;
					double angleY(0);
					double angleX(0);
					#pragma omp critical
					{
						vdefocus = vvdefocus[na];
						vverr[na].resize(vndefocus[na]);
						vvnorm0[na].resize(vndefocus[na]);
						vvnorm1[na].resize(vndefocus[na]);
						voutfilenamesCAmp = vvoutfilenamesCAmp[na];
						vinfilenames1 = vvinfilenames1[na];
						angleY = v2angles[na].a * PI180;
						angleX = v2angles[na].b * PI180;
						if (!bVerboseOutput) printf("\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleY / PI180, angleX / PI180);
					}

					// multislice propagate through the reconstructed 3D distribution of the electrostatic potential
					XArray2D<dcomplex> campOut(ny, nx, dcomplex(1.0, 0.0)); // "incident" complex amplitude
					campOut.SetHeadPtr(pHead->Clone());
					spln3d.Multislice_eV(campOut, angleY, angleX, sliceTh, k2maxo);

					// propagate to defocuse planes, rotate around the illumination axis and save the defocused intensities
					XArray2D<double > vint0n, K0n, vint1n, K1n, ampRe0, ampIm0, ampRe, ampIm;
					for (int n = 0; n < vndefocus[na]; n++)
					{
						XArray2D<dcomplex> campOut1(campOut);
						xar::XArray2DFFT<double> xafft(campOut1);
						xafft.Fresnel(vdefocus[n].b, false, k2maxo, Cs3, Cs5); // propagate to the current defocus distance
						if (vdefocus[n].a != 0) // rotate around Z" if needed
						{
							Re(campOut1, ampRe0); Im(campOut1, ampIm0);
							XArray2DSpln<double> xaSplnRe(ampRe0), xaSplnIm(ampIm0);
							xaSplnRe.Rotate(ampRe, vdefocus[n].a, yc, xc, 1.0f); // expecting uniform background with unit amplitude
							xaSplnIm.Rotate(ampIm, vdefocus[n].a, yc, xc, 0.0f); // expecting uniform background with unit amplitude
							MakeComplex(ampRe, ampIm, campOut1, false);
						}
						Abs(campOut1, vint1n);
						Abs2(campOut1, K1n); K1n -= 1.0; 
						// read in the original intensity
						if (GetFileExtension(vinfilenames1[n]) == string(".GRC"))
						{
							XArray2D<dcomplex> campOut2;
							XArData::ReadFileGRC(campOut2, vinfilenames1[n].c_str());
							Abs(campOut2, vint0n);
							Abs2(campOut2, K0n); K0n -= 1.0;
						}
						else
						{
							if (bTIFFinput)
							{
								TIFFReadFile(vint0n, vinfilenames1[n].c_str());
								vint0n.SetHeadPtr(pHead->Clone());
							}
							else
								XArData::ReadFileGRD(vint0n, vinfilenames1[n].c_str(), wl);
							K0n = vint0n; K0n -= 1.0;
							vint0n ^= 0.5;
						}
						// compare the reprojected and the original defocused intensities
						vvnorm0[na][n] = K0n.Norm(eNormC0);
						vnormaver0[na] += vvnorm0[na][n];
						vvnorm1[na][n] = K1n.Norm(eNormC0);
						vnormaver1[na] += vvnorm1[na][n];
						vint1n -= vint0n;
						vverr[na][n] = pow(vint1n.Norm(eNormL2), 2.0) / pow(vint0n.Norm(eNormL2), 2.0); // this error norm is modelled after the one in IWFR
						verraver[na] += vverr[na][n];

						//ReplaceModulus(campOut1, vint0n);
						//Conjg(campOut1); // @@@@@@@@@@@@@@@@ ?????

						#pragma omp critical
						{
							if (bVerboseOutput)
							{
								printf("\n\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleY / PI180, angleX / PI180);
								printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[n].a, vdefocus[n].b);
								printf("\nWriting output file %s ...", voutfilenamesCAmp[n].c_str());
							}
						}
						XArData::WriteFileGRC(campOut1, voutfilenamesCAmp[n].c_str(), eGRCBIN); //	write output GRC files
					}
					verraver[na] /= double(vndefocus[na]);
					vnormaver0[na] /= double(vndefocus[na]);
					vnormaver1[na] /= double(vndefocus[na]);
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

			// calculate the current reconstruction error and
			// if the difference with the previous error is smaller than the defined minimum
			// or, if the error started to increase, interrupt the iterations
			double ssej = 0.0;
			for (index_t na = 0; na < nangles; na++) ssej += verraver[na];
			ssej /= double(nangles);

			printf("\n\nAverage relative error in reprojected defocused intensities = %g", ssej);
			for (index_t na = 0; na < nangles; na++) printf("\nIllum_angle_no. = %zd, C0(orig_contrast) = %g, C0(reproj_contrast) = %g, Rel.error = %g ", na, vnormaver0[na], vnormaver1[na], verraver[na]);
		}
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