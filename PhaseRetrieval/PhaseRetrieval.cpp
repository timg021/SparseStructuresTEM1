// PhaseRetrieval.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// The IWFR algorithm is implemented here according to the description in the paper
// L. Allen, "Electron Microscope Cs Correction Using Iterative Wave - Function Reconstruction", Microscopy and Analysis 20(4):15-17 (UK), 2006

#include <complex.h>
#include <chrono>
#include <omp.h>

#include "IXAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_filt2.h"

using namespace xar;

#define BILINEAR_INTERPOLATION // if this is not defined (commented out), nearest neighbour interpolation code is used in 3D reconstruction


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting IWFR PhaseRetrieval program ...");
		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024];
		FILE* ff0 = fopen("PhaseRetrieval.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PhaseRetrieval.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 1. CT rotation span in degrees and number of rotation angles
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading rotation span and number of rotations from input parameter file.");
		double angle_span = atof(cparam) / 180.0 * PI; // total rotation span in radians 
		index_t nangles = (index_t)atoi(cparam1); // number of rotation steps 
		double angle_step = angle_span / double(nangles); // rotation step in radians
		// !!!Later on, input parameter lines 1 and 2 (rotation angles and defocus distances) should be replaced by a single name of an input text (CSV-type) file
		// which will contain one row per rotational position, each row containing a rotation angle followed by an arbitrary number of defocus distances at that angle

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2. Input defocus distances in Angstroms
		std::vector<size_t> vwhite(0); // vector of white spaces separating different defocus distances (there should be exactly one white space before each defocus distance and no spaces at the end)
		for (size_t i = 1; i < strlen(cline); i++) if (isspace(cline[i])) vwhite.push_back(i); // count the number of different defocus distances in the input parameter file
		size_t ndefocus = vwhite.size(); // number of detected defocus distances
		if (!ndefocus) throw std::exception("Error reading defocus distances from input parameter file.");
		std::vector<char[1024]> vstr_defocus(ndefocus); // vector of strings into which defocus distances will be read
		std::vector<double> vdefocus(ndefocus); // vector of defocus distances (double precision numbers)
		vwhite.push_back(strlen(cline)); // add one more entry corresponding to the end of the parameter string
		printf("\nInput defocus plane positions (%zd in total): ", ndefocus);
		for (size_t j = 0; j < ndefocus; j++)
		{
			for (size_t i = vwhite[j]; i < vwhite[j + 1]; i++)
				vstr_defocus[j][i - vwhite[j]] = cline[i];
			vstr_defocus[j][vwhite[j + 1] - vwhite[j]] = '\0'; // string terminator
			vdefocus[j] = atof(vstr_defocus[j]);
			printf("%g ", vdefocus[j]);
		}
		double zmiddle(0.0); // "middle z plane" position
		for (size_t j = 0; j < ndefocus; j++) zmiddle += vdefocus[j];
		zmiddle /= double(ndefocus);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3. Input_filename_base_of_defocus_series_of_the_sample_in_GRD_format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		printf("\nInput defocus series file name base = %s", filenamebaseIn.c_str());
		vector<string> vinfilenamesTot;
		FileNames(nangles, ndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in Angstroms
		printf("\nWavelength = %g (A)", wl);
		if (wl < 0 || wl > 1)
			throw std::exception("Error: wavelength value appears to be wrong.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam); 
		printf("\nObjective aperture = %g (mrad)", aobj);
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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Maximal number of IWFR iterations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading maximal number of iterations from input parameter file.");
		int kmax = atoi(cparam);
		printf("\nMaximal number of iterations = %d", kmax);
		if (kmax < 1)
			throw std::exception("Error: the maximal number of iterations should be >= 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Minimal phase reconstruction error
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading minimal phase reconstruction error from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nMinimal phase reconstruciton error = %g", epsilon);
		if (epsilon < 0)
			throw std::exception("Error: minimal phase reconstruction error must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9. Output defocus distances min max and step in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading output defocus distances from input parameter file.");
		double zlo = atof(cparam); // minimum output defocus in Angstroms 
		double zhi = atof(cparam1); // maximum output defocus in Angstroms 
		double zst = abs(atof(cparam2)); // output defocus step in Angstroms 
		if (zlo > zhi) std::swap(zlo, zhi);
		printf("\nOutput defocus distances: min = %g, max = %g, step = %g (Angstroms)", zlo, zhi, zst);
		int noutdefocus = int((zhi - zlo) / zst + 0.5); // number of defocus planes to propagate to
		if (noutdefocus <= 0)
			throw std::exception("Error: number of output defocus planes must be positive.");
		vector<double> voutdefocus(noutdefocus); // vector of output defocus distances
		printf("\nOutput defocus plane positions (%d in total): ", noutdefocus);
		for (int j = 0; j < noutdefocus; j++)
		{
			voutdefocus[j] = zlo + zst * j;
			printf("%g ", voutdefocus[j]);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10. Output intensity(0), phase(1), complex_amplitude(2) or 3D contrast(3)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output files format from input parameter file.");
		int noutformat = atoi(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11. Output file name base in GRD or GRC format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		printf("\nOutput file name base = %s", filenamebaseOut.c_str());
		switch (noutformat)
		{
		case 0:
			printf("\nThe program will output defocused intensity distributions in GRD format.");
			if (GetFileExtension(filenamebaseOut) != string(".GRD"))
				throw std::exception("Error: output filename extension for intensity distributions must be grd or GRD.");
			break;
		case 1:
			printf("\nThe program will output defocused phase distributions in GRD format.");
			if (GetFileExtension(filenamebaseOut) != string(".GRD"))
				throw std::exception("Error: output filename extension for phase distributions must be grd or GRD.");
			break;
		case 2:
			printf("\nThe program will output defocused complex amplitudes in GRC format.");
			if (GetFileExtension(filenamebaseOut) != string(".GRC"))
				throw std::exception("Error: output filename extension for complex amplitides must be grc or GRC.");
			break;
		case 3:
			printf("\nThe program will output 3D defocused contrast in a series of GRD files.");
			if (GetFileExtension(filenamebaseOut) != string(".GRD"))
				throw std::exception("Error: output filename extension for complex amplitides must be grd or GRD.");
			break;
		default:
			throw std::exception("Error: unknown value for output file format in input parameter file.");
		}
		vector<string> voutfilenamesTot;
		if (noutformat == 3)
			FileNames(1, noutdefocus, filenamebaseOut, voutfilenamesTot); // create "total 1D array" of output filenames
		else
			FileNames(nangles, noutdefocus, filenamebaseOut, voutfilenamesTot); // create "total 2D array" of output filenames
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file
		double xlo, xhi, xst;
		XArray3D<double> K3Out; // big 3D reconstructed array (needs to fit into RAM alongside with with everything else)

		// start of cycle over rotation angles
		for (index_t na = 0; na < nangles; na++) 
		{
 			double angle = angle_step * double(na);
			double cosangle = cos(angle);
			double sinangle = sin(angle);
			printf("\n\n*** Rotation angle[%zd] = %g (degrees)", na, angle / PI * 180.0);
			// We use the fact that currently the number of defocus planes, ndefocus, is assumed to be the same for all angles, and so we don't change it here.
			// Similarly, the vectors of input and output defocus distances, vdefocus[] and voutdefocus[], are assumed to be the same for all angles.
			// Filenames for the input and output defocus images are different for each angle, and so they need to be adjusted here.
			vector<string> vinfilenames(ndefocus);
			for (index_t n = 0; n < ndefocus; n++) vinfilenames[n] = vinfilenamesTot[na * ndefocus + n];
			vector<string> voutfilenames(noutdefocus);
			if (noutformat != 3) 
				for (index_t n = 0; n < noutdefocus; n++) voutfilenames[n] = voutfilenamesTot[na * noutdefocus + n];

			// define main work objects
			double dtemp; // auxilliary variable
			double ssej(0.0), ssejm1(0.0); // current and previous average reconstruction errors
			vector<double> verr(ndefocus); // reconstruction errors in individual defocus planes
			vector<XArray2D<double>> vint0(ndefocus), vint(ndefocus); // initial and iterated defocused intensities
			vector<XArray2D<dcomplex>> vcamp(ndefocus); // defocused complex amplitudes
			XArray2D<dcomplex> campOut; // complex amplitude in the plane z=0
			vector<double> vint0_L1(ndefocus); // L1 norms of the initial defocused intensities

			// start point of IWFR iterations
			for (int k = 0; k < kmax; k++)
			{
				// apply initial or newly reconstructed phases (the second case is equal to restoring the original moduli)
				// and propagate each defocused amplitude to the "middle" plane z = zmiddle
				#pragma omp parallel for
				for (int n = 0; n < ndefocus; n++)
				{
					try
					{
						if (k == 0) // read input defocused intensity and create initial defocused complex amplitude
						{
							XArData::ReadFileGRD(vint0[n], vinfilenames[n].c_str(), wl); //	read input GRD files
							vint0_L1[n] = vint0[n].Norm(eNormL1);
							printf("\nL1 norm of input defocused intensity no. %d = %g", n, vint0_L1[n]);
							if (vint0_L1[n] == 0) throw std::exception("Error: input intensity file is empty");
							vint0[n] ^= 0.5; // intensity --> real amplitude
							MakeComplex(vint0[n], 0.0, vcamp[n], true); // apply initial zero phases
						}
						else // apply phases obtained on the previous iteration
							for (index_t j = 0; j < vcamp[n].GetDim1(); j++)
								for (index_t i = 0; i < vcamp[n].GetDim2(); i++)
								{
									dtemp = abs(vcamp[n][j][i]);
									if (dtemp) vcamp[n][j][i] *= vint0[n][j][i] / dtemp;
									else vcamp[n][j][i] = std::polar(vint0[n][j][i], 0.0);
								}
						xar::XArray2DFFT<double> xafft(vcamp[n]);
						xafft.Fresnel(zmiddle - vdefocus[n], false, k2maxo, Cs3, Cs5); // propagate to z = 0
					}
					catch (std::exception& E)
					{
						printf("\n\n!!!Exception: %s\n", E.what());
						throw;
					}
				}

				if (k == 0 && noutformat == 3)
				{
					IXAHWave2D* ph2 = GetIXAHWave2D(vcamp[0]);
					xlo = ph2->GetXlo();
					xhi = ph2->GetXhi();
					xst = (xhi - xlo) / vcamp[0].GetDim2();
					if (xst != zst)
						throw std::exception("Error: the 3D reconstructed object is supposed to have square pixels");
				}

				// average complex amplitudes in the middle plane
				campOut = vcamp[0];
				for (index_t n = 1; n < ndefocus; n++) campOut += vcamp[n];
				campOut /= double(ndefocus);

				// propagate the averaged complex amplitude from the middle plane to the individual defocus planes
				#pragma omp parallel for shared(campOut)
				for (int n = 0; n < ndefocus; n++)
				{
					try
					{
						vcamp[n] = campOut;
						xar::XArray2DFFT<double> xafft(vcamp[n]);
						xafft.Fresnel(vdefocus[n] - zmiddle, false, k2maxo, Cs3, Cs5); // propagate to z = z[n]
						Abs(vcamp[n], vint[n]);
						vint[n] -= vint0[n];
						verr[n] = pow(vint[n].Norm(eNormL2), 2.0) / vint0_L1[n];
					}
					catch (std::exception& E)
					{
						printf("\n\n!!!Exception: %s\n", E.what());
						throw;
					}
				}

				// calculate the current reconstruction error and
				// if the difference with the previous error is smaller than the defined minimum
				// or, if the error started to increase, interrupt the iterations
				ssej = 0.0;
				for (index_t n = 0; n < ndefocus; n++) ssej += verr[n];
				ssej /= double(ndefocus);
				
				if (k == 0) printf("\nIteration number %d; SSE_aver error at 0th iteration = %g\n", k, ssej);
				else printf("\nIteration number %d; SSE_aver error difference with previous iteration = %g\n", k, ssejm1 - ssej);
				
				for (index_t n = 0; n < ndefocus; n++) printf("SSE(%zd) = %g ", n, verr[n]);
				
				if (k > 0 && (ssejm1 - ssej) < epsilon) break;
				else ssejm1 = ssej;
			}
			// end point of IWFR iterations

			// now start calculations of the output defocused images
			printf("\nPropagating to output defocus planes ...");
			
			vcamp.resize(noutdefocus); // repurpose the vector of 2D complex XArrays

			#pragma omp parallel for
			for (int n = 0; n < noutdefocus; n++) 
			{
				try
				{
					// propagate to the output defocused plane
					vcamp[n] = campOut;
					xar::XArray2DFFT<double> xafft(vcamp[n]);
					xafft.Fresnel(voutdefocus[n] - zmiddle, false, k2maxo, Cs3, Cs5); // propagate to z_out[n]
					
					// write the defocused intensity, phase or complex amplitude into a separate 2D file (unless noutformat == 3)
					XArray2D<double> ipOut;
					switch (noutformat)
					{
					case 0: // intensity out
						Abs2(vcamp[n], ipOut);
						printf("\nOutput defocus distance = %g; output file = %s", voutdefocus[n], voutfilenames[n].c_str());
						XArData::WriteFileGRD(ipOut, voutfilenames[n].c_str(), eGRDBIN);
						break;
					case 1: // phase out
						CArg(vcamp[n], ipOut);
						printf("\nOutput defocus distance = %g; output file = %s", voutdefocus[n], voutfilenames[n].c_str());
						XArData::WriteFileGRD(ipOut, voutfilenames[n].c_str(), eGRDBIN);
						break;
					case 2: // complex amplitude out
						printf("\nOutput defocus distance = %g; output file = %s", voutdefocus[n], voutfilenames[n].c_str());
						XArData::WriteFileGRC(vcamp[n], voutfilenames[n].c_str(), eGRCBIN);
						break;
					case 3: // 3D output of the contrast function
						// in this case the output takes place only once, at the last rotation angle
						break;
					default:
						throw std::exception("Error: unknown output format");
					}
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					throw;
				}
			}


			if (noutformat == 3) // add the output defocused data obtained at the current rotational angle to the 3D object that is being reconstructed
			{
				printf("\nUpdating 3D reconstructed object ...");
				index_t ny = vcamp[0].GetDim1();
				index_t nx = vcamp[0].GetDim2();

				double xc = (xhi + xlo) / 2.0; // x-coordinate of the centre of rotation
				double zc = (zhi + zlo) / 2.0; // z-coordinate of the centre of rotation

				// calculate the coordinate rotation angle parameters
				double xxx;
				vector<double> x_sinangle(nx), x_cosangle(nx);
				for (index_t i = 0; i < nx; i++)
				{
					xxx = xlo + xst * i - xc;
					x_sinangle[i] = xxx * sinangle;
					x_cosangle[i] = xxx * cosangle;
				}

				// allocate the large 3D output array
				if (na == 0) K3Out.Resize(noutdefocus, ny, nx, 0.0);

				// rotate the defocus plane around the "vertical" y axis by -angle, instead of rotating the 3D object by the angle
				index_t ii, nn;
				double zzz, z_sinangle, z_cosangle;
#if defined(BILINEAR_INTERPOLATION)
				index_t nx2 = nx - 2, noutdefocus2 = noutdefocus - 2;
				double dK, dx0, dx1, dz0, dz1;
#else
				index_t nx1 = nx - 1, noutdefocus1 = noutdefocus - 1;
#endif
				for (int n = 0; n < noutdefocus; n++)
				{
					z_sinangle = (voutdefocus[n] - zc) * sinangle;
					z_cosangle = (voutdefocus[n] - zc) * cosangle;
					for (index_t i = 0; i < nx; i++)
					{
						xxx = xc + x_cosangle[i] - z_sinangle; // x coordinate with respect to the rotated 3D sample
						zzz = zc + x_sinangle[i] + z_cosangle; // z coordinate with respect to the rotated 3D sample
#if defined(BILINEAR_INTERPOLATION)
						dx0 = abs(xxx - xlo) / xst; ii = (index_t)dx0; dx0 -= ii; dx0 *= 0.5; dx1 = 0.5 - dx0; // "bilinear" interpolation variant
						dz0 = abs(zzz - zlo) / zst; nn = (index_t)dz0; dz0 -= nn; dz0 *= 0.5; dz1 = 0.5 - dz0; // "bilinear" interpolation variant
						if (ii > nx2 || nn > noutdefocus2) continue;
#else
						if (xxx < xlo) xxx = xlo;
						if (zzz < zlo) zzz = zlo;
						ii = (index_t)((xxx - xlo) / xst + 0.5); // nearest neighbour interpolation variant
						nn = (index_t)((zzz - zlo) / zst + 0.5); // nearest neighbour interpolation variant
						if (ii > nx1) ii = nx1;
						if (nn > noutdefocus1) nn = noutdefocus1;
#endif
						for (index_t j = 0; j < ny; j++)
						{
#if defined(BILINEAR_INTERPOLATION)
							dK = 1.0 - std::norm(vcamp[n][j][i]); // contrast value at this point of the defocused plane
							K3Out[nn][j][ii] += (dx0 + dz0) * dK; // "bilinear" interpolation variant
							K3Out[nn][j][ii + 1] += (dx1 + dz0) * dK; // "bilinear" interpolation variant
							K3Out[nn + 1][j][ii] += (dx0 + dz1) * dK; // "bilinear" interpolation variant
							K3Out[nn + 1][j][ii + 1] += (dx1 + dz1) * dK; // "bilinear" interpolation variant
#else
							K3Out[nn][j][ii] += 1.0 - std::norm(vcamp[n][j][i]); // nearest neigbour interpolation variant
							//K3Out[nn][j][ii] += std::norm(vcamp[n][j][i]); // nearest neigbour interpolation variant @@@@@ temporary
#endif
						}
					}
				}

				// output the 3D array
				if (na == (nangles - 1))
				{
					K3Out /= double(nangles);
					XArray2D<double> ipOut(ny, nx);
					ipOut.SetHeadPtr(vint0[0].GetHeadPtr() ? vint0[0].GetHeadPtr()->Clone() : 0);
					XArray2DFilt<double> xafilt(ipOut);
					for (int n = 0; n < noutdefocus; n++)
					{
						for (index_t j = 0; j < ny; j++)
							for (index_t i = 0; i < nx; i++)
								ipOut[j][i] = K3Out[n][j][i];
						xafilt.AverageFilt(5, 0);
						printf("\nOutput defocus slice no. = %d; output file = %s", n, voutfilenamesTot[n].c_str());
						XArData::WriteFileGRD(ipOut, voutfilenamesTot[n].c_str(), eGRDBIN);
					}
				}
			}

		} // end of cycle over rotation angles
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