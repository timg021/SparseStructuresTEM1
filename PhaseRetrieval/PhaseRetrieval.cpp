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
#include "fftwd3drc.h"

using namespace xar;

#define TRILINEAR_INTERPOLATION // if this is not defined (commented out), nearest neighbour interpolation code is used in 3D reconstruction

//void TriangularFilter(vector<double>& xarr, int nfilt2);

int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting IWFR PhaseRetrieval program ...");
		vector<Pair> v2angles;
		vector<vector <double> > vvdefocus;

		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024];
		FILE* ff0 = fopen("PhaseRetrieval.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PhaseRetrieval.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 1. Input_file_with_rotation_angles_and_defocus_distances
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading file name with rotation angles and defocus distances from input parameter file.");
		ReadDefocusParamsFile(cparam, v2angles, vvdefocus);
		index_t nangles = v2angles.size(); // number of rotation steps 
		vector<index_t> vndefocus(nangles); // vector of numbers of defocus planes at different rotation angles
		for (index_t i = 0; i < nangles; i++) vndefocus[i] = vvdefocus[i].size();

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2. Input_filename_base_of_defocus_series_of_the_sample_in_GRD_format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		printf("\nInput defocus series file name base = %s", filenamebaseIn.c_str());
		vector<string> vinfilenamesTot;
		// FileNames2(vndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames - this has to be delayed until the 3D Laplacian filter mode is known
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3. Wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in Angstroms
		printf("\nWavelength = %g (A)", wl);
		if (wl < 0 || wl > 1)
			throw std::exception("Error: wavelength value appears to be wrong.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam); 
		printf("\nObjective aperture = %g (mrad)", aobj);
		if (aobj < 0 || aobj > 1000)
			throw std::exception("Error: objective aperture value appears to be wrong.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. Spherical aberrations Cs3 and Cs5 in mm
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading spherical aberrations from input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> Angstroms
		Cs5 *= 1.e+7; // mm --> Angstroms

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6. Maximal number of IWFR iterations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading maximal number of iterations from input parameter file.");
		int kmax = atoi(cparam);
		printf("\nMaximal number of iterations = %d", kmax);
		if (kmax < 1)
			throw std::exception("Error: the maximal number of iterations should be >= 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Minimal phase reconstruction error
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading minimal phase reconstruction error from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nMinimal phase reconstruciton error = %g", epsilon);
		if (epsilon < 0)
			throw std::exception("Error: minimal phase reconstruction error must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Output defocus distances min max and step in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading output defocus distances from input parameter file.");
		double zlo = atof(cparam); // minimum output defocus in Angstroms - !!! will be corrected with dzextra below
		double zhi = atof(cparam1); // maximum output defocus in Angstroms - !!! will be corrected with dzextra below 
		double zst = abs(atof(cparam2)); // output defocus step in Angstroms 
		if (zlo > zhi) std::swap(zlo, zhi);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //9. Extra defocus for 3D reconstruction in_Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading extra defocus for 3D reconstruction parameter from input parameter file.");
		double dzextra = atof(cparam);
		printf("\nExtra defocus for 3D reconstruction = %g (Angstroms)", dzextra);
		zlo += dzextra; zhi += dzextra;

		printf("\nOutput defocus distances: min = %g, max = %g, step = %g (Angstroms)", zlo, zhi, zst);
		if (zst == 0)
			throw std::exception("Error: output defocus step cannot be zero.");
		int noutdefocus = int((zhi - zlo) / zst + 0.5); // number of defocus planes to propagate to
		if (noutdefocus <= 0)
			throw std::exception("Error: number of output defocus planes must be positive.");
		vector<double> voutdefocus(noutdefocus); // vector of output defocus distances
		printf("\nOutput defocus plane positions (%d in total): ", noutdefocus);
		for (index_t n = 0; n < noutdefocus; n++)
		{
			voutdefocus[n] = zlo + zst * n;
			printf("%g ", voutdefocus[n]);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10. Regularization parameter for inverse 3D Laplacian
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
		if (imode3Dfilter != 2)
			FileNames2(vndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11. Regularization parameter for inverse 3D Laplacian
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading regularization parameter for inverse 3D Laplacian from input parameter file.");
		double alpha = atof(cparam);
		printf("\nRegularization parameter for inverse 3D Laplacian = %g", alpha);
		if (alpha < 0 && (imode3Dfilter == 1 || imode3Dfilter == 2))
			throw std::exception("Error: regularization parameter for 3D Laplacian filter must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12. Output intensity(0), phase(1), complex_amplitude(2) or 3D contrast(3)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output files format from input parameter file.");
		int noutformat = atoi(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13. Output file name base in GRD or GRC format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		printf("\nOutput file name base = %s", filenamebaseOut.c_str());
		switch (noutformat)
		{
		case 0:
			if (imode3Dfilter == 1 || imode3Dfilter == 2)
				throw std::exception("Error: when using 3D inverse Laplacian filter, output format must be 3D contrast(3).");
			printf("\nThe program will output defocused intensity distributions in GRD format.");
			if (GetFileExtension(filenamebaseOut) != string(".GRD"))
				throw std::exception("Error: output filename extension for intensity distributions must be grd or GRD.");
			break;
		case 1:
			if (imode3Dfilter == 1 || imode3Dfilter == 2)
				throw std::exception("Error: when using 3D inverse Laplacian filter, output format must be 3D contrast(3).");
			printf("\nThe program will output defocused phase distributions in GRD format.");
			if (GetFileExtension(filenamebaseOut) != string(".GRD"))
				throw std::exception("Error: output filename extension for phase distributions must be grd or GRD.");
			break;
		case 2:
			if (imode3Dfilter == 1 || imode3Dfilter == 2)
				throw std::exception("Error: when using 3D inverse Laplacian filter, output format must be 3D contrast(3).");
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
			if (imode3Dfilter == 2) // only read "output" files, filter them and write in different files
			{
				FileNames(1, noutdefocus, filenamebaseOut, vinfilenamesTot); // create 1D array of input filenames to read the input 2D slices of a previously reconstructed 3D object
				string filenamebaseOutNew = filenamebaseOut.insert(filenamebaseOut.find_last_of("."), "L");
				FileNames(1, noutdefocus, filenamebaseOutNew, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the filtered 3D object
			}
			else
				FileNames(1, noutdefocus, filenamebaseOut, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
		else
			// note that the number of output defocus distances is assumed to be the same at all rotation angles
			FileNames(nangles, noutdefocus, filenamebaseOut, voutfilenamesTot); // create total 2D array of output filenames to save output 2D defocused images at different rotation angles and output defocus distances
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		bool bAbort(false);
		double xlo, xhi, xst;
		double ylo, yhi, yst;
		XArray2D<double> aaa;
		IXAHead* pHead(0);
		XArray3D<double> K3Out; // big 3D reconstructed array (needs to fit into RAM alongside with with everything else)

		if (imode3Dfilter != 2) // do phase retrieval and backpropagation prior to 3D filtering and output
		{
			// start of cycle over rotation angles
			index_t ndefcurrent(0);
			for (index_t na = 0; na < nangles; na++)
			{
				double angleY = v2angles[na].y * PI180;
				double angleX = v2angles[na].x * PI180;
				double cosangleY = cos(angleY);
				double sinangleY = sin(angleY);
				double cosangleX = cos(angleX);
				double sinangleX = sin(angleX);
				printf("\n\n*** Rotation angle[%zd] = (%g, %g) (degrees)", na, angleY / PI180, angleX / PI180);

				index_t ndefocus = vndefocus[na]; // number of defocus planes at the current rotation angle
				vector<double> vdefocus = vvdefocus[na]; // vector of input defocus positions at the current defocus angle
				double zmiddle(0.0); // "middle z plane" position
				for (size_t j = 0; j < ndefocus; j++) zmiddle += vdefocus[j];
				zmiddle /= double(ndefocus);
				// The vector of output defocus distances, voutdefocus[], is assumed to be the same for all angles.
				// Filenames for the input and output defocus images are different for each angle, and so they need to be adjusted here.
				vector<string> vinfilenames(ndefocus);
				for (index_t n = 0; n < ndefocus; n++) vinfilenames[n] = vinfilenamesTot[ndefcurrent++];
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
								if (noutformat == 3)
								{
									IXAHWave2D* ph2 = GetIXAHWave2D(vint0[n]);
									xlo = ph2->GetXlo();
									xhi = ph2->GetXhi();
									xst = (xhi - xlo) / vint0[n].GetDim2();
									ylo = ph2->GetYlo();
									yhi = ph2->GetYhi();
									yst = (yhi - ylo) / vint0[n].GetDim1();
									if (xst != yst || xst != zst)	throw std::exception("Error: the 3D reconstructed object is supposed to have qubic voxels");
								}
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
							bAbort = true;
						}
					}
					if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

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
							bAbort = true;
						}
					}
					if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

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
						bAbort = true;
					}
				}
				if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

				if (noutformat == 3) // add the output defocused data obtained at the current rotational angle to the 3D object that is being reconstructed
				{
					printf("\nUpdating 3D reconstructed object ...");
					index_t ny = vcamp[0].GetDim1();
					index_t nx = vcamp[0].GetDim2();

					double xc = (xhi + xlo) / 2.0; // x-coordinate of the centre of rotation
					double yc = (yhi + ylo) / 2.0; // y-coordinate of the centre of rotation
					double zc = (zhi + zlo) / 2.0; // z-coordinate of the centre of rotation

					// calculate the coordinate rotation angle parameters
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
						K3Out.Resize(noutdefocus, ny, nx, 0.0);
						pHead = vint0[0].GetHeadPtr();
					}

					// rotate the defocus plane around the "vertical" y axis by -angle, instead of rotating the 3D object by the angle
					index_t ii, jj, nn;
					double zz;
#if defined(TRILINEAR_INTERPOLATION)
					index_t nx2 = nx - 2, ny2 = ny - 2, noutdefocus2 = noutdefocus - 2;
					double dK, dx0, dx1, dy0, dy1, dz0, dz1, dz0K, dz1K;
#else
					index_t nx1 = nx - 1, ny1 = ny - 1, noutdefocus1 = noutdefocus - 1;
#endif
					for (int n = 0; n < noutdefocus; n++)
					{
						for (index_t j = 0; j < ny; j++)
						{
							// inverse rotation around X' axis
							yyy = yc + y_cosangleX[j] - z_sinangleX[n]; // y coordinate with respect to the rotated 3D sample
							zz = y_sinangleX[j] + z_cosangleX[n];
							for (index_t i = 0; i < nx; i++)
							{
								// inverse rotation around Y axis
								xxx = xc + x_cosangleY[i] - zz * sinangleY; // x coordinate with respect to the rotated 3D sample
								zzz = zc + x_sinangleY[i] + zz * cosangleY; // z coordinate with respect to the rotated 3D sample
#if defined(TRILINEAR_INTERPOLATION)
								dx1 = abs(xxx - xlo) / xst; ii = (index_t)dx1;
								dy1 = abs(yyy - ylo) / yst; jj = (index_t)dy1;
								dz1 = abs(zzz - zlo) / zst; nn = (index_t)dz1;
								if (ii > nx2 || jj > ny2 || nn > noutdefocus2) continue;
								dx1 -= ii; dx0 = 1.0 - dx1;
								dy1 -= jj; dy0 = 1.0 - dy1;
								dz1 -= nn; dz0 = 1.0 - dz1;
#else
								if (xxx < xlo) xxx = xlo;
								if (yyy < ylo) yyy = ylo;
								if (zzz < zlo) zzz = zlo;
								ii = (index_t)((xxx - xlo) / xst + 0.5); // nearest neighbour interpolation variant
								jj = (index_t)((yyy - ylo) / yst + 0.5); // nearest neighbour interpolation variant
								nn = (index_t)((zzz - zlo) / zst + 0.5); // nearest neighbour interpolation variant
								if (ii > nx1) ii = nx1;
								if (jj > ny1) jj = ny1;
								if (nn > noutdefocus1) nn = noutdefocus1;
#endif
#if defined(TRILINEAR_INTERPOLATION)
								dK = 1.0 - std::norm(vcamp[n][j][i]);
								dz0K = dz0 * dK;
								dz1K = dz1 * dK;
								K3Out[nn][jj][ii] += dx0 * dy0 * dz0K;
								K3Out[nn][jj][ii + 1] += dx1 * dy0 * dz0K;
								K3Out[nn + 1][jj][ii] += dx0 * dy0 * dz1K;
								K3Out[nn + 1][jj][ii + 1] += dx1 * dy0 * dz1K;
								K3Out[nn][jj + 1][ii] += dx0 * dy1 * dz0K;
								K3Out[nn][jj + 1][ii + 1] += dx1 * dy1 * dz0K;
								K3Out[nn + 1][jj + 1][ii] += dx0 * dy1 * dz1K;
								K3Out[nn + 1][jj + 1][ii + 1] += dx1 * dy1 * dz1K;
#else
								K3Out[nn][jj][ii] += 1.0 - std::norm(vcamp[n][j][i]); // nearest neigbour interpolation variant
								//K3Out[nn][jj][ii] += std::norm(vcamp[n][j][i]); // nearest neigbour interpolation variant @@@@@ temporary
#endif
							}
						}
					}
				}

			} // end of cycle over rotation angles
		} // end of case if imode3Dfilter != 2
		else
		{
			// read input 2D slices
			vector<XArray2D<double>> vint(noutdefocus); // input K3 slices
			for (int n = 0; n < noutdefocus; n++)
				XArData::ReadFileGRD(vint[n], vinfilenamesTot[n].c_str(), wl); //	read input GRD files
			IXAHWave2D* ph2 = GetIXAHWave2D(vint[0]);
			xlo = ph2->GetXlo();
			xhi = ph2->GetXhi();
			xst = (xhi - xlo) / vint[0].GetDim2();
			ylo = ph2->GetYlo();
			yhi = ph2->GetYhi();
			yst = (yhi - ylo) / vint[n].GetDim1();
			if (xst != yst || xst != zst)	throw std::exception("Error: the input 3D object is supposed to have qubic voxels");
			pHead = vint[0].GetHeadPtr();

			// fill the input 3D object
			index_t ny = vint[0].GetDim1();
			index_t nx = vint[0].GetDim2();
			K3Out.Resize(noutdefocus, ny, nx, 0.0);
			for (index_t n = 0; n < noutdefocus; n++)
				for (index_t j = 0; j < ny; j++)
					for (index_t i = 0; i < nx; i++)
						K3Out[n][j][i] = vint[n][j][i];
		}

		// apply inverse 3D Laplace filter to the reconstructed 3D array
		if (noutformat == 3)
		{
			printf("\n\n*** Saving the reconstructed 3D object into output files ...");
			K3Out /= double(nangles);

			// regularized inverse 3D Laplacian filter
			if (imode3Dfilter == 1 || imode3Dfilter == 2)
			{
				printf("\nInverse Laplace filtering 3D reconstructed object ...");
				//allocate space for FFT transform and create FFTW plans
				Fftwd3drc fftf((int)noutdefocus, (int)ny, (int)nx);

				// FFT of test array
				fftf.SetRealXArray3D(K3Out);
				fftf.ForwardFFT();

				/// multiply FFT of K3Out arrays by the FFT version of regularized inverse 3D Laplacian
				double fact = 2.0 * PI * wl * (dzextra * 2.0);
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
							dtemp = sin(i * i * dksi2 + djk2);
							dtemp != 0 ? dtemp1 = 1.0 / dtemp : dtemp1 = 0.0; // protection against division by zero
							pout[m][0] *= dtemp1;
							pout[m][1] *= dtemp1;
							m++;
						}
					}
				}

				// inverse FFT of the product
				fftf.InverseFFT();
				fftf.GetRealXArray3D(K3Out);
			}

			// output the 3D array
			#pragma omp parallel for shared(K3Out)
			for (int n = 0; n < noutdefocus; n++)
			{
				try
				{
					XArray2D<double> ipOut(K3Out.GetDim2(), K3Out.GetDim3());
					ipOut.SetHeadPtr(pHead ? pHead->Clone() : 0);

					for (index_t j = 0; j < K3Out.GetDim2(); j++)
						for (index_t i = 0; i < K3Out.GetDim3(); i++)
							ipOut[j][i] = K3Out[n][j][i];

					printf("\nOutput defocus slice no. = %d; output file = %s", n, voutfilenamesTot[n].c_str());
					XArData::WriteFileGRD(ipOut, voutfilenamesTot[n].c_str(), eGRDBIN);
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
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


/*
void TriangularFilter(vector<double>& xarr, int nfilt2)
// Averages 1D array over each nfilt2 * 2 + 1 successive elements using triangular-shaped weight distribution
{
	if (xarr.size() < nfilt2 * 2 + 1) throw std::exception("Unsuitable arguments in TriangleFilter().");
	double fact = 1.0 / double((nfilt2 + 1) * (nfilt2 + 1));
	vector<double> vtemp(xarr.size(), 0.0);
	vector<double> vweight(nfilt2 * 2 + 1);
	for (int j = -nfilt2; j <= nfilt2; j++)	vweight[j + nfilt2] = double(1 + (nfilt2 - abs(j))) * fact;
	
	for (int i = nfilt2; i < xarr.size() - nfilt2; i++)
		for (int j = -nfilt2; j <= nfilt2; j++) 
			vtemp[i] += xarr[i + j] * vweight[j + nfilt2];
	for (int i = nfilt2; i < xarr.size() - nfilt2; i++) xarr[i] = vtemp[i];
}
*/
