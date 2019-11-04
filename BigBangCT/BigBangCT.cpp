// BigBangCT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define CORRELATION_BASED_METHOD 1

#ifdef CORRELATION_BASED_METHOD

#include <complex.h>
#include <chrono>
#include <omp.h>
#include <fftw3.h>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "XA_move3.h"
#include "XA_file.h"
#include "fftwf3drc.h"

//!!! NOTE that fftw3 and XArray3D have the same notation for the order of the array dimensions. Both use C-style array structure, i.e. the last index changes the fastest, 
//i.e. in XArray3D(dim1, dim2, dim3) the fastest changing dimension is dim3, and in fftw_r2c_3d(n0, n1, n2) the fastest dimension is n2.
//Note also that these dimensions are associated with the following physical coordinates by default: dim1(n0) <-> nz, dim2(n1) <-> ny, dim3(n2) <-> nx.


using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting BigBangCT program ...");
		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];
		FILE* ff0 = fopen("BigBangCT.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file BigBangCT.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2nd line: Defocus_distance_MIN,MAX,STEP_in_Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading defocus parameters from input parameter file.");
		double zmin = atof(cparam); // minimum defocus in Angstroms 
		double zmax = atof(cparam1); // maximum defocus in Angstroms 
		double zstep = atof(cparam2); // defocus step in Angstroms
		printf("\nDefocus planes parameters: min = %g, max = %g, step = %g", zmin, zmax, zstep);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // background intensity value
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading background intensity parameter from input parameter file.");
		float finten0 = (float)atof(cparam); // minimum defocus in Angstroms 
		printf("\nBackground intensity value = %g", finten0);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for defocus series of the whole sample
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base for the whole sample from input parameter file.");
		string filenamebaseIn1 = cparam;
		printf("\nDefocus series file name base for the whole sample = %s", filenamebaseIn1.c_str());
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of different atom types
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of atom types from input parameter file.");
		index_t natomtypes = (index_t)atoi(cparam);
		printf("\nNumber of different atom types = %zd", natomtypes);
		vector<index_t> natom(natomtypes);
		vector<string> filenamebaseIn2(natomtypes); // file name bases for defocus series of different single atoms
		vector< vector<double> > rpos2(natomtypes); // vector of XYZ positions of template atoms
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			rpos2[nat].resize(3);
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of atoms of this type and file base name for defocus series of a single atom of this type
			if (sscanf(cline, "%s %s %s %s %s %s", ctitle, cparam, cparam1, cparam2, cparam3, cparam4) != 6) throw std::exception("Error reading atom type %d parameters from input parameter file.", int(nat));
			natom[nat] = index_t(atoi(cparam)); // number of atoms of this type
			rpos2[nat][0] = atof(cparam1); // X-coordinate of template atom no. 'nat'
			rpos2[nat][1] = atof(cparam2); // Y-coordinate of template atom no. 'nat'
			rpos2[nat][2] = atof(cparam3); // Z-coordinate of template atom no. 'nat'
			filenamebaseIn2[nat] = cparam4; // file name base for defocus series of single atom of this type
			printf("\nDefocus series file name base for the single atom no.%zd = %s", nat, filenamebaseIn2[nat].c_str());
		}
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength parameter from input parameter file.");
		double wl = atof(cparam); // wavelength in input file units (usually, Angstroms). Unfortunately, it is not saved in the GRD files
		printf("\nWavelength = %g (A)", wl);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // average atom size for masking out in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading atom size parameter from input parameter file.");
		double atomsize = atof(cparam); // atom diameter in physical units
		printf("\nAtom diameter = %g (A)", atomsize);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // average atom trace Z-length for masking out in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading atom trace Z-length parameters from input parameter file.");
		double atomlength = atof(cparam); // atom "trace" length in the defocus direction in Angstroms to mask "in" the template atom
		double atomsizeZ0 = atof(cparam1); // atom "trace" length in the defocus direction in Angstroms to mask "out" when searching for atoms of the same type
		double atomsizeZ1 = atof(cparam2); // atom "trace" length in the defocus direction in Angstroms to mask "out" when searching for atoms of the next type
		printf("\nAtom trace lengths in Angstroms: 'in' = %g, 'out_same' = %g, 'out_previous' = %g", atomlength, atomsizeZ0, atomsizeZ1);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // high-frequency bandpath radius
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading low frequency limit from input parameter file.");
		int iHPathRad = index_t(atoi(cparam)); // high-frequency bandpath radius: all (abs.)frequencies lower than this one will be zet to zero
		printf("\nHigh frequency lower threshold = %d", iHPathRad);
		index_t iHPathRad2 = index_t(iHPathRad * iHPathRad);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // output file in Vesta XYZ format for detected atom locations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name for detected atom locations from input parameter file.");
		string filenameOutXYZ = cparam;
		printf("\nOutput XYZ file name = %s", filenameOutXYZ.c_str());
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // optional auxillary data output mode 
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading parameter for selecting output mode from input parameter file.");
		int iCorrArrayOut = atoi(cparam); // if this parameter is 1, the 1st masked array is output, 2 -> 2nd masked array output, 3 -> 3D correlation output is created.
		printf("\nOptional file output mode = %d", iCorrArrayOut);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // optional output file name base
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base for saving auxilliary data.");
		string filenamebaseOut = cparam;
		printf("\nAuxilliary output file name baze = %s", filenamebaseOut.c_str());
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("The number of parallel threads in input parameter file should be >= 1.");

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		// calculate some useful parameters
		index_t ndefocus = index_t((zmax - zmin) / zstep + 0.5); // number of defocus planes, it determines the number of input files to read
		if (ndefocus < 1) throw std::exception("Error: the number of defocus planes is less than 1.");
		index_t nz = ndefocus;
		printf("\nNumber of defocus planes = %zd.", nz);
		index_t ny, nx, nx2; // x and y sizes of the arrays
		index_t nangles = 1; // !!! nangles values other than 1 are currently not supported in the code below
		double xmin = 0.0, ymin = 0.0;  // default values - may be overwritten below by data read from input files
		double xstep = 1.0, ystep = 1.0; // default values - may be overwritten below by data read from input files
		index_t karad = index_t(atomsize / zstep / 2.0 + 0.5); // atom radius in the number of physical z-step units
		index_t jarad = index_t(atomsize / ystep / 2.0 + 0.5); // atom radius in the number of physical y-step units - may be overwritten below by data read from input files
		index_t iarad = index_t(atomsize / xstep / 2.0 + 0.5); // atom radius in the number of physical x-step units - may be overwritten below by data read from input files
		index_t karadt = int(atomlength / zstep / 2.0 + 0.5); // 1/2 length of the template atom image "trace" in the defocus direction to mask "in", in the number of physical z-step units
		index_t karad0 = index_t(atomsizeZ0 / zstep / 2.0 + 0.5); // 1/2 length of an atom image "trace" in the defocus direction to mask when searching for atoms of the same type, in the number of physical z-step units
		index_t karad1 = index_t(atomsizeZ1 / zstep / 2.0 + 0.5); // 1/2 length of an atom image "trace" in the defocus direction to mask when searching for atoms of the next type, in the number of physical z-step units
			printf("\nAtom trace lengths in grid steps: 'in' = %zd, 'out_same' = %zd, 'out_previous' = %zd", 2 * karadt + 1, 2 * karad0 + 1, 2 * karad1 + 1);
		index_t natomtotal = 0; // total number of found atoms

		// allocate storage for detected atom positions
		vector< vector< vector<index_t> > > vvvatompos(natomtypes); // positions of all atoms
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			vvvatompos[nat].resize(natom[nat]); // vvvatompos[nat] is a nat-size vector of vatompos
			for (index_t na = 0; na < natom[nat]; na++) vvvatompos[nat][na].resize(4); // each vector vvvatompos[nat][na] stores (k,j,i) indexes of the position of one atom
		}

		// make a vector of atom type names (extracted from input 1-atom defocus series file names and used for output only)
		vector<string> strAtomNames(natomtypes); // array of atom type names
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			index_t ii0 = filenamebaseIn2[nat].rfind("\\") + 1;
			index_t ii1 = filenamebaseIn2[nat].find('.', ii0);
			strAtomNames[nat] = filenamebaseIn2[nat].substr(ii0, ii1 - ii0);
		}
		
		//**************************************************** start searching for atom positions

		for (index_t nat = 0; nat < natomtypes; nat++) // the cycle over the atom type
		{
			// load the test 3D defocus series array
			XArray3D<float> aaa;
			XArray3DMove<float> aaamove(aaa); // the associated XArray class for applying masks to aaa later
			printf("\n\nNow searching for atoms type no. %zd (%s) ...", nat + 1, strAtomNames[nat].c_str());
			printf("\nReading sample defocus series files %s ...", filenamebaseIn1.c_str());
			IXAHWave2D* ph2new = CreateWavehead2D();
			XArray2D<float> inten;
			inten.SetHeadPtr(ph2new);
			vector<string> infiles;
			FileNames(nangles, ndefocus, filenamebaseIn1, infiles);
			for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
			{
				for (index_t kk = 0; kk < ndefocus; kk++)
				{
					XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
					if (nn == 0 && kk == 0)
					{
						nx = inten.GetDim2();
						ny = inten.GetDim1();
						xmin = GetXlo(inten);
						ymin = GetYlo(inten);
						xstep = GetXStep(inten);
						ystep = GetYStep(inten);
						jarad = index_t(atomsize / ystep / 2.0 + 0.5);
						iarad = index_t(atomsize / xstep / 2.0 + 0.5);
						nx2 = nx / 2 + 1;
						aaa.Resize(nz, ny, nx, 0.0f);
						IXAHWave3D* ph3new = CreateWavehead3D();
						ph3new->SetData(wl, zmin, zmax, ymin, ymin + ystep * ny, xmin, xmin + xstep * nx);
						aaa.SetHeadPtr(ph3new);
					}
					else
					{
						if (inten.GetDim1() != ny) throw std::runtime_error("different ny dimension in input file");
						if (inten.GetDim2() != nx) throw std::runtime_error("different nx dimension in input file");
					}
					for (index_t jj = 0; jj < ny; jj++)
						for (index_t ii = 0; ii < nx; ii++)
						{
							aaa[kk][jj][ii] = inten[jj][ii];
							//aaa[kk][jj][ii] = inten[jj][ii] - finten0;
							//aaa[kk][jj][ii] = ::fabs(aaa[kk][jj][ii]);
						}
				}
			}
			printf("\nSize of input images: (nx,ny,nz) = (%zd, %zd, %zd); minimums = (%g, %g, %g); steps = (%g, %g, %g).", nx, ny, nz, xmin, ymin, zmin, xstep, ystep, zstep);

			// subtract the background, take modulus, then subtrac the mean
			aaa -= finten0;
			aaa.Abs();
			aaa -= (float)aaa.Norm(eNormAver);

			// mask with zeros the vicinity of locations of previously found atoms of other types
			for (index_t natprev = 0; natprev < nat; natprev++)
				for (index_t na = 0; na < natom[natprev]; na++)
					aaamove.FillCylinderPeriodic(vvvatompos[natprev][na][0], vvvatompos[natprev][na][1], vvvatompos[natprev][na][2], karad1, jarad, iarad, 0.0f);

			// optional auxilliary data output
			if (iCorrArrayOut == 1 && nat == natomtypes - 1) // output the masked 1st input array
			{
				printf("\nWriting masked 1st input array in output files %s ...", filenamebaseOut.c_str());
				FileNames(nangles, ndefocus, filenamebaseOut, infiles);
				for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
				{
					for (index_t kk = 0; kk < ndefocus; kk++)
					{
						for (index_t jj = 0; jj < ny; jj++)
							for (index_t ii = 0; ii < nx; ii++)
								inten[jj][ii] = aaa[kk][jj][ii];
						XArData::WriteFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), xar::eGRDBIN);
					}
				}
			}

			// create 3D array of local norms of the test array
			XArray3D<float> abL2(nz, ny, nx, std::numeric_limits<float>::max()); // 3D array of L2 norms of template-size subarrays of aaa; will be later multiplied by L2 norm of the template array
			abL2.SetHeadPtr(aaa.GetHeadPtr()->Clone());
			{ // create a scope so that the bbb array below would be destroyed immediately after it has been used
				XArray3D<float> bbb(aaa); // temporary 3D array containing squares of elements of aaa
				bbb *= bbb;
				omp_set_num_threads(nThreads);
				#pragma omp parallel default(none) shared(abL2, bbb, karadt, jarad, iarad)
				{
					float bsum;
					#pragma omp for schedule(dynamic) nowait
					for (int k = (int)karadt; k < int(nz - karadt); k++)
						for (index_t j = jarad; j < ny - jarad; j++)
							for (index_t i = iarad; i < nx - iarad; i++)
							{
								bsum = 0.0f;
								for (index_t k1 = k - karadt; k1 <= k + karadt; k1++)
									for (index_t j1 = j - jarad; j1 <= j + jarad; j1++)
										for (index_t i1 = i - iarad; i1 <= i + iarad; i1++)
											bsum += bbb[k1][j1][i1];
								abL2[k][j][i] = sqrt(bsum);
							}
				}
			}

			//allocate space for FFT transform and create FFTW plans
			XArray3D<xar::fcomplex> ccc(nz, ny, nx2);
			Fftwf3drc fftf((int)nz, (int)ny, (int)nx);

			// FFT of test array
			printf("\nFFT of the sample defocus series ...");
			fftf.SetRealXArray3D(aaa);
			fftf.ForwardFFT();

			// store away the result of the FFT of the test array
			fftf.GetComplexXArray3D(ccc);


			// load 3D template array (defocus series of a single atom)
			aaa.Fill(0.0);
			printf("\nReading single atom defocus series files %s ...", filenamebaseIn2[nat].c_str());
			FileNames(nangles, ndefocus, filenamebaseIn2[nat], infiles);
			for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
			{
				for (index_t kk = 0; kk < ndefocus; kk++)
				{
					XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
					if (inten.GetDim1() != ny) throw std::runtime_error("different ny dimension in input file");
					if (inten.GetDim2() != nx) throw std::runtime_error("different nx dimension in input file");
					for (index_t jj = 0; jj < ny; jj++)
						for (index_t ii = 0; ii < nx; ii++)
						{
							aaa[kk][jj][ii] = inten[jj][ii];
							//aaa[kk][jj][ii] = inten[jj][ii] - finten0;
							//aaa[kk][jj][ii] = ::fabs(aaa[kk][jj][ii]);
						}
				}
			}

			// subtract the background, take modulus, then subtrac the mean
			aaa -= finten0;
			aaa.Abs();
			aaa -= (float)aaa.Norm(eNormAver);

			// find the centre of gravity of the 3D template array, i.e. the position of the template atom
			double integ = 0.0, xpos = 0.0, ypos = 0.0, zpos = 0.0, dtemp;
			for (index_t kk = 0; kk < ndefocus; kk++)
				for (index_t jj = 0; jj < ny; jj++)
					for (index_t ii = 0; ii < nx; ii++)
					{
						dtemp = abs(aaa[kk][jj][ii]);
						integ += dtemp;
						xpos += dtemp * ii;
						ypos += dtemp * jj;
						zpos += dtemp * kk;
					}
			xpos /= integ; ypos /= integ; zpos /= integ;
			index_t ipos2 = index_t(xpos + 0.5), jpos2 = index_t(ypos + 0.5), kpos2 = index_t(zpos + 0.5);
			double xpos2 = xmin + xstep * ipos2, ypos2 = ymin + ystep * jpos2, zpos2 = zmin + zstep * kpos2;
			printf("\nCentre of mass position of single atom array in pixels = (%zd, %zd, %zd), and in physical units = (%g, %g, %g).", ipos2, jpos2, kpos2, xpos2, ypos2, zpos2);
			xpos2 = rpos2[nat][0]; ypos2 = rpos2[nat][1]; zpos2 = rpos2[nat][2];
			ipos2 = index_t((xpos2 - xmin) / xstep + 0.5); jpos2 = index_t((ypos2 - ymin) / ystep + 0.5); kpos2 = index_t((zpos2 - zmin) / zstep + 0.5);
			printf("\nPosition of this single atom in the parameter file in pixels = (%zd, %zd, %zd), and in physical units = (%g, %g, %g).", ipos2, jpos2, kpos2, xpos2, ypos2, zpos2);
			printf("\nDimensions of the 3D template array in pixels are: nx = %zd, ny = %zd, nz = %zd.", iarad * 2 + 1, jarad * 2 + 1, karadt * 2 + 1);

			// set to zero the values of all pixels outside atomsize vicinity of the centre of mass of the template 1-atom pattern
			aaamove.FillRectComplementPeriodic(kpos2, jpos2, ipos2, karadt, jarad, iarad, 0.0f);

			// optional auxilliary data output
			if (iCorrArrayOut == 2 && nat == natomtypes - 1) // output the masked template array
			{
				printf("\nWriting masked 2nd input array in output files %s ...", filenamebaseOut.c_str());
				FileNames(nangles, ndefocus, filenamebaseOut, infiles);
				for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
				{
					for (index_t kk = 0; kk < ndefocus; kk++)
					{
						for (index_t jj = 0; jj < ny; jj++)
							for (index_t ii = 0; ii < nx; ii++)
								inten[jj][ii] = aaa[kk][jj][ii];
						XArData::WriteFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), xar::eGRDBIN);
					}
				}
			}

			// multiply the L2 of subarrays of the test array by the L2 norm of the template array
			abL2 *= (float)aaa.Norm(eNormL2);

			// FFT of the template array
			printf("\nFFT of the single atom defocus series ...");
			fftf.SetRealXArray3D(aaa);
			fftf.ForwardFFT();

			/// multiply FFTs of 2 arrays, taking the conjugate of the second one
			printf("\nMultiplying FFT of the first 3D array by the conjugate of the FFT of the second ...");
			printf("\nFourier space high-pass filter radius = %g", sqrt(iHPathRad2));
			float ftemp;
			fftwf_complex* pout = fftf.GetComplex();
			int m = 0;
			for (index_t k = 0; k < nz; k++)
				for (index_t j = 0; j < ny; j++)
					for (index_t i = 0; i < nx2; i++)
					{
						if (k * k + j * j + i * i < iHPathRad2)
						{
							pout[m][1] = pout[m][0] = 0.0f; // high-frequency bandpath filter of the correlation array
						}
						else
						{
							ftemp = pout[m][0] * ccc[k][j][i].real() + pout[m][1] * ccc[k][j][i].imag();
							pout[m][1] = -pout[m][1] * ccc[k][j][i].real() + pout[m][0] * ccc[k][j][i].imag();
							pout[m][0] = ftemp;
							//@@@@@ normalization by the modulus (leaving the phase only)
							//ftemp = sqrt(pout[m][0] * pout[m][0] + pout[m][1] * pout[m][1]);
							//pout[m][1] /= ftemp;
							//pout[m][0] /= ftemp;
						}
						m++;
					}

			// inverse FFT of the product
			printf("\nInverse FFT of the product ...");
			fftf.InverseFFT();

			// obtain the correlation array of the test array with the sliding template array
			fftf.GetRealXArray3D(aaa);

			// calculate the array of (local) correlation coefficients
			int kk, jj;
			for (index_t k = 0; k < nz; k++)
			{
				kk = nmodm(int(kpos2 + k), double(nz));
				for (index_t j = 0; j < ny; j++)
				{
					jj = nmodm(int(jpos2 + j), double(ny));
					for (index_t i = 0; i < nx; i++)
						aaa[k][j][i] /= abL2[kk][jj][nmodm(int(ipos2 + i), double(nx))];
				}
			}

			// optional auxilliary data output
			if (iCorrArrayOut == 3 && nat == 0 /*natomtypes - 1*/) // output the 3D correlation distribution array
			{
				printf("\nWriting 3D correlation array in output files %s ...", filenamebaseOut.c_str());
				FileNames(nangles, ndefocus, filenamebaseOut, infiles);
				for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
				{
					for (index_t kk = 0; kk < ndefocus; kk++)
					{
						for (index_t jj = 0; jj < ny; jj++)
							for (index_t ii = 0; ii < nx; ii++)
								inten[jj][ii] = aaa[kk][jj][ii];
						XArData::WriteFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), xar::eGRDBIN);
					}
				}
			}

			// find the maximums
			printf("\nFinding maximums in the correlation array ...");
			index_t kmax = 0, jmax = 0, imax = 0;
			for (index_t na = 0; na < natom[nat]; na++)
			{
				natomtotal++;
				float amax = aaa.Max3D(kmax, jmax, imax);

				vvvatompos[nat][na][0] = nmodm(int(kpos2 + kmax), double(nz)); // absolute z position of the located atom
				vvvatompos[nat][na][1] = nmodm(int(jpos2 + jmax), double(ny)); // absolute y position of the located atom
				vvvatompos[nat][na][2] = nmodm(int(ipos2 + imax), double(nx)); // absolute x position of the located atom
				vvvatompos[nat][na][3] = index_t((amax * 1.0e+8) + 0.5); // "correlation coefficient" of the located atom

				double xmaxA = xmin + vvvatompos[nat][na][2] * xstep;
				double ymaxA = ymin + vvvatompos[nat][na][1] * ystep;
				double zmaxA = zmin + vvvatompos[nat][na][0] * zstep;

				printf("\nAtom type %zd, atom number %zd:", nat + 1, na + 1);
				printf("\nOptimal shift (i,j,k) of the 2nd array to the 1st one in pixels = (%zd, %zd, %zd).", imax, jmax, kmax);
				printf("\nAbsolute position (x,y,z) of the detected atom in physical units = (%g, %g, %g).", xmaxA, ymaxA, zmaxA);
				printf("\nCorrelation coefficient = %g.", amax);

				// fill the atomsize vicinity of the found maximum by zeros, in order to make possible the search for the next largest maximum
				if (na < natom[nat] - 1) 
					aaamove.FillCylinderPeriodic(kmax, jmax, imax, karad0, jarad, iarad, 0.0f);
			}
		} // end of cycle over different atom types


		// bubble-sort the found locations of atoms of each type separately, in accordance with the increasing z-coordinate 
		// (this is done just for increased convenience of checking the output data manually later on)
		vector<index_t> vtemp(3);
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			index_t n = vvvatompos[nat].size();
			for (index_t i = 0; i < n - 1; i++)
				// Last i elements are already in place    
				for (index_t j = 0; j < n - i - 1; j++)
					if (vvvatompos[nat][j][0] > vvvatompos[nat][j + 1][0])
						std::swap<vector<index_t> >(vvvatompos[nat][j], vvvatompos[nat][j + 1]);
		}

		// write the locations of the detected atoms into an XYZ file compatible with Vesta XYZ format
		FILE* ff = fopen(filenameOutXYZ.c_str(), "wt");
		printf("\n\nWriting output file %s in Vesta XYZ format ...\n", filenameOutXYZ.c_str());
		fprintf(ff, "%zd\n", natomtotal); // number of detected atoms
		fprintf(ff, "%s\n", "Atom positions detected by BigBangCT"); // free-form file info line
		for (index_t nat = 0; nat < natomtypes; nat++)
			for (index_t na = 0; na < natom[nat]; na++)
				fprintf(ff, "%s %f %f %f %f\n", strAtomNames[nat].c_str(), xmin + vvvatompos[nat][na][2] * xstep, ymin + vvvatompos[nat][na][1] * ystep, zmin + vvvatompos[nat][na][0] * zstep, vvvatompos[nat][na][3] / 1.0e+8);
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


#endif // ifdef CORRELATION_BASED_METHOD