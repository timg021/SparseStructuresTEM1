// 1atomwiseMSimage.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <thread>
#include <chrono>

#include "XA_ini.h"
#include "XArray2D.h"
#include "XA_data.h"

#include "autosliccmd.h"
#undef TEG_MULTITHREADED // this program is intended to run single-threaded only

#include "pdb.h"
extern "C" int pdb_read(pdbdata* ppd, int nfiletype, char* pdbfile);
extern "C" int xyz_Kirck_write1(int i, pdbdata* ppd, char* outfile, char* cfileinfo, double ctblength, double xminx0, double yminy0, double zminz0);

extern xar::XArray2D<float> intenTot; // accumulated intensity data - created in 1autosliccmd.cpp, wrtten into a file at the end in this module

using namespace xar;

int Counter_Obj::counter = 0; // active worker thread counter (this main thread is not counted)
bool Counter_Obj::isUpdated = false; // update status of the counter object

int main(void)
{
#ifdef TEG_MULTITHREADED
	Counter_Obj thread_counter; // increments the thread counter on construction and decrements it on destruction
#endif // TEG_MULTITHREADED
	vector<string> autoslictxt(29); // 29 is the current number of input parameters; if it is changed, the corresponding changes need to be applied in autosliccmd.cpp too.

	try
	{
		printf("\nStarting 1atomwiseMSimage program ...");
		// read input parameter file
		FILE* ff0 = fopen("MsctKirkland.txt", "rt");
			if (!ff0) throw std::exception("Error: cannot open parameter file MsctKirkland.txt.");
			char cline[1024], ctitle[1024], cparam[1024];
			// The ordering of these parameters is 'historic', it can be changed, but then the corresponding changes need to be applied in autosliccmd.cpp too.
			fgets(cline, 1024, ff0); // 1st line - comment
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2nd line: Input_file_with_atomic_numbers_and_coordinates_in_XYZ_format
			autoslictxt[0] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3rd line: Output_GRD/GRC_filename
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 3 of input parameter file.");
			string outfilename(cparam); // output filename stub
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4th line: Output_intensity(0),_phase(1)_or_complex_amplitude(2)
			autoslictxt[27] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5th line: Use_multislice(0),_projection(1)_or_1st_Born(2)_approximation
			autoslictxt[25] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6th line: Incident__electron_beam_energy_in_keV
			autoslictxt[10] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7th line: Wavefunction_size_in_pixels,_Nx,Ny
			autoslictxt[11] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8th line: Slice_thickness_in_Angstroms
			autoslictxt[13] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9th line: Propagation(defocus)_distance_for_exit_wave_in_Angstroms !!!!!!!!!!!!!! CHANGE
			autoslictxt[26] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10th line: Total_CT_rotation_span_in_degrees
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 10 of input parameter file.");
			double angle_max = atof(cparam); // total rotation span in degrees 
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11th line: Number_of_CT_rotation_angles
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 11 of input parameter file.");
			size_t nangles = (size_t)atoi(cparam); // total rotation span in degrees 
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12th line: Number_of_worker_threads_to_launch_in_CT_simulation_mode
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 12 of input parameter file.");
			unsigned int ncores = (unsigned int)atoi(cparam) + 1; // number of threads to use (expected to be equal to the number of cores) 
			fgets(cline, 1024, ff0); // 13st line - comment
			fgets(cline, 1024, ff0); // 14st line - comment
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15th line: Replicate_unit_cell_by_NCELLX,NCELLY,NCELLZ
			autoslictxt[1] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 16th line: Do_you_want_to_include_partial_coherence
			autoslictxt[3] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17th line: ____Illumination_angle_min,_max_in_mrad
			autoslictxt[4] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 18th line: ____Spherical_aberration_Cs3,_Cs5_in_mm
			autoslictxt[5] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19th line: ____Defocus_mean,_standard_deviation,_and_sampling_size_in_Angstroms
			autoslictxt[6] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 20th line: ____Objective_aperture_in_mrad
			autoslictxt[7] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 21th line: Do_you_want_to_start_from_previous_result
			autoslictxt[8] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 22nd line: ____Name_of_file_to_start_from
			autoslictxt[9] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 23nd line: Crystal_tilt_x,y_in_mrad
			autoslictxt[12] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 24th line: Do_you_want_to_record_the_(real,imag)_value_of_selected_beams_vs._thickness
			autoslictxt[14] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 25th line: ____Name_of_file_for_beams_info
			autoslictxt[15] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 26th line: ____Number_of_beams
			autoslictxt[16] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 27th line: Do_you_want_to_include_thermal_vibrations
			autoslictxt[17] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 28th line: ____Temperature_in_degrees_K
			autoslictxt[18] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 29th line: ____Number_of_configurations_to_average_over
			autoslictxt[19] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 30th line: ____Initial_seed_for_random_number_generator
			autoslictxt[20] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 31st line: Do_you_want_to_output_intensity_vs._depth_cross_section
			autoslictxt[21] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 32nd line: ____Name_of_file_to_get_depth_profile_image
			autoslictxt[22] = cline;
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 33rd line: ____Y_position_of_depth_cross_section_in_Angstroms
			autoslictxt[23] = cline;
		fclose(ff0); // close input parameter file

		// find out the number of CPU cores available in the computer
		//unsigned int ncores = std::thread::hardware_concurrency();
		//printf("\nNumber of CPU cores detected: %d", ncores);

		char inpdbfile[256]; // input file name
		char infiletype[256]; // input file type
		char strctblength[256]; // CT qube side length
		char outfile[256]; // output file name
		char cfileinfo[256]; // freeform info line

		// read PDB file translate parameter file
		FILE* ffpar = fopen("pdb.txt", "rt");
			if (!ffpar) throw std::exception("Error: cannot open parameter file pdb.txt.");
			fgets(inpdbfile, 256, ffpar); strtok(inpdbfile, "\n");
			fgets(infiletype, 256, ffpar); strtok(infiletype, "\n");
			fgets(strctblength, 256, ffpar); strtok(strctblength, "\n");
			fgets(outfile, 256, ffpar); strtok(outfile, "\n");
			fgets(cfileinfo, 256, ffpar); strtok(cfileinfo, "\n");
		fclose(ffpar);

		int nfiletype = 0;
		nfiletype = (int)atof(infiletype);
		double ctblength = 0;
		ctblength = atof(strctblength);

		// read PDB or VestaXYZ file to get atomic coordinates
		pdbdata pd;
		pdb_read(&pd, nfiletype, inpdbfile);
		double xmin, xmax, ymin, ymax, zmin, zmax;
		xmin = xmax = pd.adata[0].x;
		ymin = ymax = pd.adata[0].y;
		zmin = zmax = pd.adata[0].z;
		for (size_t i = 1; i < pd.natoms; i++)
		{
			if (pd.adata[i].x < xmin) xmin = pd.adata[i].x;
			else if (pd.adata[i].x > xmax) xmax = pd.adata[i].x;
			if (pd.adata[i].y < ymin) ymin = pd.adata[i].y;
			else if (pd.adata[i].y > ymax) ymax = pd.adata[i].y;
			if (pd.adata[i].z < zmin) zmin = pd.adata[i].z;
			else if (pd.adata[i].z > zmax) zmax = pd.adata[i].z;
		}
		double xminx0 = 0.5 * ctblength - 0.5 * (xmax - xmin) - xmin;
		double yminy0 = 0.5 * ctblength - 0.5 * (ymax - ymin) - ymin;
		double zminz0 = 0.5 * ctblength - 0.5 * (zmax - zmin) - zmin;
		if (xminx0 < -xmin || yminy0 < -ymin || zminz0 < -zmin)
			throw std::exception("Error: sample's x, y or z extent is larger than the defined CT sample qube side length!!!");

		// define output filename
		nangles = pd.natoms; // !!! define "number of angles" equal to natoms
		size_t i_dot = outfilename.rfind('.');
		size_t nfield_length = (nangles == 1) ? 1 : 1 + size_t(log10(double(nangles - 1))); //maximum number of digits in the output file name
		char ndig[8];
		sprintf(ndig, "%zd", nfield_length); //convert the calculated maximum number of digits into a string, e.g. 3 into "3"
		string myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded numbers into file names - see usage below

		char buffer[128], bufangle[128];
		string strAngle, outfilename_i;
		//double angle;
		//double angle_step = angle_max / 180.0 * PI / double(nangles); // rotation step in radians

		// start the execution timer
		std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

		// start the cycle over projection angles
		for (size_t i = 0; i < nangles; i++)
		{
			printf("\nAngle = %zd", i);
			
			//write new XYZ file with a single atom no. i
			xyz_Kirck_write1((int)i, &pd, outfile, cfileinfo, ctblength, xminx0, yminy0, zminz0);

			outfilename_i = outfilename;
			sprintf(buffer, myformat.data(), i);
			outfilename_i.insert(i_dot, buffer);
			//angle = angle_step * double(i);
			sprintf(bufangle, "%f", 0.0); strAngle = bufangle;
			//Here we call Kirkland's autoslic at each angle
			autoslictxt[2] = "3.Name_of_file_to_get_binary_output_of_multislice_result: " + outfilename_i;
			autoslictxt[24] = "25.Sample_(xz)_rotation_angle_in_radians: " + strAngle;
			autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 1"; // the first thread must initialize the FFTW plan, subsequent ones can copy it
#ifdef TEG_MULTITHREADED
			if (i > 0) autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 0";
			thread_counter.SetUpdated(false);
			std::thread threadObj(autosliccmd, autoslictxt);
			if (i == 0) threadObj.join(); // we need to let the first worker thread finish execution, so that it can create the FFTW "plan" to be shared  with other threads
			if (threadObj.joinable()) threadObj.detach(); // if we don't do this, threadObj will call Terminate() on the attached thread when the threadObj goes out of scope
			while (!thread_counter.GetUpdated() || thread_counter.GetCount() >= ncores)
				std::this_thread::sleep_for(std::chrono::milliseconds(10)); // we allow ncores of threads to be launched
#else
			//autosliccmd(autoslictxt); // single-threaded execution mode
#endif // TEG_MULTITHREADED
		}

		// save the accumulated total intensity file
		intenTot += 1.0f;
		XArData::WriteFileGRD(intenTot, outfilename.c_str(), xar::eGRDBIN);

#ifdef TEG_MULTITHREADED
		while (thread_counter.GetCount() > 1)
			std::this_thread::sleep_for(std::chrono::milliseconds(100)); // wait for all worker threads to finish or fail
#endif // TEG_MULTITHREADED
		std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
		printf("\nMain program finished. Execution time = %I64d s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	printf("\nPress any key to exit..."); getchar();
	return 0;
}
