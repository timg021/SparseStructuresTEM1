// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <thread>
#include <chrono>

#include "XA_ini.h"
#include "autosliccmd.h"

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
		// read input parameter file
		FILE* ff0 = fopen("MsctKirkland.txt", "rt");
			if (!ff0) throw std::exception("Error: cannot open parameter file MsctKirkland.txt.");

		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024];
	
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
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9th line: Defocus_distance_MIN,MAX,STEP_in_Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading line 9 of input parameter file.");
		double defocus_min = atof(cparam); // minimum defocus in Angstroms 
		double defocus_max = atof(cparam1); // maximum defocus in Angstroms 
		double defocus_step = atof(cparam2); // defocus step in Angstroms 
		autoslictxt[26] = ""; // this parameter is not used any more, defocus values are passed as a separate argument vdefocus
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10th line:Objective_aperture_in_mrad
		autoslictxt[7] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11th line: Include_thermal_vibrations(1)_or_not(0)
		autoslictxt[17] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12th line: ____Temperature_in_degrees_K
		autoslictxt[18] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13th line: ____Number_of_configurations_to_average_over
		autoslictxt[19] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14th line: ____Initial_seed_for_random_number_generator
		autoslictxt[20] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15th line: Total_CT_rotation_span_in_degrees
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 10 of input parameter file.");
		double angle_max = atof(cparam); // total rotation span in degrees 
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 16th line: Number_of_CT_rotation_angles
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 11 of input parameter file.");
		size_t nangles = (size_t)atoi(cparam); // total rotation span in degrees 
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17th line: Number_of_worker_threads_to_launch_in_CT_simulation_mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading line 12 of input parameter file.");
		unsigned int ncores = (unsigned int)atoi(cparam) + 1; // number of threads to use (expected to be equal to the number of cores) 
		fgets(cline, 1024, ff0); // 18st line - comment
		fgets(cline, 1024, ff0); // 119st line - comment
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 20th line: Replicate_unit_cell_by_NCELLX,NCELLY,NCELLZ
		autoslictxt[1] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 21th line: Do_you_want_to_include_partial_coherence
		autoslictxt[2] = ""; // this parameter is not used any more, output filenames are passed as a separate argument vstrfileout
		autoslictxt[3] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 22th line: ____Illumination_angle_min,_max_in_mrad
		autoslictxt[4] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 23th line: ____Spherical_aberration_Cs3,_Cs5_in_mm
		autoslictxt[5] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 24th line: ____Defocus_mean,_standard_deviation,_and_sampling_size_in_Angstroms
		autoslictxt[6] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 25th line: Do_you_want_to_start_from_previous_result
		autoslictxt[8] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 26nd line: ____Name_of_file_to_start_from
		autoslictxt[9] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 27nd line: Crystal_tilt_x,y_in_mrad
		autoslictxt[12] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 28th line: Do_you_want_to_record_the_(real,imag)_value_of_selected_beams_vs._thickness
		autoslictxt[14] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 29th line: ____Name_of_file_for_beams_info
		autoslictxt[15] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 30th line: ____Number_of_beams
		autoslictxt[16] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 31st line: Do_you_want_to_output_intensity_vs._depth_cross_section
		autoslictxt[21] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 32nd line: ____Name_of_file_to_get_depth_profile_image
		autoslictxt[22] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 33rd line: ____Y_position_of_depth_cross_section_in_Angstroms
		autoslictxt[23] = cline;
		fclose(ff0); // close input parameter file

		char buffer[128], bufangle[128];
		string strAngle, outfilename_i, outfilename_j;
		double angle;
		double angle_step = angle_max / 180.0 * PI / double(nangles); // rotation step in radians
		
		size_t ndefocus = size_t((defocus_max - defocus_min) / defocus_step + 0.5); // number of defocus planes to propagate to at each rotation angle		
		printf("\nNumber of defocus planes = %zd.", ndefocus);
		vector<string> vstrfileout(ndefocus); // vector of full output filenames
		vector<double> vdefocus(ndefocus); // vector of defocus distances
		for (size_t j = 0; j < ndefocus; j++) vdefocus[j] = defocus_min + defocus_step * j;
		
		// create formatting string to add properly formatted indexes at the end of the output file names
		size_t i_dot = outfilename.rfind('.'), nfieldA_length, nfieldB_length;
		char ndig[8];
		string myformat("");
		if (ndefocus > 1)
		{
			nfieldA_length = 1 + size_t(log10(double(ndefocus - 1))); //maximum number of digits corresponding to defocuses in the output file name
			sprintf(ndig, "%zd", nfieldA_length); //convert the calculated maximum number of digits corresponding to defocuses into a string, e.g. 5 into "5"
			myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded defocus indexes into file names
		}
		if (nangles > 1)
		{
			nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the output file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat += "_%0" + string(ndig) + "d"; //construct format string for inserting two 0-padded angle indexes into file names - see usage below
		}
		
		// find out the number of CPU cores available in the computer
		//unsigned int ncores = std::thread::hardware_concurrency();
		//printf("\nNumber of CPU cores detected: %d", ncores);

		// start the execution timer
		std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
	
		// start the cycle over projection angles
		printf("\nStarting TEG MultisliceK program ...");
		for (size_t i = 0; i < nangles; i++)
		{
			printf("\nAngle = %zd", i);
			angle = angle_step * double(i);
			sprintf(bufangle, "%f", angle); strAngle = bufangle;
			autoslictxt[24] = "25.Sample_(xz)_rotation_angle_in_radians: " + strAngle;

			// start the cycle over defocus distances (we only create output file names in this inner cycle)
			for (size_t j = 0; j < ndefocus; j++)
			{
				outfilename_j = outfilename;
				if (ndefocus == 1 && nangles > 1) sprintf(buffer, myformat.data(), i);
				else if (ndefocus > 1 && nangles == 1) sprintf(buffer, myformat.data(), j);
				else sprintf(buffer, myformat.data(), j, i);
				outfilename_j.insert(i_dot, buffer);
				vstrfileout[j] = outfilename_j;
			}

			//Here we call Kirkland's autoslic at each angle
			autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 1"; // the first thread must initialize the FFTW plan, subsequent ones can copy it
#ifdef TEG_MULTITHREADED
			// NOTE that in this multithreaded model, exceptions thrown by worker threads will NOT be caught in the master thread
			if (i > 0) autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 0";
			thread_counter.SetUpdated(false);
			std::thread threadObj(autosliccmd, autoslictxt, vdefocus, vstrfileout);
			if (i == 0) threadObj.join(); // we need to let the first worker thread finish execution, so that it can create the FFTW "plan" to be shared  with other threads
			if (threadObj.joinable()) threadObj.detach(); // if we don't do this, threadObj will call Terminate() on the attached thread when the threadObj goes out of scope
			while (!thread_counter.GetUpdated() || thread_counter.GetCount() >= (int)ncores)
				std::this_thread::sleep_for(std::chrono::milliseconds(10)); // we allow ncores of threads to be launched
			if (thread_counter.GetCount() < 0) // meaning that a thread has requested the whole program to terminate
				throw std::exception("A thread has requested the whole program to terminate.");
#else
			autosliccmd(autoslictxt, vdefocus, vstrfileout); // single-threaded execution mode
#endif // TEG_MULTITHREADED

		}

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
