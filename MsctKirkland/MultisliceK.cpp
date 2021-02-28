// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <thread>
#include <chrono>

#include "XA_ini.h"
#include "XA_file.h"
#include "autosliccmd.h"

using namespace xar;

int Counter_Obj::counter = 0; // active worker thread counter (this main thread is not counted)
bool Counter_Obj::isUpdated = false; // update status of the counter object

int main(void)
{
#ifdef TEG_MULTITHREADED
	Counter_Obj thread_counter; // increments the thread counter on construction and decrements it on destruction
#endif // TEG_MULTITHREADED
	vector<string> autoslictxt(30); // 30 is the current number of input parameters; if it is changed, the corresponding changes need to be applied in autosliccmd.cpp too.

	printf("\nStarting MsctKirkland program ...");
	try
	{
		// read input parameter file
		constexpr char iparfile[] = "MsctKirkland.txt";
		FILE* ff0 = fopen(iparfile, "rt");
		if (!ff0) throw std::exception((std::string("Error opening parameter file ") + iparfile + ".").c_str());
		else printf("\nReading input parameters from %s file ...", iparfile);

		char cline[1024], ctitle[1024], cparam[1024];
	
		// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
		while (true)
		{
			fgets(cline, 1024, ff0); 
			if (!(cline[0] == '/' && cline[1] == '/')) break;
		}

		strtok(cline, "\n"); // 1st parameter: Input_file_with_atomic_numbers_and_coordinates_in_XYZ_format
		autoslictxt[0] = cline;
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2nd parameter: Output_files_shall_contain_intensity(0),_phase(1),_complex_amplitude(2)_or_3D_potential(3)
		autoslictxt[27] = cline; // The numbering of these parameters is 'historic', it can be changed, but then the corresponding changes need to be applied in autosliccmd.cpp too.
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output type from the input parameter file.");
		int nOutputType = atoi(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3rd parameter: Output_TIFF/GRD/GRC_filename_template
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file template from the input parameter file.");
		string outfilename(cparam); // output filename stub
		if ((GetFileExtension(outfilename) != string(".TIFF")) && (GetFileExtension(outfilename) != string(".TIF")) && (GetFileExtension(outfilename) != string(".GRD")) && (GetFileExtension(outfilename) != string(".GRC")))
			throw std::exception("Error: output filename extension must be TIF, GRD or GRC.");
		if ((nOutputType == 0 || nOutputType == 1 || nOutputType == 3) && !(GetFileExtension(outfilename) == string(".GRD") || GetFileExtension(outfilename) == string(".TIF") || GetFileExtension(outfilename) == string(".TIFF")))
			throw std::exception("Output type (parameter 2) is inconsistent with output filename extension in the input parameter file.");
		if (nOutputType == 2 && (GetFileExtension(outfilename) != string(".GRC")))
			throw std::exception("Output type (parameter 2) is inconsistent with output filename extension in the input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4th paramter: Use_multislice(0),_projection(1),_or_1st_Born(2)_approximation
		autoslictxt[25] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading calculation mode (parameter 4) from the input parameter file.");
		int nCalculationMode = atoi(cparam);
		if ((nCalculationMode != 0) && (nCalculationMode != 1) && (nCalculationMode != 2))
			throw std::exception("Unknown calculation mode (parameter 4) in the input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5th paramter: Incident__electron_beam_energy_in_keV
		autoslictxt[10] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading electron beam energy from the input parameter file.");
		double ev = atof(cparam) * 1000; // accelerating voltage in eV
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6th parameter: Wavefunction_size_in_pixels,_Nx,Ny
		autoslictxt[11] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7th parameter: Slice_thickness_in_Angstroms
		autoslictxt[13] = cline;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8th parameter:Objective_aperture_in_mrad
		autoslictxt[7] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9th parameter: Spherical_aberration_Cs3_and_Cs5_in_mm
		autoslictxt[5] = cline;
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10th paramter: Include_thermal_vibrations(1)_or_not(0)
		autoslictxt[17] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading thermal vibration option from the input parameter file.");
		int iThermalVibr = atoi(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11th parameter: ____Temperature_in_degrees_K
		autoslictxt[18] = cline;
		double dTemperature(-1.0);
		if (iThermalVibr != 0)
		{
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading temperature from the input parameter file.");
			dTemperature = atof(cparam);
			if (dTemperature < 0) throw std::exception("Temperature in degrees K cannot be negative.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12th parameter: ____Number_of_configurations_to_average_over
		autoslictxt[19] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading the number of configurations parameter from the input parameter file.");
		int iNumConfig = atoi(cparam);
		if (iNumConfig < 0) throw std::exception("Number of configurations cannot be negative.");
		if (iThermalVibr != 0 && iNumConfig > 1 && !(nOutputType == 0 || nOutputType == 3))
			throw std::exception("Output type can only be 0 (intensity) or 3 (3D potential) when the thermal vibrations option with multiple configurations is enabled in the input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13th parameter: Ice_layer_thickness_in_Angstroms
		autoslictxt[29] = cline;
		if (iThermalVibr != 0 && dTemperature > 273.15)
		{
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading ice layer thickness from the input parameter file.");
			if (atof(cparam) > 0) throw std::exception("Temperature cannot be higher than 273.15K in the presence of an ice layer.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14th parameter: Text_file_with_output_rotation_angles_in_degrees_and_defocus_distances_in_Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading the name of the file with defocus distances and rotational positions from the input parameter file.");
		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;
		ReadDefocusParamsFile(string(cparam), v2angles, vvdefocus);
		index_t nangles = v2angles.size(); // number of rotation steps 
		vector<index_t> vndefocus(nangles); // vector of numbers of defocus planes at different rotation angles
		for (index_t i = 0; i < nangles; i++) vndefocus[i] = vvdefocus[i].size();
		if (GetFileExtension(outfilename) == string(".GRC")) // check that there is only one input defocused complex amplitude file per each illumination angle
		{
			for (index_t i = 0; i < nangles; i++)
				if (vndefocus[i] != 1)
					throw std::exception("Error: only one output defocused complex amplitude file per each illumination direction is allowed.");
		}
		vector<string> voutfilenamesTot;
		FileNames2(vndefocus, outfilename, voutfilenamesTot); // create "2D array" of output filenames

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15th parameter: Number_of_worker_threads_to_launch_in_CT_simulation_mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading the number of worker threads from the input parameter file.");
		unsigned int ncores = (unsigned int)atoi(cparam) + 1; // number of threads to use (expected to be equal to the number of cores) 

		fclose(ff0); // close input parameter file

		//!!The following additional parameters are required in autosliccmd.cpp, but are not currently used by this program
		autoslictxt[1] = "2.Replicate_unit_cell_by_NCELLX,NCELLY,NCELLZ: 1 1 1";
		autoslictxt[2] = ""; // this parameter is not used any more, output filenames are passed as a separate argument vstrfileout
		autoslictxt[3] = "4.Do_you_want_to_include_partial_coherence: 0";
		autoslictxt[4] = "5.____Illumination_angle_min,_max_in_mrad: 0.0 0.0";
		autoslictxt[6] = "7.____Defocus_mean,_standard_deviation,_and_sampling_size_in_Angstroms: 0.0 0.0 0.0";
		autoslictxt[8] = "9.Do_you_want_to_start_from_previous_result: 0";
		autoslictxt[9] = "10.____Name_of_file_to_start_from: 0";
		autoslictxt[12] = "13.Crystal_tilt_x,y_in_mrad: 0.0 0.0";
		autoslictxt[14] = "15.Do_you_want_to_record_the_(real,imag)_value_of_selected_beams_vs._thickness: 0";
		autoslictxt[15] = "16.____Name_of_file_for_beams_info: 0";
		autoslictxt[16] = "17.____Number_of_beams: 0";
		autoslictxt[20] = "21.____Initial_seed_for_random_number_generator: 1";
		autoslictxt[21] = "22.Do_you_want_to_output_intensity_vs._depth_cross_section: 0";
		autoslictxt[22] = "23.____Name_of_file_to_get_depth_profile_image: 0";
		autoslictxt[23] = "24.____y_position_of_depth_cross_section_in_Angstroms: 0.0";
		autoslictxt[26] = ""; // this parameter is not used any more, defocus values are passed as a separate argument vdefocus		

		constexpr double hp = 6.62607004e-34; // Planck's constant (m2 kg / s)
		constexpr double cc = 299792458; // speed of light (m / s)
		constexpr double ee = 1.602176634e-19; // electron charge (coulomb)
		constexpr double m0 = 9.1093837015e-31; // electron rest mass (kg)
		double ewl = hp * cc / sqrt(ee * ev * (2.0 * m0 * cc * cc + ee * ev));
		printf("\nElectron wavelength = %f (pm)", ewl * 1.0e+12);

		// find out the number of CPU cores available in the computer
		//unsigned int ncores = std::thread::hardware_concurrency();
		//printf("\nNumber of CPU cores detected: %d", ncores);

		// start the execution timer
		std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
	
		// start the cycle over projection angles
		char bufangle[1024];
		index_t ndefcurrent(0);
		for (size_t i = 0; i < nangles; i++)
		{
			Pair angle = v2angles[i];
			printf("\nIllumination angle: y = %g, x' = %g (degrees)", angle.a, angle.b);
			sprintf(bufangle, "%f %f", angle.a * PI180, angle.b * PI180);
			autoslictxt[24] = "25.Sample_Y_and_X'_rotation_angles_in_radians: " + string(bufangle);

			// start the cycle over defocus distances (we only create output file names in this inner cycle)
			index_t ndefocus = vndefocus[i]; // number of defocus planes at the current illumination angle
			vector<Pair> vdefocus = vvdefocus[i]; // vector of Z" angles and defocus distances at the current illumination angle
			// Note that vdefocus[n].a should contain Z" angle in degrees, because xar::XArray2DSpln<T>::Rotate function expects the rotation angle in degrees
			vector<string> vstrfileout(ndefocus); // vector of output filenames at the current illumination angle
			for (index_t n = 0; n < ndefocus; n++) vstrfileout[n] = voutfilenamesTot[ndefcurrent++];


			//Here we call Kirkland's autoslic at each angle
			autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 1"; // the first thread must initialize the FFTW plan, subsequent ones can copy it
#ifdef TEG_MULTITHREADED
			// NOTE that in this multithreaded model, exceptions thrown by worker threads will NOT be caught in the master thread
			if (i > 0) autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 0";
			thread_counter.SetUpdated(false);
			std::thread threadObj(autosliccmd, autoslictxt, vdefocus, vstrfileout);
			if (i == 0) threadObj.join(); // we need to let the first worker thread finish execution, so that it can create the FFTW "plan" to be shared  with other threads
			if (threadObj.joinable()) threadObj.detach(); // if we don't do this, threadObj will call Terminate() on the attached thread when the threadObj goes out of scope
			while (thread_counter.GetCount() >= 0 && (!thread_counter.GetUpdated() || thread_counter.GetCount() >= (int)ncores))
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
