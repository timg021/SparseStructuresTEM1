// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define TEG_MULTISLICE //if undefined, Kirkland's multislice code is used instead

#include <stdio.h>

#include "XA_ini.h"

#include "autosliccmd.h"

using namespace xar;

int main(void)
{
	string outfilename("C:\\Users\\tgureyev\\Downloads\\C1a.grd");
	const size_t nangles = 1;

	size_t i_dot = outfilename.rfind('.');
	size_t nfield_length = (nangles == 1) ? 1 : 1 + size_t(log10(double(nangles - 1))); //maximum number of digits in the output file name
	char ndig[8];
	sprintf(ndig, "%zd", nfield_length); //convert the calculated maximum number of digits into a string, e.g. 3 into "3"
	string myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded numbers into file names - see usage below
	string autoslictxt[numaslicpars], strAngle;

	double angle;
	char buffer[128], bufangle[128];
	string outfilename_i;
	try
	{
		for (size_t i = 0; i < nangles; i++)
		{
			printf("Angle = %zd\n", i);
			outfilename_i = outfilename;
			sprintf(buffer, myformat.data(), i);
			outfilename_i.insert(i_dot, buffer);
			angle = PI * double(i) / double(nangles);
			sprintf(bufangle, "%f", angle); strAngle = bufangle;
			//Here we call Kirkland's autoslic at each angle
			//!!! One has to define at least parameters no. 0, 2, 10, 11, 13, 25 and 26 in the list below
			autoslictxt[0] = "1.Name_of_file_with_input_atomic_coordinates_in_x,y,z_format: 1Catom.xyz";
			autoslictxt[1] = "2.Replicate_unit_cell_by_NCELLX,NCELLY,NCELLZ: 1 1 1";
			autoslictxt[2] = "3.Name_of_file_to_get_binary_output_of_multislice_result: " + outfilename_i;
			autoslictxt[3] = "4.Do_you_want_to_include_partial_coherence: 0";
			autoslictxt[4] = "5.____Illumination_angle_min,_max_in_mrad: 0.0 0.0";
			autoslictxt[5] = "6.____Spherical_aberration_Cs3,_Cs5_in_mm: 0.0 0.0";
			autoslictxt[6] = "7.____Defocus_mean,_standard_deviation,_and_sampling_size_in_Angstroms: 0.0 0.0 0.0";
			autoslictxt[7] = "8.____Objective_aperture_in_mrad: 0.0";
			autoslictxt[8] = "9.Do_you_want_to_start_from_previous_result: 0";
			autoslictxt[9] = "10.____Name_of_file_to_start_from: 0";
			autoslictxt[10] = "11.Incident_beam_energy_in_keV: 200.0";
			autoslictxt[11] = "12.Wavefunction_size_in_pixels,_Nx,Ny: 512 512";
			autoslictxt[12] = "13.Crystal_tilt_x,y_in_mrad: 0.0 0.0";
			autoslictxt[13] = "14.Slice_thickness_in_Angstroms: 1.0";
			autoslictxt[14] = "15.Do_you_want_to_record_the_(real,imag)_value_of_selected_beams_vs._thickness: 0";
			autoslictxt[15] = "16.____Name_of_file_for_beams_info: 0";
			autoslictxt[16] = "17.____Number_of_beams: 0";
			autoslictxt[17] = "18.Do_you_want_to_include_thermal_vibrations: 0";
			autoslictxt[18] = "19.____Type_the_temperature_in_degrees_K: 0.0";
			autoslictxt[19] = "20.____Type_number_of_configurations_to_average_over: 0";
			autoslictxt[20] = "21.____Type_initial_seed_for_random_number_generator: 0";
			autoslictxt[21] = "22.Do_you_want_to_output_intensity_vs._depth_cross_section: 0";
			autoslictxt[22] = "23.____Type_name_of_file_to_get_depth_profile_image: 0";
			autoslictxt[23] = "24.____Type_y_position_of_depth_cross_section_in_Angstroms: 0.0";
			autoslictxt[24] = "25.Sample_(xz)_rotation_angle_in_radians: " + strAngle;
			autoslictxt[25] = "26.Coordinates_of_the_centre_of_rotation_in_Angstroms: 2.64 2.64";
			autoslictxt[26] = "27.Use_multislice(0),_projection(1)_or_1st_Born(2)_approximation: 1";
			autosliccmd(autoslictxt);
		}
	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	return 0;
}
