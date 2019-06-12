// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#define TEG_MULTISLICE //if undefined, Kirkland's multislice code is used instead

#include <stdio.h>
#include "XA_ini.h"
#include "XA_data.h"
#include "XAHWave.h"

#ifdef TEG_MULTISLICE
#include "ProjSpheres.h"
#else
#include "autosliccmd.h"
#endif // TEG_MULTISLICE

using namespace xar;

int main(void)
{
	string outfilename("C:\\Users\\tgureyev\\Downloads\\aaa.grd");
	const size_t nangles = 1;

#ifdef TEG_MULTISLICE
	const size_t nslices = 1; 
	const size_t nx(512), ny(512);
	vector<double> vHead(5); //analogue of Wavehead2D
	const double energ(1.0); //E in keV
	const double LengthScale(1.0); //scaling factor for length-type parameters
	vHead[0] = 12.398E-4 / energ; //wl in microns
	vHead[1] = -LengthScale; //ylo
	vHead[2] = LengthScale; //yhi
	vHead[3] = -LengthScale; //xlo
	vHead[4] = LengthScale; //xhi
	const double zlo(-LengthScale), zhi(LengthScale);

	const size_t nSpheres(10);
	XArray1D<double> R_in(nSpheres), x_in(nSpheres), y_in(nSpheres), z_in(nSpheres), x_in_i(nSpheres), z_in_i(nSpheres);

	R_in[0] = 0.3 * LengthScale;
	for (size_t i = 1; i < nSpheres; i++) R_in[i] = 0.1 * LengthScale;
	
	x_in[0] = 0;
	for (size_t i = 1; i < nSpheres; i++) x_in[i] = pow(-1, i) * double(i % 3 + 1) * 0.2 * LengthScale;

	y_in[0] = -0.5 * LengthScale;
	for (size_t i = 1; i < nSpheres; i++) y_in[i] = 0.5 * LengthScale;

	z_in[0] = 0;
	for (size_t i = 1; i < nSpheres; i++) z_in[i] = pow(-1, int(i / 5)) * double(i % 3 + 1) * 0.2 * LengthScale;

	vector<dcomplex> nc(nSpheres);
	vector<double> beta_in(nSpheres), delta_in(nSpheres);
	double BetaOrder = 1.E-5, DeltaOrder = 1.e-6; //scaling factors for beta and delta
	for (size_t i = 0; i < nSpheres; i++) beta_in[i] = double(i % 5 + 1) * BetaOrder;
	for (size_t i = 0; i < nSpheres; i++) delta_in[i] = double(i % 5 + 1) * DeltaOrder;
	for (size_t i = 0; i < nSpheres; i++)
		nc[i] = dcomplex(1.0 - delta_in[i], beta_in[i]);
#endif

	size_t i_dot = outfilename.rfind('.');
	size_t nfield_length = (nangles == 1) ? 1 : 1 + size_t(log10(double(nangles - 1))); //maximum number of digits in the output file name
	char ndig[8];
	sprintf(ndig, "%zd", nfield_length); //convert the calculated maximum number of digits into a string, e.g. 3 into "3"
	string myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded numbers into file names - see usage below

	double angle;
	char buffer[128], bufangle[128];
	string outfilename_i;
#ifndef TEG_MULTISLICE
	string autoslictxt[numaslicpars], strAngle;
#endif
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
#ifdef TEG_MULTISLICE
			for (size_t k = 0; k < nSpheres; k++)
			{
				x_in_i[k] = x_in[k] * cos(angle) + z_in[k] * sin(angle); //x - position of the centre of the internal sphere at the i - th rotation step
				z_in_i[k] = -x_in[k] * sin(angle) + z_in[k] * cos(angle); //z - position of the centre of the internal sphere at the i - th rotation step
			}
			MultisliceSphereNF(nx, ny, vHead, nc, R_in, x_in_i, y_in, z_in_i, zlo, zhi, nslices, outfilename_i.data());
#else
			//!!! Now instead of TEG function MultisliceSphereNF() we will call Kirkland's autoslic at each angle
			//!!! One has to define at least parameters no. 0, 2, 10, 11, 13, 25 and 26 in the list below
			autoslictxt[0] = "1.Name_of_file_with_input_atomic_coordinates_in_x,y,z_format: 3j6kLysLesKirck.xyz";
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
			autoslictxt[11] = "12.Wavefunction_size_in_pixels,_Nx,Ny: 1024 1024";
			autoslictxt[12] = "13.Crystal_tilt_x,y_in_mrad: 0.0 0.0";
			autoslictxt[13] = "14.Slice_thickness_in_Angstroms: 15.0";
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
			autoslictxt[25] = "26.Coordinates_of_the_centre_of_rotation_in_Angstroms: 56.526327 56.526327";
			autosliccmd(autoslictxt);
#endif
		}
	}
	catch (std::exception& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}

	return 0;
}
