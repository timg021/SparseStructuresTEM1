Brief description of the "inverse" simulation (3D object reconstruction) module "PhaseRetrieval.exe"
of the Differential Holographic Tomography (DHT) software package

PhaseRetrieval.exe module performs 3D object reconstruction and reprojection tasks, 
using input data (defocused images) that have been generated by the MsctKirkland.exe module. 
This program has been written in C++ and can in principle be compiled under any OS supporting 
standard C++ compilers and execution environment. However, at present the executable 
module is only available for 64-bit Windows OS. Any Windows x64 PC can in principle be 
used for running this code, but multi-core CPU systems are recommended for faster 
execution.

This program is based on the DHT algorithm described in [T.E. Gureyev, H.M. Quiney, A. 
Kozlov, D.M. Paganin, G. Schmalz and L.J. Allen, Relative roles of multiple scattering and 
Fresnel diffraction in the imaging of small molecules using electrons, Part II: Differential 
Holographic Tomography, arXiv 2012.07012].

The two main "modes" of execution of PhaseRetrieval.exe correspond to:  
(1) reconstruction of a 3D distribution of the electrostatic potential from defocused 
images (or complex amplitudes) generated by MsctKirkland.exe program for a set of 
different orientations of the imaged object and different defocus distances;
(2) reprojection through a 3D distribution of the electrostatic potential that was either 
generated by MsctKirkland.exe from an XYZ file or, more commonly, was previously reconstructed
by the PhaseRetrieval.exe program itself from defocused images. The result of the reprojection
operation consists of a set of defocused complex amplitudes calculated at the rotational orientations
and defocus distances specified in the text file from Parameter 1 for the input data.
As follows from this description, "mode 2" is not mutually exclusive with "mode 1", and the two
modes can be performed sequentially in a single execution of PhaseRetrieval.exe.

Also, the result of the "reprojection" operation ("mode 2") (a set of defocused complex amplitudes)
can in turn be used as new input for PhaseRetrieval.exe, giving the user an option for manual iteration
of the DHT reconstruction. The latter option can in priniciple be used for iterative refinement of
the DHT reconstruction, however this option has not been fully tested yet at this time.

The control of execution of PhaseRetrieval.exe program is managed via an editable text file 
PhaseRetrieval.txt which contains all modifiable input parameters for PhaseRetrieval.exe. In 
order for the program to run, the input parameter file PhaseRetrieval.txt must be present in 
the same folder where PhaseRetrieval.exe is started from. 

The format of PhaseRetrieval.txt file allows any number of "comment" lines to be present
at the beginning and/or at the end of the file, each such line starting with two forward
slash symbols, "//". The comment lines are ignored by by PhaseRetrieval.exe. Other 
(non-comment) lines must all have the same structure and the fixed sequential order with 
respect to each other. Each non-comment line consists of:
(a) an arbitrary "expression" with no white spaces, which typically represents the name of a 
parameter. The contents of this parameter name expression are ignored by 
PhaseRetrieval.exe. An example of the parameter name 
is "9.Output_defocus_distances_min_max_and_step_in_angstroms:".
(b) one or more alphanumerical entries separated by a single white space from each other 
and from the parameter name, each such entry representing one parameter value, for example 
" -30.0 0.0 0.1171875", which corresponds to three parameter values equal to -30.0, 0.0 and 
0.1171875.
Note that the number of white spaces separating the parameters must be equal to one, 
otherwise a "zero" parameter may be wrongly read in. There should be also no further white 
spaces after the last parameter in the line. Each line should be terminated by the usual "end-
of-line / new line / caret-return", etc., symbols, as common for a given OS. The number and 
types of parameters in each parameter line of the file PhaseRetrieval.txt is pre-determined and 
cannot be changed, only the values of these parameters can be edited. The parameter names 
and comment lines are supposed to give explicit hints as to what type and how many 
parameters are expected in each line. 

Here is a typical example of a valid PhaseRetrieval.txt parameter file.

//*****Input parameter file for PhaseRetrieval.exe program*****
1.Input_file_with_rotation_angles_and_defocus_distances: DefocusRand36_NEW.txt
2.Filename_template_of_defocus_series_in_TIFF,_GRD_or_GRC_format: C:\Users\Downloads\Temp\asp.grd
3.Xmin_Xmax_Ymin_Ymax_in_Angstroms_(for_TIFF_input_files_only): 0 10 0 10
4.Wavelength_in_Angstroms: 0.025
5.Objective_aperture_in_mrad: 0.0
6.Spherical_aberration_Cs3,_Cs5_in_mm: 0 0
7.Phase_retrieval_method_IWFR(1),_CTFL2(2),_MinLogAmp(3)_or_EqB7(4): 1
8.Save_phase_retrieved_defocused_complex_amplitudes_in_files_Yes(1)_or_No(0): 0
9.Maximal_number_of_iterations(IWFR): 200
10.Minimal_phase_reconstruction_error(IWFR)_or_regularization_parameter(CTFL2): 1.e-8
11.Output_defocus_distances_min_max_and_step_in_Angstroms: -10.0 0.0 0.0390625
12.Extra_defocus_for_3D_reconstruction_in_Angstroms: 1.0
13.Low_threshold_value_in_Volts_and_multiplicative_scaling_factor_for_3D_potential: 200 0.9
14.3D_Laplacian_filter_mode:_not_apply(0)_apply(1)_apply_filter_only(2): 0
15.Regularization_parameter_for_inverse_3D_Laplacian_filter: 0.01
16.Reprojection_of_3D_potential:_not_apply(0)_apply(1)_apply_to_imported_potential(2): 1
17.Slice_thickness_for_multislice_in_Angstroms: 0.1375
18.File_name_base_of_3D_potential_in_TIFF_or_GRD_format: C:\Users\Downloads\Temp\out.grd
19.Verbose_output_during_execution_Yes(1)_or_No(0): 0
20.Number_of_parallel_threads: 20
//*****
//*****!!!Don't forget to save this file after editing the parameters

Note that the "numbering" of parameters in the parameter names, such as "10" in 
"10.Output_defocus_distances_min_max_and_step_in_angstroms:" is inconsequential and is 
ignored by PhaseRetrieval.exe. On the other hand, the order of the parameter lines in this file 
is very important and must be the same as in the above example. In other words, 
PhaseRetrieval.exe interprets the input parameters according to the order of the non-comment 
lines in the parameter file, while ignoring any numbering that may or may not appear in the 
parameter names.

In the above example, Parameter 1 contains the name of an input file with rotation angles and 
defocus distances, which is usually the same file that was previously used by 
MsctKirkland.exe program to generate the set of input defocused images. The file may contain an
arbitrary number of leading comment lines, each starting with a pair of forward slashes "//". 
Non-comment lines all have the same structure, but may contain different number of entries, 
according to the following scheme. Columns one and two always contain the rotation angles (in degrees)
around the Y and X' axes, respectively (where X' axis is the X axis after the initial rotation
of the 3D space around the Y axis). These two angles are followed by an arbitrary number of pairs
of values, with each odd column (starting from the third one) containing the angle (in degrees)
of rotation around the Z" axis (which corresponds to the Z axis after the initial rotation around
the Y and X' axes). Each even column (starting from the fourth one) contains a defocus distance in 
angstroms along the Z" axis. The number of lines in this file is not limited.
An example of a valid "rotation angles and defocus distances" file is given below. This file must be present
in the same folder where PhaseRetrieval.exe is started from, or, alternatively, the filename can
include a fully specified pathname (OS specific). 

Parameter 2 contains a fully specified pathname template for the input files containing either 
defocused images (in GRD format or uncompressed 32-bit floating-point TIFF format) or defocused
complex amplitudes (in GRC format). Note that this template is used for generation of the eventual
input and output filenames which PhaseRetrieval.exe creates by inserting a single number or a pair of
numbers that correspond to the index of a defocus distance (first inserted number) and the index of a 
rotational orientation (second inserted number). As a result, when a filename template 
"ag.grd" is given in this parameter, the output filenames may look like "ag0_00.grd, 
ag0_01.grd, ..., ag0_35.grd, ag1_00.grd, ..., ag1_35.grd, ag2_00.grd, ..., ag2_35.grd" (36 x 3 
= 108 files in total), if 36 different rotational positions, with 3 defocus distances at each 
rotational position, are specified in the file whose name is given in 
Parameter 1. See the description of GRD and GRC file formats and an example 
DefocusRand36_NEW.txt file below. Note that when input defocused complex amplitudes 
are specified, there can be only one GRC file per illumination direction, and the subsequent 
2D phase retrieval (see below) is skipped. The reason for such a behaviour of the program is that
each complex amplitude is assumed to contain the correct complex wavefunction at the given orientaion
and the defocused distance. In principle, this complex amplitude can be propagated
to any other defocused distance by calculating Fresnel diffraction integrals. In particular,
no phase retrieval is required in this case. Obviously, such GRC (complex amplitude) input is intended
for testing and debugging purposes only, since only defocused intensities are expected to be available
under experimentally-relevant conditions.

Parameter 3 contains the physical dimensions (in angstroms) of the images specified in Parameter 2 in
the case when these images are given in TIFF files. This parameter is ignored when the input images
are given in GRD or GRC format, since these files contain the physical dimensions of the images internally.

Parameter 4 contains the wavelength of the illuminating wave in angstroms. 

Parameter 5 contains the objective aperture in milliradians. Zero or negative value here is 
interpreted as an infinite aperture.

Parameter 6 contains the values of the spherical aberrations Cs3 and Cs5 (in millimetres) of 
the imaging system. Zero values correspond to the absence of aberrations.

Parameter 7 can be equal to 1, 2, 3 or 4. It determines the 2D phase retrieval method to be used as 
part of the DHT algorithm internally in the PhaseRetrieval.exe program. The first two methods
can be used only in the case when there 2 or more defocused images per each illumination directon,
while the second two methods can be used only when there is 1 defocused image at each illumination 
direction. "1" corresponds to the IWFR phase retrieval method [L. J. Allen, W. McBride, N. L. O'Leary
and  M. P. Oxley, Ultramicroscopy 100 (2004) 91-104.]. "2" corresponds to the CTF-L2 phase retrieval 
method [D. Paganin, A. Barty, P. J. Mcmahon and K. A. Nugent, Quantitative phase-amplitude 
microscopy. III. The effects of noise, J. Micros. (2004) 51-61.] "3" corresponds to the phase retrieval
from a single intensity according to phase = -0.5 * log(Intensity). "4" corresponds to the phase
retrieval from a single intensity distribution according to phase = -0.5 * sqrt(Kmax^2 - K^2), where
K is the contrast function, K = 1 - Intensity, and Kmax is the maximum value of the contrast function
for this defocused image. When input defocused complex amplitudes are specified in Parameter 2, the 
2D phase retrieval step is skipped.

Parameter 8 specifies if the phase-retrieved defocused complex amplitudes are to be saved in 
GRC files. If this parameter is equal to 1, then the phase-retrieved defocused complex amplitudes
are saved to GRC files, with the names that are created by modifying the input defocused intensity
filenames by changing the file extensions to ".grc" and adding the capital letter "D" before that
file extension.

Parameter 9 specifies the maximal allowed number of iterations to be used in the IWFR 
algorithm. This parameter is ignored in the cases when there IWFR method is not used.

Parameter 9 specifies, depending on the value of Parameter 7, either the minimal phase 
reconstruction error (typical value ~ 1.e-8) which leads to the interruption of the IWFR iterations
(note that in the "verbose" mode - see below - PhaseRetrieval.exe outputs the value of the phase 
reconstruction error at each IWFR iteration, so these values can be used as a guide for setting a 
sensible value of the minimal error in this parameter), or the Tikhonov regularization 
parameter (typical value ~ 1.e-12) for the CTF-L2 algorithm (setting this parameter to zero forces the 
internal implementation of the CTF-L2 algorithm to use an estimated "optimal" value of this 
parameter which can be different at each illumination direction). When the input 
dataset contains only one defocused intensity image for a given illumination direction (orientation), 
the value of this parameter is ignored.

Parameter 10 contains three related input parameters: the minimum, maximum and the step (z-
distance between successive planes) in angstroms of the output (xy) planes at which the 3D 
potential is reconstructed. The (xy) extent and the grid steps in this volume are taken from the 
input defocused images (the GRD and GRC file formats contain the relevant information - 
see the description of the GRD and GRC formats below). Note that this program requires the 
x, y and z grid steps to be the same, i.e. it always uses cubic voxels in the reconstructed 
volume.

Parameter 11 contains the "extra defocus" distance in angstroms used in the DHT algorithm 
for 3D reconstruction. Typical values of this parameter can be between 1 and 10 angstroms. 
See details in [T.E. Gureyev et al., arXiv 2012.07012].

Parameter 12 contains two "technical" parameters used in the 3D reconstruction. These are 
the low threshold value in volts and the multiplicative scaling factor (dimensionless) for the 
reconstructed 3D potential. Any values of the reconstructed potential that are lower than this 
threshold value are replaced by zeros (a typical potential is assumed to be mostly positive). 
This operation helps to get rid of the erroneous background that may appear due to the 
limited number of available views (illumination directions) and also due to multiple 
scattering. The reconstructed potential is also multiplied by the "scaling factor", which may 
typically be between 0.1 and 1. This multiplication happens prior to the thresholding. The optimal
value of this scaling factor should be determined by the user on the basis of comparison of the
output 3D potential with a priori known realistic physical values for the atomic potentials. 
Another method for fine-tuning of this scaling factor can be based on the comparison of the contrast
in the re-projected defocused images with the contrast in the input defocused images. Larger
values of the scaling parameter result in larger 3D potential values, which in turn leads to 
larger contrast values in the defocused images (and vice versa). Note that the contrast values
for the input and re-projected defocused images are printed out by PhaseRetrieval.exe before the
end of the execution.

Parameter 13 can be equal to 0, 1 or 2. It determines if and how the inverse 3D Laplacian 
filter is applied to the 3D potential. "0" corresponds to not applying the inverse 3D Laplacian. 
"1" corresponds to applying it. "2" corresponds to skipping the initial 3D potential 
reconstruction from defocused images, reading the 3D potential distribution from the files 
defined in Parameter 16 below and applying the inverse 3D Laplacian filter to the potential 
imported from the files. A typical usage cycle of PhaseRetrieval.exe may consist of an initial 
reconstruction of the 3D potential from defocused images and saving them into files in the 
first run of PhaseRetrieval.exe, with Parameter 13 set to zero. This may be followed by a 
second run (restart) of PhaseRetrieval.exe program, with a suitably modified parameter file, 
where Parameter 13 is set to 2, in particular. In the latter case, PhaseRetrieval.exe saves the 
filtered 3D potential in the GRD or TIFF files with the names that are obtained from the names 
defined by Parameter 17, but with inserted letter "L" at the end of the filename, immediately 
preceding the file extension. 

Parameter 14 contains the regularization parameter for the inverse 3D Laplacian filter. If the 
filter is not used, this parameter is ignored.

Parameter 15 can be equal to 0, 1 or 2. When this parameter is set to 0, the multislice 
reprojection of the reconstructed 3D potential is not applied. When this parameter is set to 1, 
the multislice reprojection of the reconstructed 3D potential is applied. When this parameter 
is set to 2, the initial reconstruction of the 3D potential is skipped, the distribution of the 
potential is read from the files defined by Parameter 17 and the multislice reprojection is 
applied to the imported potential. When Parameter 15 is equal to 1 or 2, the defocused 
complex amplitudes are calculated at the angular positions defined in the file specified in 
Parameter 1. These defocused complex amplitudes are saved in the GRC files with the names 
that correspond to those in Parameter 2, but with the file extension always replaced by ".grc" 
and the letter "R" inserted at the end of the filename, immediately preceding the file 
extension. Note that the intensity of these reprojected defocused complex amplitudes may or 
may not be replaced by the intensity read from the files given by Parameter 2 (this feature has 
not been fully settled yet). The phase of the reprojected defocused complex amplitudes 
always corresponds to the phase that was determined by the multislice projection through the 
reconstructed or imported 3D potential and the subsequent free-space propagation to a given 
defocus plane.

Parameter 16 defines the slice thickness in angstroms for the multislice algorithm used for the 
reprojection of the reconstructed or imported 3D potential.

Parameter 17 contains a fully specified pathname for the output or input files (in GRD 
format  or uncompressed 32-bit floating-point TIFF format) containing 2D sections of 
the 3D electrostatic potential where the potential is either 
reconstructed by this program or was produced previously and saved in these files (the 
behaviour depends on the values of Parameters 13 and 15). Note that here the pathname 
contains a "template" for the eventual output filenames which PhaseRetrieval.exe creates by 
inserting numbers that correspond to the z-index of a (xy) plane with a 2D section of the 3D 
potential. As a result, when a filename template "out.grd" is given in this parameter, the 
actual filenames with the 3D potential may look like "out000.grd, out001.grd,..., out255.grd", 
if 256 (xy) planes are specified as a consequence of the values in Parameter 9.

Parameter 18 determines the amount of diagnostic information that PhaseRetrieval.exe 
outputs during the execution. 

Parameter 19 contains the desired number of worker threads to launch during the execution of 
PhaseRetrieval.exe. The recommend number is equal to the hyper-threading capacity of one's 
computer (often, this number is equal to the number of available CPU cores times 2), minus 1 
(one thread may be reserved for the "main user thread" that leaves the computer responsive to 
user interactions during prolonged calculations). For example, if you run PhaseRetrieval.exe 
on a PC with an Intel Core i5 CPU with 6 cores, the recommended value of this parameter is 11.
 
============================================================================================
** Example of a text file whose name may appear in Parameter 1 (this is typically the same 
file as used in MsctKirkland.exe for the generation of the defocused images that are used here 
as the input dataset):

//Y rotation	X' rotation	Z" rotation1	Defocus1	Z" rotation2"	Defocus2	Z" 
rotation3	Defocus3
//theta=360.0*rand()	fi=180/PI()*ACOS(2*RAND()-1)-90		psi=360.0*rand()
					
26.00810694	-8.522712116	129.6892203	6.743773574	288.416648	14.6618478
	166.2917111	9.563806875
172.5031608	33.49812809	299.8940784	10.41893232	297.2303384	8.923491161
	77.87911711	8.478823205
190.4084923	35.29255135	149.6711892	8.466044745	317.3861666	11.12645956
	282.8326037	9.782168557
281.8396809	31.86876039	247.383095	6.684336937	282.8837994	10.91649724
	129.0044657	9.249695759
121.8589718	-15.73317918	166.9345344	9.128775043	197.1203744	9.460768742
	104.9953255	7.066174276
13.23660864	-35.52467996	180.3408877	11.1417195	332.2446587	9.793296018
	254.2326199	6.010249766
97.94685887	33.81912113	256.5517521	14.53473129	77.00788592	5.398747561
	124.534224	8.885638777
287.9758582	-39.55341516	120.7328669	10.44514973	297.0706126	13.33282855
	196.1524186	6.334542628
168.1278053	-6.837024217	243.30537	8.865746818	177.9092585	7.48798294
	341.794898	5.270372678
261.4319794	-41.61580786	347.6080615	9.856249747	57.23544936	8.711746802
	44.70407416	6.903986789
185.5423696	-4.107074411	105.0992619	12.23325712	12.71098296	10.50607993
	15.66657027	7.327044709
192.7022474	82.89589771	162.0607763	14.5610955	218.6979113	8.718024357
	264.5740903	6.05555455
20.17007635	42.58198525	93.13592458	7.840112874	165.1363809	13.92938044
	157.8398247	10.68574609
333.2362849	77.1940075	246.2594822	13.44612092	347.7083846	5.461197781
	46.82300175	13.84238007
303.113637	-21.34731924	180.9085991	6.730433274	67.02125086	13.13993847
	187.0145653	11.43568727
10.56068821	-53.2860474	31.31217606	8.425675115	92.43879993	7.163533536
	87.90094468	12.1309391
186.9984386	54.28811639	12.14141647	14.16775017	59.77924412	12.19133307
	200.6417596	8.112978462
110.2506822	27.40563189	216.0672076	12.38859577	325.9786977	9.611240633
	225.0603636	9.003978108
184.9159323	-3.625488046	135.932916	6.802278071	212.5119307	13.14448392
	281.4830319	6.647482517
77.89779902	-13.56512848	113.0782808	13.97810233	204.9260522	11.93188928
	325.8090956	12.27359692
300.9113031	54.43709417	116.1180724	5.379896969	275.4608852	7.65517426
	88.88008278	9.000631351
154.8188017	-43.3613752	179.4017215	8.567061926	134.5644319	13.26329474
	207.5203022	5.367051864
259.9839366	-19.1904972	322.2485636	11.84008175	131.3935473	14.8144175
	310.6652163	10.97259488
101.5203286	-41.29122375	115.9810789	13.94433717	169.6087817	5.449410631
	350.9521245	13.40051361
94.33638037	12.49272118	147.4698991	8.735070032	254.3671792	10.97841944
	133.5247864	9.517436808
299.5040164	-6.041900305	79.58697258	10.81277137	69.33081869	14.59612232
	316.4086498	12.79662275
88.96933464	-1.45922459	315.1460236	6.104582973	351.6848002	12.30385844
	326.8911979	11.15975927
247.509192	43.01790127	288.8351795	9.368473515	285.0750141	6.44152244
	6.20089253	12.52895663
135.4598683	32.12838575	328.9430274	9.001046462	277.5515331	8.743884584
	164.0206539	5.660730802
5.416267431	-30.52809612	243.9488404	10.70592495	29.84028504	8.439657455
	335.1056686	8.661439169
247.3367065	21.0140705	125.6417566	9.213996425	98.33586814	7.68789426
	173.1604265	11.58097361
250.1751286	-7.94078762	314.7840398	12.15535647	264.4284848	9.84649179
	218.397115	5.142187996
130.1482017	54.30189826	148.5037788	13.91551908	91.78206487	14.57023721
	32.04432909	7.856108786
222.4045395	45.62699219	349.3502007	5.718044183	8.507734407	14.31396298
	235.7349004	11.52369882
18.22662693	1.242340419	3.991048553	13.13944878	350.0059348	8.593798703
	320.1890816	11.54488428
230.2144586	24.8655988	209.3497261	7.902525523	50.12716802	12.11207933
	28.4281332	8.219013905

=============================================================================== 
Specifications of the GRD file format:

GRD files contain:
(1) string "DSAA" or "DSBB" (in ASC and BIN GRD files, respectively)
(2) nx ny - array dimensions as "short"s (in ASC format, separated by a white space)
(3) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a 
white space)
(4) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a 
white space)
(5) zlo zhi - minimum and maximum array value as "double"s (in ASC format, separated by a 
white space)
(6) u[i][j] - array values as "float"s (in ASC format, separated by a white space), j index 
changes most rapidly
 
=================================================================================
Specifications of the GRC file format

GRC files contain:
(1) -5 - as a "short" value (this is the GRC file format identifier)
(2) wl - wavelength as a "double"
(3) nx ny - array dimensions as "short"s (in ASC format, separated by a white space)
(4) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a 
white space)
(5) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a 
white space)
(6) real(u[i][j]) imag(u[i][j]) - real and imaginary parts of array values as "float"s (in ASC 
format, separated by a white space), j index changes most rapidly.