#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "pdb.h"
#include "Seb.h"

// This program reads a PDB, Vesta XYZ or Kirkland XYZ file, 
// optionally centers the position of the "molecule" within a given "CT cube" [0, ctblength] x [0, ctblength] x [0, ctblength],
// optionally rotates the "molecule" by a given angle around the y axis, 
// optionally removes "duplicate" atomes, and
// outputs the data in the form of a muSTEM, Vesta XYZ or Kirkland XYZ file, while it also
// sorts all atoms in descending order with respect to atom weights, and sorts all atoms with the same weight in ascending order with respect to z coordinate
//
// This programs uses the "miniball" code written by Martin Kutz <kutz@math.fu-berlin.de> and Kaspar Fischer <kf@iaeth.ch> for finding the minimal ball enclosing the molecule
// See https://github.com/hbf/miniball for authors, theory and copyright information


int main(void)
{
	int i, j, k;
	char pdbfile[1024]; // input file name
	char outfile[1024]; // output file name
	char cfileinfo[1024]; // free-form 1st line in the output file

	char cline[1024], ctitle[1024], cparam[1024]; // auxiliary storage


	// Read input parameter file pdb.txt
	printf("\nStarting ...");
	FILE* ffpar = fopen("pdb.txt", "rt");
	if (!ffpar)
	{
		printf("\nInput parameter file pdb.txt not found!!!\n");
		return -1;
	}

	// line 0 in input parameter file
	fgets(cline, 1024, ffpar); // 0 line - comment

	// line 1
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 1nd line: input file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile) != 2)
	{
		printf("\n!!!Error reading input file name from input parameter file.");
		return -1;
	}

	// line 2
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 2nd line: input file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input file type from input parameter file.");
		return -1;
	}
	int nfiletype = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype = atoi(cparam);
	if (nfiletype != 0 && nfiletype != 1 && nfiletype != 2)
	{
		printf("\n!!!Unknown input file type in input parameter file.");
		return -1;
	}

	// line 3
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //3rd line: make all coordinates positive and centre them, or not 
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading 'make all coordinates positive' switch from input parameter file.");
		return -1;
	}
	int npositive = 0; // 1 - make all output coordinates positive and centre the sample inside the CT sample cube, 0 - do not do this
	npositive = atoi(cparam);
	if (npositive != 0 && npositive != 1)
	{
		printf("\n!!!Unknown value of the 'make all output coordinates positive' switch in input parameter file.");
		return -1;
	}

	// line 4
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //4th line: CT cube side length
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading CT cube side length from input parameter file.");
		return -1;
	}
	double ctblength = 0; // CT box length
	ctblength = atof(cparam);
	if (npositive == 1)
	{
		if (ctblength <= 0)
		{
			printf("\n!!!CT cube side length %g is not positive in pdb.txt!!!", ctblength);
			return -1;
		}
	}
	else
		printf("\n!!!CT cube side length parameter from the input file will be ignored!!!");

	// line 5
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //5th line: rotation angle around y
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading rotation angle around y axis from input parameter file.");
		return -1;
	}
	double angley = 0; // rotation angle 
	angley = atof(cparam) * 3.141592653589793 / 180.0;
	double cosangley = cos(angley);
	double sinangley = sin(angley);

	// line 6
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //6th line: rotation angle around x' (i.e. around x axis after the initial rotation around y axis)
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading rotation angle around x axis from input parameter file.");
		return -1;
	}
	double anglex = 0; // rotation angle 
	anglex = atof(cparam) * 3.141592653589793 / 180.0;
	double cosanglex = cos(anglex);
	double sinanglex = sin(anglex);

	// ... line 
	//fgets(cline, 1024, ffpar); strtok(cline, "\n");  //...th line: rotation angle around z // this rotation is redundant
	//if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	//{
	//	printf("\n!!!Error reading rotation angle around z axis from input parameter file.");
	//	return -1;
	//}
	//double anglez = 0; // rotation angle 
	//anglez = atof(cparam) * 3.141592653589793 / 180.0;
	//double cosanglez = cos(anglez);
	//double sinanglez = sin(anglez);

	// line 7
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //7th line: maximum distance to remove duplicates
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading maximum distance to remove duplicates from input parameter file.");
		return -1;
	}
	double dupdist = 0.0;
	dupdist = atof(cparam);
	if (dupdist < 0) dupdist = 0.0;
	double dupdist2 = dupdist * dupdist;

	// line 8
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 8th line: output data sort
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading output data sort type from input parameter file.");
		return -1;
	}
	int noutsort = 0; // output file sort 0 - no sort, 1 - sort by ascending order of z coordinate, 2- sort by descending order of the occupancy values
	noutsort = atoi(cparam);
	if (noutsort != 0 && noutsort != 1 && noutsort != 2)
	{
		printf("\n!!!Unknown output data sort type in the input parameter file.");
		return -1;
	}

	// line 9
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 9th line: output file name
	if (sscanf(cline, "%s %s", ctitle, outfile) != 2)
	{
		printf("\n!!!Error reading output file name from input parameter file.");
		return -1;
	}

	// line 10
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 10th line: output file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading output file type from input parameter file.");
		return -1;
	}
	int noutfiletype = 0; // output file type 0 -  for Kirkland-format XYZ file, 1 - for muSTEM-foram XTL file
	noutfiletype = atoi(cparam);
	if (noutfiletype != 0 && noutfiletype != 1 && noutfiletype != 2)
	{
		printf("\n!!!Unknown output file type in the input parameter file.");
		return -1;
	}

	// line 11
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 11th line: free form first line for the output file
	if (sscanf(cline, "%s %s", ctitle, cfileinfo) != 2)
	{
		printf("\n!!!Error reading free-form line for the output file from input parameter file.");
		return -1;
	}

	fclose(ffpar); // close input parameter file


	// read input PDB, Vesta or Kirkland XYZ file
	pdbdata pd;
	pdbdata_init(&pd);
	if (nfiletype == 0)
	{
		printf("\nReading input file %s in PDB format ...", pdbfile);
		if (read_pdb(pdbfile, &pd) == -1) return -1; // read PDB file
	}
	else if (nfiletype == 1)
	{
		printf("\nReading input file %s in Vesta export XYZ format ...", pdbfile);
		if (data_from_VestaXYZfile(pdbfile, &pd) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype == 2)
	{
		printf("\nReading input file2 %s in Kirkland XYZ format ...\n", pdbfile);
		if (data_from_KirklandXYZfile(pdbfile, &pd) == -1) return -1; // read Kirkland XYZ file
	}

	//print_pdb_data(&pd);

	// find the spatial bounds of the molecule
	double xmin, xmax, ymin, ymax, zmin, zmax, xc0 = 0, yc0 = 0, zc0 = 0;
	xmin = xmax = pd.adata[0].x;
	ymin = ymax = pd.adata[0].y;
	zmin = zmax = pd.adata[0].z;

	for (i = 1; i < pd.natoms; i++)
	{
		if (pd.adata[i].x < xmin) xmin = pd.adata[i].x;
		else if (pd.adata[i].x > xmax) xmax = pd.adata[i].x;
		if (pd.adata[i].y < ymin) ymin = pd.adata[i].y;
		else if (pd.adata[i].y > ymax) ymax = pd.adata[i].y;
		if (pd.adata[i].z < zmin) zmin = pd.adata[i].z;
		else if (pd.adata[i].z > zmax) zmax = pd.adata[i].z;
		xc0 += pd.adata[i].x;
		yc0 += pd.adata[i].y;
		zc0 += pd.adata[i].z;
	}
	xc0 /= pd.natoms;
	yc0 /= pd.natoms;
	zc0 /= pd.natoms;
	double xsize = xmax - xmin;
	double ysize = ymax - ymin;
	double zsize = zmax - zmin;
	double diam = floor(sqrt(xsize * xsize + ysize * ysize + zsize * zsize)) + 1.0; // we want to round the diameter up and ensure that atoms won't stick out after rotation due to numerical precision

	printf("\nThe x, y and z extent of the molecule in the input file was:");
	printf("\n  [%g, %g], [%g, %g], and [%g, %g], respectively.", xmin, xmax, ymin, ymax, zmin, zmax);
	printf("\nThe diagonal of the minimal parallelpiped enclosing the molecule is %g.", diam);

	if (npositive == 1) // if the shift and centring of the output atomic coordinates is required
	{
		double r2;
		double r2max = (pd.adata[0].x - xc0) * (pd.adata[0].x - xc0) + (pd.adata[0].y - yc0) * (pd.adata[0].y - yc0) + (pd.adata[0].z - zc0) * (pd.adata[0].z - zc0);
		for (i = 1; i < pd.natoms; i++)
		{
			r2 = (pd.adata[i].x - xc0) * (pd.adata[i].x - xc0) + (pd.adata[i].y - yc0) * (pd.adata[i].y - yc0) + (pd.adata[i].z - zc0) * (pd.adata[i].z - zc0);
			if (r2 > r2max) r2max = r2;
		}
		diam = floor(2.0* sqrt(r2max)) + 1.0; // we want to round the diameter up and ensure that atoms won't stick out after rotation due to numerical precision
		printf("\nThe average-coordinate centre of the molecule was located at the point x = %g, y = %g, z = %g.", xc0, yc0, zc0);
		printf("\nThe minimal diameter of the ball with this centre and enclosing the molecule is %g.", diam);

		// Miniball code
		// See https://github.com/hbf/miniball for authors, theory and copyright information
		typedef double FT;
		typedef Seb::Point<FT> Point;
		typedef std::vector<Point> PointVector;
		typedef Seb::Smallest_enclosing_ball<FT> Miniball;
		PointVector S;
		std::vector<double> coords(3);
		for (int i = 0; i < pd.natoms; ++i) 
		{
			coords[0] = pd.adata[i].x;
			coords[1] = pd.adata[i].y;
			coords[2] = pd.adata[i].z;
			S.push_back(Point(3, coords.begin()));
		}
		Miniball mb(3, S);
		diam = floor(2.0 * mb.radius()) + 1.0; // we want to round the diameter up and ensure that atoms won't stick out after rotation due to numerical precision
		Miniball::Coordinate_iterator center_it = mb.center_begin();
		xc0 = center_it[0];
		yc0 = center_it[1];
		zc0 = center_it[2];
#ifdef _DEBUG
		mb.verify();
#endif
		printf("\nThe centre of the minimal ball enclosing the molecule was located at the point x = %g, y = %g, z = %g.", xc0, yc0, zc0);
		printf("\nThe diameter of this minimal ball is %g.", diam);

		if (ctblength < diam)
		{
			ctblength = diam;
			printf("\n!!! CT_sample_cube_side_length parameter in the input parameter file was too small and has been replaced with %g.", diam);
		}

		xsize = ysize = zsize = ctblength;

		double xc = 0.5 * ctblength;
		double yc = 0.5 * ctblength;
		double zc = 0.5 * ctblength;

		printf("\nThe molecule will be shifted to the cube with side length %g and positive coordinates.", diam);
		printf("\nThe new centre of the molecule will be located at the point x = %g, y = %g, z = %g.\n", xc, yc, zc);

		double xshift = xc - xc0;
		double yshift = yc - yc0;
		double zshift = zc - zc0;

		// TEG: we make all coordinates non-negative, as it seems that Kirkland's software expects this,
		// and we also centre each of XYZ coordinates within the interval [0, ctblength], so that the sample can remain within the cube during rotation
		// and we rotate the molecule around the y axis by the angley and around the z axis by the anglez set in the input parameter file
		double xk, yk, zk;
		for (i = 0; i < pd.natoms; i++)
		{
			// centering
			xk = pd.adata[i].x + xshift;
			yk = pd.adata[i].y + yshift;
			zk = pd.adata[i].z + zshift;

			//rotation around y axis
			pd.adata[i].x = xc + (xk - xc) * cosangley + (zk - zc) * sinangley;
			pd.adata[i].y = yk;
			pd.adata[i].z = zc + (-xk + xc) * sinangley + (zk - zc) * cosangley;

			//rotation around x axis
			yk = pd.adata[i].y;
			zk = pd.adata[i].z;
			pd.adata[i].y = yc + (yk - yc) * cosanglex + (zk - zc) * sinanglex;
			pd.adata[i].z = zc + (-yk + yc) * sinanglex + (zk - zc) * cosanglex;

			//rotation around z axis
			//xk = pd.adata[i].x;
			//yk = pd.adata[i].y;
			//pd.adata[i].x = xc + (xk - xc) * cosanglez + (yk - yc) * sinanglez;
			//pd.adata[i].y = yc + (-xk + xc) * sinanglez + (yk - yc) * cosanglez;
		}
	}

	//translate element symbols into atomic weights
	int* ia = (int*) malloc(pd.natoms * sizeof(int));
	if (nfiletype == 2)
	{
		for (i = 0; i < pd.natoms; i++) ia[i] = pd.adata[i].serial;
		pdb_symbols(&pd, ia);
	}
	else if (pdb_atomnumbers(&pd, ia))
	{
		printf("\n!!!Error encountered while finding atomic weight for a given element name!!!\n");
		return -1; // the function pdb_atomnumbers will print its error messages before exiting
	}

	// sort entries by atom weight in descending order
	pdb_bubbleSort1(&pd, ia);

	// find the number of different atom types (this assumes that all entries have been sorted in descending order)
	int natypes = 1;
	for (i = 0; i < pd.natoms - 1; i++)
		if (ia[i + 1] != ia[i]) natypes++;

	// find the number of atoms for each type
	int* vtypes = (int*) malloc(natypes * sizeof(int));
	j = 0;
	vtypes[0] = 1;
	for (i = 0; i < pd.natoms - 1; i++)
	{
		if (ia[i + 1] != ia[i]) { j++; vtypes[j] = 1; }
		else vtypes[j]++; // count in the next(!) atom
	}

	// print out the numbers of atoms for each atom type present in the input file
	printf("\nInput file contains %d different types of atoms", natypes);
	k = 0;
	for (j = 0; j < natypes; j++)
		for (i = 0; i < vtypes[j]; i++)
		{
			if (i == 0) printf("\nThere are %d %s atoms (with atomic weight = %d) in the input file.", vtypes[j], pd.adata[k].element, ia[k]);
			k++;
		}
	printf("\n");

	// sort atoms of each type by z coordinate in ascending order or occupancy values in descending order
	if (noutsort == 1)
	{
		k = 0;
		for (j = 0; j < natypes; j++)
		{
			if (pdb_bubbleSort2(&pd, k, k + vtypes[j]) != 0) return -1;
			k += vtypes[j];
		}
	}
	else if (noutsort == 2)
	{
		k = 0;
		for (j = 0; j < natypes; j++)
		{
			if (pdb_bubbleSort3(&pd, k, k + vtypes[j]) != 0) return -1;
			k += vtypes[j];
		}
	}

    // find and remove duplicates
	double r2;
	for (i = 0; i < pd.natoms; i++)
	{
		if (pd.adata[i].tempFactor == -777) continue; // if this atom has already been marked as duplicate, skip it
		for (j = i + 1; j < pd.natoms; j++)
		{
			if (pd.adata[j].tempFactor == -777) continue; // don't compare atom with those already classified as duplicates
			r2 = (pd.adata[i].x - pd.adata[j].x) * (pd.adata[i].x - pd.adata[j].x) + (pd.adata[i].y - pd.adata[j].y) * (pd.adata[i].y - pd.adata[j].y) + (pd.adata[i].z - pd.adata[j].z) * (pd.adata[i].z - pd.adata[j].z);
			if (r2 < dupdist2) // atom j is too close to atom i
			{
				if (pd.adata[i].occupancy < pd.adata[j].occupancy) 
				{
					pd.adata[i].tempFactor = -777; 
					printf("\n!!! Found duplicate atoms entries no. %d and %d (0 based), distance = %g Angstroms.", i, j, sqrt(dupdist2));
					printf("\n!!!   Atom no. %d is removed", i);
					break;
				}  // mark this atom as duplicate
				else
				{
					pd.adata[j].tempFactor = -777; // mark this atom as duplicate
					printf("\n!!! Found duplicate atoms entries no. %d and %d (0 based), distance = %g Angstroms.", i, j, sqrt(dupdist2));
					printf("\n!!!   Atom no. %d is removed", j);
				}
			}
		}
	}

	// count the number of non-duplicate atoms
	int natoms1 = 0;
	for (i = 0; i < pd.natoms; i++) if (pd.adata[i].tempFactor != -777) natoms1++;
	printf("\n In total, %d duplicate atoms found in the input file.", pd.natoms - natoms1);

	// create structure containing the output molecule with removed duplicates
	pdbdata pd1;
	pdbdata_init(&pd1);
	pd1.natoms = natoms1;
	pdballocate(&pd1);
	int* ia1 = (int*) malloc(pd1.natoms * sizeof(int));
	j = 0;
	for (i = 0; i < pd.natoms; i++) 
		if (pd.adata[i].tempFactor != -777)
		{
			ia1[j] = ia[i];
			pd1.adata[j] = pd.adata[i];
			j++;
		}

	// find the number of different atom types (this assumes that all entries have been sorted in descending order)
	int natypes1 = 1;
	for (i = 0; i < pd1.natoms - 1; i++)
		if (ia1[i + 1] != ia1[i]) natypes1++;

	// find the number of atoms for each type
	int* vtypes1 = (int*) malloc(natypes1 * sizeof(int));
	j = 0;
	vtypes1[0] = 1;
	for (i = 0; i < pd1.natoms - 1; i++)
	{
		if (ia1[i + 1] != ia1[i]) { j++; vtypes1[j] = 1; }
		else vtypes1[j]++; // count in the next(!) atom
	}

	// print out the numbers of atoms for each atom type present in the output file
	k = 0;
	for (j = 0; j < natypes1; j++)
		for (i = 0; i < vtypes1[j]; i++)
		{
			if (i == 0) printf("\nThere are %d %s atoms (with atomic weight = %d) in the output file.", vtypes1[j], pd1.adata[k].element, ia1[k]);
			k++;
		}
	printf("\n");


	// output to the target file
	FILE* ff = fopen(outfile, "wt");
	if (noutfiletype == 0) // muSTEM format output
	{
		printf("\nWriting output file %s in muSTEM XTL format ...", outfile);

		fprintf(ff, "%s\n", cfileinfo); // free-form file info line
		fprintf(ff, "%f %f %f 90.0 90.0 90.0\n", xsize, ysize, zsize); // if npositive == 1, ctblength is written as the unit cell dimension
		fprintf(ff, "%f \n", 200.0); // probe energy
		fprintf(ff, "%d \n", natypes1); // number of atoms types

		k = 0;
		for (j = 0; j < natypes1; j++)
		{
			fprintf(ff, "%s\n", pd1.adata[k].element); // next atom type info
			fprintf(ff, "%d %f 1.0 0.0\n", vtypes1[j], (double)ia[k]); // number of atoms of this type, atom weight for this type
			for (i = 0; i < vtypes1[j]; i++)
			{
				fprintf(ff, "%f %f %f\n", pd1.adata[k].x / ctblength, pd1.adata[k].y / ctblength, pd1.adata[k].z / ctblength); // fractional coordinates of an atom
				k++;
			}
		}
		
		free(vtypes1);
	}
	else if (noutfiletype == 1) // Vesta XYZ output
	{
		printf("\nWriting output file %s in Vesta XYZ format ...", outfile);
		fprintf(ff, "%d\n", pd1.natoms);
		fprintf(ff, "%s\n", cfileinfo); // free-form file info line
		for (i = 0; i < pd1.natoms; i++)
		{
			if (pd1.adata[i].occupancy <= 0 || pd1.adata[i].occupancy > 1) pd1.adata[i].occupancy = 1.0;
			fprintf(ff, "%s %f %f %f %f\n", pd1.adata[i].element, pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z, pd1.adata[i].occupancy);
		}
	}
	else if (noutfiletype == 2) // Kirkland XYZ output
	{
		printf("\nWriting output file %s in Kirkland's XYZ format ...", outfile);
		fprintf(ff, "%s\n", cfileinfo); // free-form file info line
		fprintf(ff, "%f %f %f\n", xsize, ysize, zsize); // if npositive == 1, ctblength is written as the unit cell dimension
		for (i = 0; i < pd1.natoms; i++)
		{
			if (pd1.adata[i].occupancy <= 1.e-7 || pd1.adata[i].occupancy > 1) pd1.adata[i].occupancy = 1.0;
			fprintf(ff, "%d %f %f %f %f\n", ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z, pd1.adata[i].occupancy);
		}
		fprintf(ff, "%i\n\n\n", -1);
	}

	fclose(ff);
	free(ia);

	printf("\nFinished!\n");

	return 0;
}