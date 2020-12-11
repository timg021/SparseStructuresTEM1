#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "pdb.h"


int main(void)
{
	char pdbfile[1024], pdbfile1[1024];
	char outfile[1024];
	char cfileinfo[1024];
	char cline[1024], ctitle[1024], cparam[1024]; // auxiliary storage

	// Read input parameter file pdb.txt
	printf("\nStarting ...\n");
	FILE* ffpar = fopen("pdb-compare.txt", "rt");
	if (!ffpar)
	{
		printf("\nInput parameter file pdb-compare.txt not found!!!\n");
		return -1;
	}
	
	// line 0 in input parameter file
	fgets(cline, 1024, ffpar); // 0 line - comment

	// line 1
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 1nd line: input test file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile) != 2)
	{
		printf("\n!!!Error reading input test file name from input parameter file.");
		return -1;
	}

	// line 2
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 2nd line: input test file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input test file type from input parameter file.");
		return -1;
	}
	int nfiletype = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype = atoi(cparam);
	if (nfiletype != 0 && nfiletype != 1 && nfiletype != 2)
	{
		printf("\n!!!Unknown input test file type in input parameter file.");
		return -1;
	}

	// line 3
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 3nd line: input reference file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile1) != 2)
	{
		printf("\n!!!Error reading input reference file name from input parameter file.");
		return -1;
	}

	// line 4
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 4nd line: input reference file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input reference file type from input parameter file.");
		return -1;
	}
	int nfiletype1 = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype1 = atoi(cparam);
	if (nfiletype1 != 0 && nfiletype1 != 1 && nfiletype1 != 2)
	{
		printf("\n!!!Unknown input reference file type in input parameter file.");
		return -1;
	}

	// line 5
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 5th line: output file name
	if (sscanf(cline, "%s %s", ctitle, outfile) != 2)
	{
		printf("\n!!!Error reading output file name from input parameter file.");
		return -1;
	}

	// line 8
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 8th line: free form first line for the output file
	if (sscanf(cline, "%s %s", ctitle, cfileinfo) != 2)
	{
		printf("\n!!!Error reading free-form line for the output file from input parameter file.");
		return -1;
	}

	fclose(ffpar);


	// read input file1
	pdbdata pd;
	pdbdata_init(&pd);
	if (nfiletype == 0)
	{
		printf("Reading input file1 %s in PDB format ...\n", pdbfile);
		if (read_pdb(pdbfile, &pd) == -1) return -1; // read PDB file
	}
	else if (nfiletype == 1)
	{
		printf("Reading input file1 %s in Vesta XYZ export format ...\n", pdbfile);
		if (data_from_VestaXYZfile(pdbfile, &pd) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype == 2)
	{
		printf("Reading input file1 %s in Kirkland XYZ format ...\n", pdbfile);
		if (data_from_KirklandXYZfile(pdbfile, &pd) == -1) return -1; // read Kirkland XYZ file
	}
	//print_pdb_data(&pd);

	// sort entries by z coordinate in ascending order
	//pdb_bubbleSort2(&pd);

	// read input file2
	pdbdata pd1;
	pdbdata_init(&pd1);
	if (nfiletype1 == 1)
	{
		printf("Reading input file2 %s in Vesta XYZ export format ...\n", pdbfile1);
		if (data_from_VestaXYZfile(pdbfile1, &pd1) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype1 == 0)
	{
		printf("Reading input file2 %s in PDB format ...\n", pdbfile1);
		if (read_pdb(pdbfile1, &pd1) == -1) return -1; // read PDB file
	}
	else if (nfiletype1 == 2)
	{
		printf("Reading input file2 %s in Kirkland XYZ format ...\n", pdbfile1);
		if (data_from_KirklandXYZfile(pdbfile1, &pd1) == -1) return -1; // read Kirkland XYZ file
	}
	//print_pdb_data(&pd);


	//translate element symbols into atomic weights
	int i, j;
	int* ia = malloc(pd.natoms * sizeof(int));
	if (nfiletype != 2)
	{
		if (pdb_atomnumbers(&pd, ia))
		{
			printf("\n!!!Error encountered while finding atomic weight for a given element name!!!\n");
			return -1;
		}
	}
	else
		for (i = 0; i < pd.natoms; i++) ia[i] = pd.adata[i].serial;

	int* ia1 = malloc(pd1.natoms * sizeof(int));
	if (nfiletype1 != 2)
	{
		if (pdb_atomnumbers(&pd1, ia1))
		{
			printf("\n!!!Error encountered while finding atomic weight for a given element name!!!\n");
			return -1;
		}
	}
	else
		for (i = 0; i < pd1.natoms; i++) ia1[i] = pd1.adata[i].serial;


	printf("\nComparing atoms positions from %s test file with those in %s reference file ...", pdbfile, pdbfile1);
	printf("\n");
//	if (pd.natoms > pd1.natoms)
//	{
//		printf("\n!!!Error: the test file contains more atom entries than the reference file.!!!\n");
//		return -1;
//	}

	
	// open the output file early in order to start writing some results there as they are obtained
	FILE* ff = fopen(outfile, "wt");
	if (ff == NULL)
	{
		printf("\nERROR: cannot open %s file for writing!!!", outfile);
		return -2;
	}
	fprintf(ff, "For each atom from %s file, locating the closest atom in %s file.\n", pdbfile, pdbfile1);
	fprintf(ff, "Each entry contains the following data:\n");
	fprintf(ff, "1st column contains atomic weights (test data first, matched reference data second),\n");
	fprintf(ff, "2nd, 3rd, and 4th columns contain x, y and z atomic coordinates, respectively (test data first, matched reference data second),\n");
	fprintf(ff, "5th column contains original data for the test file and the L2 distance between atoms for the reference file data.\n");
	fprintf(ff, "%s", cfileinfo); // free-form file info line
	fprintf(ff, "\n");

	// for each atom in pd find the closest atom in pd1 
	pdbdata pd2;
	pdbdata_init(&pd2);
	pd2.natoms = pd.natoms * 2;
	pdballocate(&pd2);
	int* ia2 = malloc(pd.natoms * 2 * sizeof(int));
	int jmin;
	double r2, r2min;
	for (i = 0; i < pd.natoms; i++)
	{
		jmin = 0;
		r2min = (pd.adata[i].x - pd1.adata[0].x) * (pd.adata[i].x - pd1.adata[0].x) + (pd.adata[i].y - pd1.adata[0].y) * (pd.adata[i].y - pd1.adata[0].y) + (pd.adata[i].z - pd1.adata[0].z) * (pd.adata[i].z - pd1.adata[0].z);
		for (j = 1; j < pd1.natoms; j++)
		{
			r2 = (pd.adata[i].x - pd1.adata[j].x) * (pd.adata[i].x - pd1.adata[j].x) + (pd.adata[i].y - pd1.adata[j].y) * (pd.adata[i].y - pd1.adata[j].y) + (pd.adata[i].z - pd1.adata[j].z) * (pd.adata[i].z - pd1.adata[j].z);
			if (r2 < r2min) { r2min = r2;  jmin = j; }
		}
		pd2.adata[i * 2] = pd.adata[i];
		ia2[i * 2] = ia[i];
		pd2.adata[i * 2 + 1] = pd1.adata[jmin];
		ia2[i * 2 + 1] = ia1[jmin];
		// mark the found "nearest neighbour" atom in the reference file and flag if it was already marked
		if (pd1.adata[jmin].tempFactor == -777)
		{
			printf("\n!!! Duplicate localization for input atom no. %d, atom type = %d, atom positions = (%g, %g, %g).", i, ia[i], pd.adata[i].x, pd.adata[i].y, pd.adata[i].z);
			fprintf(ff, "\n!!! Duplicate localization for input atom no. %d, atom type = %d, atom positions = (%g, %g, %g).", i, ia[i], pd.adata[i].x, pd.adata[i].y, pd.adata[i].z);
			pd2.adata[i * 2].tempFactor = -777;
		}
		else
			pd1.adata[jmin].tempFactor = -777;
		if (ia1[jmin] == 1)
		{
			printf("\n??? Input atom no. %d, atom type = %d, atom positions = (%g, %g, %g), has been matched with a Hydrogen atom!!!.", i, ia[i], pd.adata[i].x, pd.adata[i].y, pd.adata[i].z);
			fprintf(ff, "\n??? Input atom no. %d, atom type = %d, atom positions = (%g, %g, %g), has been matched with a Hydrogen atom!!!.", i, ia[i], pd.adata[i].x, pd.adata[i].y, pd.adata[i].z);
		}
		// replace the occupancy entry in the reference array by the calculated distance to the nearest neighbour
		pd2.adata[i * 2 + 1].occupancy = sqrt(r2min);
	}

	// calculate average distance and std
	int nodup = 0;
	double adist = 0.0f;
	for (i = 0; i < pd.natoms; i++)
	{
		if (pd2.adata[i * 2].tempFactor == -777) continue;
		nodup++;
		adist += pd2.adata[i * 2 + 1].occupancy;
	}
	if (nodup == 0) adist = 0.0; else adist /= (double)nodup;
	double stddist = 0.0;
	for (i = 0; i < pd.natoms; i++)
	{
		if (pd2.adata[i * 2].tempFactor == -777) continue;
		stddist += (pd2.adata[i * 2 + 1].occupancy - adist) * (pd2.adata[i * 2 + 1].occupancy - adist);
	}
	if (nodup <= 1) stddist = 0.0; else stddist = sqrt(stddist / ((double)(nodup - 1)));
	printf("\nThere are %d non-duplicate atoms in the output file in total.", nodup);
	printf("\nAverage distance between non-duplicate test and template atoms = %g.", adist);
	printf("\nStandard deviation of the distance between non-duplicate test and template atoms = %g.", stddist);
	printf("\n");
	fprintf(ff, "\nThere are %d non-duplicate atoms in the output file in total.", nodup);
	fprintf(ff, "\nAverage distance between non-duplicate test and template atoms = %g.", adist);
	fprintf(ff, "\nStandard deviation of the distance between non-duplicate test and template atoms = %g.", stddist);
	fprintf(ff, "\n");

	// list all atoms from the reference file which have not been matched by any of the atoms from the test file
	for (j = 1; j < pd1.natoms; j++)
		if (pd1.adata[j].tempFactor != -777)
		{
			printf("\n@@@ Reference atom no. %d, atom type = %d, atom positions = (%g, %g, %g) has not been matched.", j, ia1[j], pd1.adata[j].x, pd1.adata[j].y, pd1.adata[j].z);
			fprintf(ff, "\n@@@ Reference atom no. %d, atom type = %d, atom positions = (%g, %g, %g) has not been matched.", j, ia1[j], pd1.adata[j].x, pd1.adata[j].y, pd1.adata[j].z);
		}
	
	// output to the target file
	printf("\nWriting output comparison file %s ...\n", outfile);
	fprintf(ff, "\n");
	for (i = 0; i < pd.natoms * 2; i++)
	{
		//if (ia[i] == 16) 
			// for even i the last column contains the occupancy or correlation coefficient from the first (original) file
			// for odd i the last column contains -1 times the radial distance between atoms i and i+1
			if (i % 2 == 0) fprintf(ff, "\nEntry no. %d", (int)(i / 2));
			fprintf(ff, "\n%d %f %f %f %f", ia2[i], pd2.adata[i].x, pd2.adata[i].y, pd2.adata[i].z, pd2.adata[i].occupancy);
	}
		
	fclose(ff);
	free(ia);

	printf("Finished!");

	return 0;
}