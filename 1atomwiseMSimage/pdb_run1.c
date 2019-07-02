#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "pdb.h"


int pdb_read(pdbdata* ppd, int nfiletype, char* pdbfile)
{
	//pdbdata pd;
	printf("\nReading atomic coordinates file ...\n");
	pdbdata_init(ppd);
	if (nfiletype)
	{
		printf("\nReading input XYZ file in Vesta export format ...\n");
		data_from_VestaXYZfile(pdbfile, ppd); // read Vesta export XYZ file
	}
	else
	{
		printf("\nReading input PDB file ...\n");
		read_pdb(pdbfile, ppd); // read PDB file
	}
	//print_pdb_data(&pd);

	printf("Finished!");
	return 0;
}


int xyz_Kirck_write1(int i, pdbdata* ppd, char* outfile, char* cfileinfo, double ctblength, double xminx0, double yminy0, double zminz0)
{
	if (i < 0 || i >= ppd->natoms)
	{
		printf("!!!Bad atom index %d in xyz_Kirck_write1() function.", i);
		return -1;
	}

	int ia, j;

	printf("\nWriting output XYZ file in Kirkland's format ...\n");
	FILE* ff = fopen(outfile, "wt");
	fprintf(ff, "%s\n", cfileinfo); // free-form file info line
	fprintf(ff, "%f %f %f\n", ctblength, ctblength, ctblength); // ctblength is written as the unit cell dimension

	// translate the element symbol into atomic number (weight)
	ia = 0; j = 0;
	while (isspace(ppd->adata[i].element[j])) j++; // skip leading white spaces
	if (j > 2)
	{
		printf("!!!Too many white spaces %s!!!", ppd->adata[i].element);
		return -1;
	}
	if (ppd->adata[i].element[j] == 'H') ia = 1; // H = hydrogen
	else if (ppd->adata[i].element[j] == 'C')
	{
		ia = 6; // C = carbon
		if (ppd->adata[i].element[j + 1] == 'A') ia = 20; // CA = calcium
		else if (ppd->adata[i].element[j + 1] == 'O') ia = 27; // CO = cobalt
	}
	else if (ppd->adata[i].element[j] == 'N')
	{
		ia = 7; // N = nitrogen
		if (ppd->adata[i].element[j + 1] == 'A') ia = 11; // NA = sodium
		else if (ppd->adata[i].element[j + 1] == 'I') ia = 28; // NI = nickel
	}
	else if (ppd->adata[i].element[j] == 'O') ia = 8; // O = oxygen
	else if (ppd->adata[i].element[j] == 'P') ia = 15; // P = posphorus
	else if (ppd->adata[i].element[j] == 'S') ia = 16; // S = sulphur
	else if (ppd->adata[i].element[j] == 'F') ia = 26; // assuming FE = ferrum
	else if (ppd->adata[i].element[j] == 'M') ia = 25; // assuming MN = manganese
	else if (ppd->adata[i].element[j] == 'Z') ia = 30; // assuming ZN = zinc
	else if (ia == 0)
	{
		printf("!!!Unknown atom type %s!!!", ppd->adata[i].element);
		return -1;
	}
	// TEG: we make all coordinates non-negative, as it seems that Kirkland's software expects this,
	// and we also centre each of XYZ coordinates within the interval [0, ctblength], so that the sample can remain within the cube during rotation
	fprintf(ff, "%d %f %f %f %f\n", ia, ppd->adata[i].x + xminx0, ppd->adata[i].y + yminy0, ppd->adata[i].z + zminz0, ppd->adata[i].occupancy);

	fprintf(ff, "%i\n\n\n", -1);
	fclose(ff);

	printf("Finished!");

	return 0;
}