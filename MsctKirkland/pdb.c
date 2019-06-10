
/*
 *  pdb.h   : read a pdb file and store data in a structure
 *  A. Martin January 2013
 *
 */

#include "pdb.h"



/*
 * Routines to write:
 * read the pdb
 * check pdb for number of atoms (I have the feeling I wanted to check something else in this routine too
 * Initialize the pdbdata struct
 * Initilize pdbdata arrays (after checking for the number of atoms)
 */


/*
 *  Initilize values of a pdbdata structure
 *  Values that require pdbcheck are set to 0
 *  Arrays that require pdbcheck are set to NULL
 */
void pdbdata_init( pdbdata * pd)
{

  pd->natoms = 0;
  pd->adata = NULL;
}



/*
 *  Read all the atom information into a pdbdata structure
 */
void read_pdb( char * pdbname, pdbdata * pd )
{

  int natoms;

  /*
   *  do a pdb check to get the number of atoms
   */
  pdb_check( pdbname, &natoms);
  
  /*
   *  Allocate memory
   */
  printf("Number of atoms %d\n", natoms);
 
  pd->natoms = natoms;

  pdballocate( pd );


  /*
   * Read atom data into the pdbdata
   */
  pdb_atomdata_from_file(pdbname, pd);

  /*
   *  Print the data in pd
   */
  //print_pdb_data( pd );

}



/*
 *   Put info for each atom line int pdbdata
 */
void pdb_atomdata_from_file(char * pdbfname, pdbdata * pd)
{

  FILE * file = fopen ( pdbfname, "r" );
  char line[256];
  char line_temp[256];
  int i=0;
  line[0] = '\0';

  int count =0;

  if (file != NULL){
    while( fgets(line, sizeof line, file) != NULL) {
     
      //fputs ( line, stdout ); /* write the line */
      strcpy( line_temp, line );

      if (isATOMline(line_temp) == 1) {
	readATOMline( line, &(pd->adata[count]) );
	count ++;
      }
      i++;
      //printf("Line %d ; ATOM count %i\n", i, count);  
      
    }
   
    fclose(file);
  } 
  else {
    printf("pdb file was not found (for pdb_check).\n");
  }
}




/*
 * Counts the number of ATOM lines in the file
 */

void pdb_check(char * pdbfname, int * natoms){

  FILE * file = fopen ( pdbfname, "r" );
  char line[256];
  int i=0;
  line[0] = '\0';

  int count =0;

  if (file != NULL){
    while( fgets(line, sizeof line, file) != NULL) {
     
      //fputs ( line, stdout ); /* write the line */
      
      if (isATOMline(line) == 1) {
	count ++;
      }
      i++;
      //printf("Line %d ; ATOM count %i\n", i, count);
           
    }
    
    *natoms = count;

    fclose(file);
  } 
  else {
    printf("pdb file was not found (for pdb_check).\n");
  }
}


/*
 *  allocate an array of atomdata structures in a pdbdata structure
 */
void pdballocate(pdbdata * pd)
{

  
  if (pd->adata != NULL){
    free(pd->adata);
    pd->adata = NULL;
    
    } 
 
  pd->adata = malloc( pd->natoms * sizeof(atomdata) );
  int i;
  for (i=0;i<pd->natoms;i++){
    initialise_atomdata( &(pd->adata[i]) );
  }

}

/*
 *  free data for pdbdata.adata
 */
void pdbfree(pdbdata * pd)
{

  free(pd->adata);
  pd->adata = NULL;
 
}

void initialise_atomdata( atomdata * ad ){

  ad->rname[0] = '\0';
  ad->atomName[0] = '\0';
  ad->resName[0] = '\0';
  ad->element[0] = '\0';
  ad->charge[0] = '\0';
  ad->serial = 0;
  ad->altLoc = '\0';
  ad->chainID = '\0';
  ad->x = 0;
  ad->y = 0;
  ad->z = 0;
  ad->occupancy = 0.0;
  ad->tempFactor = 0.0;
}


/*
 *  Returns 1 if the first word of the line is ATOM, and 0 otherwise
 */

int isATOMline(char * line)
{

  char * token;
  token = "\0";

  token = strtok(line," =:\n"); 

  if (strcmp(token,"ATOM") == 0) {
    return 1;
  }
  else {
    return 0;
  }
 
}


/*
 *  Reads a line by tokenizing using strtok
 *  If it is an ATOM line, it prints each token
 */

void readline( char * line){

  char * token;
  token = "\0";

  token = strtok(line," =:\n"); 

  if (strcmp(token,"ATOM") == 0) {
  
    while ( token != NULL){
      
      fputs ( token, stdout ); /* write the next part of the line */
      fputs ( " ", stdout ); /* space between tokens */
      token = strtok(NULL," =:\n"); 
    }
    fputs ( "\n", stdout ); /* new line */

  }
 

}




/*
 *  Reads an ATOM line and stores relevant information
 */
void readATOMline( char * line, atomdata * ad)
{

  int i;
  int pos, len;

  /*
   *  Copy the record name
   */
  pos = 0;
  len = 6;
  char rname[7];
  for (i=0;i<len;i++){
    rname[i] = line[pos+i];
  }
  pos += len;
  rname[len] = '\0';
  strcpy( ad->rname, rname );

  /*
   *  Copy the atom serial number
   */
  len = 5;
  char serial[6];
  for (i=0;i<len;i++){
    serial[i] = line[pos+i];
  }
  pos += len;
  serial[len] = '\0';
  ad->serial = atoi(serial);


  /*
   *  12th position is blank
   */
  pos ++;


  /*
   *  Copy the atom name
   */
  len = 4;
  char aname[5];
  for (i=0;i<len;i++){
    aname[i] = line[pos+i];
  }
  pos += len;
  aname[len] = '\0';
  strcpy( ad->atomName, aname );

  /*
   *  Copy the alternative location indicator
   */
  ad->altLoc = line[pos];
  pos ++;


  /*
   *  Copy the residual name
   */
  len = 3;
  char resName[4];
  for (i=0;i<len;i++){
    resName[i] = line[pos+i];
  }
  pos += len;
  resName[len] = '\0';
  strcpy( ad->resName, resName );

  /*
   *  21st position is blank
   */
  pos ++;
  
  /*
   *  Copy the chain ID
   */
  ad->chainID = line[pos];
  pos ++;
  

  /*
   *  Copy the residue sequence number
   */
  len = 4;
  char resSeq[5];
  for (i=0;i<len;i++){
    resSeq[i] = line[pos+i];
  }
  pos += len;
  resSeq[len] = '\0';
  ad->resSeq = atoi( resSeq );


  /*
   *  Code for the insertion of residues
   */
  ad->iCode = line[pos];
  pos ++;

  /*
   *  28-30 is blank
   */
  pos += 3;



  /*
   *  the x co-ordinate of the atom
   */
  len = 8;
  char x[9];
  for (i=0;i<len;i++){
    x[i] = line[pos+i];
  }
  pos += len;
  x[len] = '\0';
  ad->x = atof( x );


  /*
   *  the y co-ordinate of the atom
   */
  len = 8;
  char y[8];
  for (i=0;i<len;i++){
    y[i] = line[pos+i];
  }
  pos += len;
  ad->y = atof( y );

  
  /*
   *  the z co-ordinate of the atom
   */
  len = 8;
  char z[8];
  for (i=0;i<len;i++){
    z[i] = line[pos+i];
  }
  pos += len;
  ad->z = atof( z );


  /*
   *  occupancy
   */
  len = 6;
  char occ[6];
  for (i=0;i<len;i++){
    occ[i] = line[pos+i];
  }
  pos += len;
  ad->occupancy = atof( occ );


  /*
   *  temperature factor
   */
  len = 6;
  char tempfact[6];
  for (i=0;i<len;i++){
    tempfact[i] = line[pos+i];
  }
  pos += len;
  ad->tempFactor = atof( tempfact );


  /*
   *  67-76 is blank
   */
  pos += 10;


  /*
   *  element
   */
  len = 2;
  char element[3];
  for (i=0;i<len;i++){
    element[i] = line[pos+i];
  }
  pos += len;
  element[len] = '\0';
  strcpy( ad->element, element );



  /*
   *  charge
   */
  len = 2;
  char charge[3];
  for (i=0;i<len;i++){
    charge[i] = line[pos+i];
  }
  pos += len;
  charge[len] = '\0';
  strcpy( ad->charge, charge );


}


/*
 * print all the information contained in pdfdata.adata
 */
void print_pdb_data( pdbdata * pd )
{

  int i;

  for (i=0;i<pd->natoms;i++){
    
    printf("just resName %s\n", pd->adata[i].resName);
    
    printf("atom %d : %s | %d | %s | %c | %s | %c | %d | %c | %f | %f | %f | %f | %f | %s | %s\n",
	   i, pd->adata[i].rname, pd->adata[i].serial, pd->adata[i].atomName, 
	   pd->adata[i].altLoc, pd->adata[i].resName, pd->adata[i].chainID, 
	   pd->adata[i].resSeq, pd->adata[i].iCode, 
	   pd->adata[i].x, pd->adata[i].y, pd->adata[i].z, 
	   pd->adata[i].occupancy, pd->adata[i].tempFactor, 
	   pd->adata[i].element, pd->adata[i].charge);

  }

}


void pdb_chemical_composition( pdbdata * pd ){

  int i;
  int n[6];
  for(i=0;i<6;i++){
    n[i] = 0;
  }

  for (i=0;i<pd->natoms;i++){
    //printf("%s\n", pd->adata[i].element);
    if(pd->adata[i].element[1] == 'H'){
      n[0]++;
    }
    if(pd->adata[i].element[1] == 'C'){
      n[1]++;
    }

    if(pd->adata[i].element[1] == 'N'){
      n[2]++;
    }
    if(pd->adata[i].element[1] == 'O'){
      n[3]++;
    }
    if(pd->adata[i].element[1] == 'P'){
      n[4]++;
    }
    if(pd->adata[i].element[1] == 'S'){
      n[5]++;
    }
  }

  for(i=0;i<6;i++){
    printf("%d ", n[i]);
  }
  printf("\n");
}
