/*
      *** autosliceCMD.cpp ***

Front end for autoslic.cpp containing a command line user interface.

------------------------------------------------------------------------
Copyright 1998-2012 Earl J. Kirkland

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 
------------------------------------------------------------------------

  ANSI C and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  gcc -O -fopenmp -o autoslic autoslic.c slicelib.o
                       tiffsubs.o  cfpix.o -lfftw3f

  Transmit an electron wave through a specimen using the
  multislce method with automatic slicing.  Read in the (x,y,z)
  coordinates of the whole specimen and break into slices
  on-the-fly.

  started 24-july-1996 E. Kirkland
  working 19feb-1997 ejk
  last revised 19-feb-1997 ejk
  added look-up-table vzatomLUT() for 3X-4X increase 
        in speed 23-may-1997 ejk
  put bandwith limit inside trlayer() 1-oct-1997 ejk
  added Gaussian thermal displacements 1-oct-1997 ejk
  removed /sqrt(3) from Thermal rms displacements 
    to be consistent with Int'l X-ray tables 22-dec-1997 ejk
  corrected zmin/max error with thermal displac. 24-dec-1997 ejk
  fixed small aliasing problem 5-jan-1998 ejk
  added unit cell replication option and moved ReadXYZcoord()
    into slicelib.c  11-jan-1998 ejk
  added astigmatism and modify to use different set of
    random offsets on each illum. angle with partial coherence
         5-feb-1998 ejk
  fix typo in z range message with partial coherence and
    thermal vibrations 9-july-1998 ejk
  update memory allocation routines 19-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  fixed bug in zmin/zmax calculation in coherent mode
     (move to after sortByZ() - it was before ) 8-jan-2002 ejk
  add cross section option (in non-partial coherence mode only)
        27-may-2005 ejk
  convet to faster sortByZ() 8-feb-2006 ejk
  move sortbyz() to slicelib.c 5-sep-2006 ejk
  add echo on y position in pixels for xz mode 4-may-2007 ejk
  update data type of nxl,nyl to be consistent with new tiffsubs
     17-jul-2007 ejk
  move xz depthpix save to be after transmit+propagate to get a
     full slice and proper anti-aliasing and also be consisten
     with what you get doing it by hand  and increase possible
     slices output (nz was off by one) 24-jan-2008 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better (consistent with autostem) 23-mar-2008 ejk
  take small things out of loop in trlayer() 14-may-2008 ejk
  parameterize vzatomLUT() vs r^2 instead of r to avoid a lot of sqrt()
      calls (a little faster)  6-jun-2008 ejk
  move vzatomLUT() to slicelib.c  11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  add Cs5 (and Cs->Cs3) 15-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4
     and convert to 4 char TAB size formatting 21-feb-2010 ejk
  add parallel computing of a few parts 21-feb-2010 ejk
  start conversion to faster FFTW 24-feb-2010 ejk
  move some things into slicelibW.c to share 6-mar-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add option to average over many frozen phonon
      configurations 3-aug-2010 ejk
  add multipole aberrations to probe 12-may-2011 ejk
  start conversion to floatTIFF.cpp and C++ 28-may-2012 ejk
  working 3-jun-2012 ejk
  convert to cfpix/fftw class from raw fftw 13-nov-2012 to 21-nov-2012 ejk
  move calculation into a class with separate command line front end
      29-may-2013 ejk
  change RNG seed argument to referenece so it get updated for 
      successive calls 21-sep-2013 ejk
  fix %ld to %d format in printout of aslice.nillum  12-oct-2013 ejk
  add param[] for mode 15-oct-2013 ejk
  convert to streams and strings 22-mar-2014 ejk
  fix minor formatting issue 2-jun-2014 ejk

  ax,by,cz  = unit cell size in x,y
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 4 characters 
  
*/

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output

#include "XArray2D.h"
#include "XA_data.h"

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice
#include "floatTIFF.hpp"   // file I/O routines in TIFF format
#include "autoslic.hpp"    //  the calculation engine

#include "autosliccmd.h" //declaration of the function autosliccmd()

#define MANY_ABERR      //  define to include many aberrations

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

const int NSMAX= 1000;   // max number of slices
const int NCMAX= 1024;   // max characters in file names
const int NZMAX= 103;    // max atomic number Z


int autosliccmd(string params[numaslicpars])
{
    string filein, fileout, filestart, filebeam, filecross, cline, description;
  
    const char version[] = "2-jun-2014 (ejk)";

    int lstart=0, lpartl=0, lbeams=0, lwobble=0, lcross=0, nwobble=1;
    int ix, iy, iz, nx, ny, nzout, i, nslic0, islice,
        ndf, nbout(0), ib, ncellx, ncelly, ncellz, NPARAM;
    int nillum, nzbeams;
    int *hbeam(0), *kbeam(0);
    int natom, *Znum, done, status, multiMode(0);
    long  ltime;

    unsigned long iseed, iseed1;

    float v0, mm0, wavlen, rx, ry, ax, by, cz, pi,
        rmin, rmax, aimin, aimax, ctiltx, ctilty,
        acmin, acmax, Cs3, Cs5, df0, sigmaf, dfdelt(0), aobj(0),
        temperature, ycross(0), dx, dy;
   
    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;

    float *x, *y, *z, *occ, *wobble;
    float *param, *sparam;

    double timer, deltaz, vz;

	//@@@@@ start temporary code
	int nmode; // the switch between multislice(0), projection(1) and 1st Born(2) approximations
	float angle(0); // sample rotation angle in radians (in xz plane, i.e. around y axis)
	float xc(0), zc(0); // coordinates of the centre of rotation in Angstroms
	//@@@@@ end temporary code
        
    ofstream fp1;

    cfpix pix;		// to get results of calculation
    cfpix wave0;	// initial wavefunction (if used)
    cfpix depthpix;	// to get xz cross section results
    cfpix beams;	// to get valuse of requested beams during propagation

    floatTIFF myFile;	//  file input/output
    autoslic aslice;	// has the calculation engine

/*  echo version date and get input file name */

    cout << "autoslic(e) version dated " << version << endl;
    cout << "Copyright (C) 1998-2013 Earl J. Kirkland"  << endl;
    cout << "This program is provided AS-IS with ABSOLUTELY NO WARRANTY"  << endl;
    cout <<     " under the GNU general public license" << endl << endl;

    cout << "perform CTEM multislice with automatic slicing and FFTW" << endl;
#ifdef USE_OPENMP
    cout << "and multithreaded using openMP" << endl;
#endif
    cout <<  " "  << endl;

    pi = (float) (4.0 * atan( 1.0 ));
    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
    sparam = (float*) malloc1D( NPARAM, sizeof(float), "sparam" );
    for( ix=0; ix<NPARAM; ix++ ) param[ix] = 0.0F;

    //@@@@@ start temporary code
	//FILE* ff0 = fopen("autoslic.txt", "rt");
	char chaa[1024], cinarg[1024];
	//if (!ff0)
	//{
	//	cout << "!!!Cannot open input parameter file!!! Input any character to exit...";
	//	cin >> chaa;
	//	exit(-1);
	//}
	//cout << "Name of file with input atomic "
    //       << "coord. in x,y,z format:"  << endl;
    //cin >> filein;
	//filein = "1grl.xyz";
	if (sscanf(params[0].data(), "%s %s", chaa, cinarg) != 2)
	{
		cout << "!!!Error reading line 1 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	filein = cinarg;

/*  get simulation options */

    //cout << "Replicate unit cell by NCELLX,NCELLY,NCELLZ :" << endl;
    //cin >> ncellx >> ncelly >> ncellz;
    //if( ncellx < 1 ) ncellx = 1;
    //if( ncelly < 1 ) ncelly = 1;
    //if( ncellz < 1 ) ncellz = 1;
	//ncellx = 1; ncelly = 1; ncellz = 1;
	if (sscanf(params[1].data(), "%s %d %d %d", chaa, &ncellx, &ncelly, &ncellz) != 4)
	{
		cout << "!!!Error reading line 2 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}

    //cout << "Name of file to get binary output of multislice result:" << endl;
    //cin.getline( fileout, NCMAX );   //  ????
    //cin >> fileout ;
	//fileout = "1grl.tif";
	if (sscanf(params[2].data(), "%s %s", chaa, cinarg) != 2)
	{
		cout << "!!!Error reading line 3 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	fileout = cinarg;

    //lpartl = askYN("Do you want to include partial coherence");
	//lpartl = 0;
	if (sscanf(params[3].data(), "%s %d", chaa, &lpartl) != 2)
	{
		cout << "!!!Error reading line 4 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[4].data(), "%s %g %g", chaa, &acmin, &acmax) != 3)
	{
		cout << "!!!Error reading line 5 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[5].data(), "%s %g %g", chaa, &Cs3, &Cs5) != 3)
	{
		cout << "!!!Error reading line 6 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[6].data(), "%s %g %g %g", chaa, &df0, &sigmaf, &dfdelt) != 4)
	{
		cout << "!!!Error reading line 7 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[7].data(), "%s %g", chaa, &aobj) != 2)
	{
		cout << "!!!Error reading line 8 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[8].data(), "%s %d", chaa, &lstart) != 2)
	{
		cout << "!!!Error reading line 9 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[9].data(), "%s %s", chaa, cinarg) != 2)
	{
		cout << "!!!Error reading line 10 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	filestart = cinarg;
	if (sscanf(params[10].data(), "%s %g", chaa, &v0) != 2)
	{
		cout << "!!!Error reading line 11 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[11].data(), "%s %d %d", chaa, &nx, &ny) != 3)
	{
		cout << "!!!Error reading line 12 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[12].data(), "%s %g %g", chaa, &ctiltx, &ctilty) != 3)
	{
		cout << "!!!Error reading line 13 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[13].data(), "%s %lg", chaa, &deltaz) != 2)
	{
		cout << "!!!Error reading line 14 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[14].data(), "%s %d", chaa, &lbeams) != 2)
	{
		cout << "!!!Error reading line 15 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[15].data(), "%s %s", chaa, cinarg) != 2)
	{
		cout << "!!!Error reading line 16 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	filebeam = cinarg;
	if (sscanf(params[16].data(), "%s %d", chaa, &nbout) != 2)
	{
		cout << "!!!Error reading line 17 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[17].data(), "%s %d", chaa, &lwobble) != 2)
	{
		cout << "!!!Error reading line 18 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[18].data(), "%s %g", chaa, &temperature) != 2)
	{
		cout << "!!!Error reading line 19 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[19].data(), "%s %d", chaa, &nwobble) != 2)
	{
		cout << "!!!Error reading line 20 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[20].data(), "%s %d", chaa, &iseed1) != 2)
	{
		cout << "!!!Error reading line 21 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[21].data(), "%s %d", chaa, &lcross) != 2)
	{
		cout << "!!!Error reading line 22 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[22].data(), "%s %s", chaa, cinarg) != 2)
	{
		cout << "!!!Error reading line 23 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	filecross = cinarg;
	if (sscanf(params[23].data(), "%s %g", chaa, &ycross) != 2)
	{
		cout << "!!!Error reading line 24 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[24].data(), "%s %g", chaa, &angle) != 2)
	{
		cout << "!!!Error reading line 25 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[25].data(), "%s %g %g", chaa, &xc, &zc) != 3)
	{
		cout << "!!!Error reading line 26 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	if (sscanf(params[26].data(), "%s %d", chaa, &nmode) != 2)
	{
		cout << "!!!Error reading line 27 of input parameter array!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	//fclose(ff0);
	//cout << "Input parameter file has been read successfully!\n";

    //acmin = acmax = 0;
    if( lpartl == 1 ) {
        //cout << "Illumination angle min, max in mrad:" << endl;
        //cin >> acmin >> acmax;
        acmin  = acmin  * 0.001F;
        acmax  = acmax  * 0.001F;
        //cout << "Spherical aberration Cs3, Cs5(in mm.):" << endl;
        //cin >>  Cs3 >> Cs5;
		param[pCS]  = (float) ( Cs3*1.0e7 );
        param[pCS5] = (float) ( Cs5*1.0e7 );
        //cout << "Defocus, mean, standard deviation, and"
        //       " sampling size (in Angstroms) =" << endl;
        //cin >> df0 >> sigmaf >> dfdelt;
        param[pDEFOCUS] = (float) df0;
        param[pDDF] = (float) sigmaf;
        //cout << "Objective aperture (in mrad) =" << endl;
        //cin >> aobj;
		//@@@@@ end temporary code
        aobj = (float) fabs( aobj * 0.001F);
#ifdef MANY_ABERR 
        /*   get higher order aberrations if necessary */
        cout << "type higher order aber. name (as C32a, etc.) followed\n"
            << " by a value in mm. (END to end)" << endl;
        done = multiMode = 0;
        do{
            cin >> cline ;
            if( cline.compare( "END" ) == 0  ) {
                done = 1;
            } else {
                cin >> vz;
                status = readCnm( cline, param, vz );        
                if( status < 0 ) {
                    cout << "unrecognized aberration, exit..." << endl;
                    exit( EXIT_SUCCESS );
                } else multiMode = 1;
            }
        } while( !done );

#endif
        lstart = 0;
    } else {
        cout << "NOTE, the program image must also be run." << endl;
		//@@@@@ start temporary code
        //lstart = askYN("Do you want to start from previous result");
		//lstart = 0;
		//@@@@@ end temporary code
    }
	
	//@@@@@ start temporary code
    if ( lstart == 1 ) {
        //cout << "Name of file to start from:" << endl;
        //cin >> filestart;
    } else {
        //cout << "Incident beam energy in kev:" << endl;
        //cin >> v0;
		//v0 = 200;
        //cout << "Wavefunction size in pixels, Nx,Ny:" << endl;
        //cin >> nx >> ny;
		//nx = ny = 1024;
		//@@@@@ end temporary code
    }

	//@@@@@ start temporary code
    //cout << "Crystal tilt x,y in mrad.:" << endl;
    //cin >> ctiltx >> ctilty;
	//ctiltx = ctilty = 0;
	//@@@@@ end temporary code
    ctiltx = ctiltx /1000;
    ctilty = ctilty /1000;

    /*  remember that the slice thickness must be > atom size
        to use projected atomic potential */
	//@@@@@ start temporary code
    //cout << "Slice thickness (in Angstroms):" << endl;
    //cin >> deltaz;
	//deltaz = 1.3575;
	//@@@@@ end temporary code
    if( deltaz < 1.0 ) {
        cout << "WARNING: this slice thickness is probably too thin"
            << " for autoslice to work properly." << endl;
    }

    if( lpartl == 0 ) {
		//@@@@@ start temporary code
        //lbeams = askYN("Do you want to record the (real,imag) value\n"
        //" of selected beams vs. thickness");
		//lbeams = 0;
		//@@@@@ end temporary code
        if( lbeams == 1 ) {
            //cout << "Name of file for beams info:" << endl;
            //cin >> filebeam;
            //cout << "Number of beams:" << endl;
            //cin >> nbout;
            if( nbout<1 ) nbout = 1;
            hbeam = (int*) malloc1D( nbout, sizeof(int), "hbeam" );
            kbeam = (int*) malloc1D( nbout, sizeof(int), "kbeam" );
            for( ib=0; ib<nbout; ib++) {
                cout << "Beam " << ib+1 << ", h,k=" << endl;
                cin >> hbeam[ib] >> kbeam[ib];
            }
        }
    }

	//@@@@@ start temporary code
	//lwobble = askYN("Do you want to include thermal vibrations");
	//lwobble = 0;
	//@@@@@ end temporary code
    if( lwobble == 1 ) {
        //cout << "Type the temperature in degrees K:" << endl;
        //cin >> temperature ;
        //cout << "Type number of configurations to average over:" << endl;
        //cin >>  nwobble;
        if( nwobble < 1 ) nwobble = 1;
        /* get random number seed from time if available 
            otherwise ask for a seed */
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            //cout << "Type initial seed for random number generator:" << endl;
            //cin >>  iseed;
			iseed = iseed1;
        } else
            cout << "Random number seed initialized to " << iseed << endl;
    } else temperature = 0.0F;

    if( lpartl == 0 ) {
		//@@@@@ start temporary code
        //lcross = askYN("Do you want to output intensity vs. depth cross section");
		//lcross = 0;
		//@@@@@ end temporary code
        if( lcross == 1 ){
            //cout << "Type name of file to get depth profile image:" << endl;
            //cin >> filecross;
            //cout << "Type y position of depth cross section (in Ang.):" << endl;
            //cin >> ycross;
        }
    }

/* start timing the actual computation just for fun */

    timer = cputim();
#ifdef USE_OPENMP
    walltimer = walltim();  /* wall time for openMP */
#endif

/*  get starting value of transmitted wavefunction if required
   (this can only be used in coherent mode)
    remember to save params for final output pix  */

    if ( lstart == 1 ) {
        if( myFile.read( filestart.c_str() ) != 1 ) {
            cout << "Cannot open input file: " << filestart << endl; 
            exit( 0 );
        }

        if( myFile.getnpix() != 2 ) {
           cout << "Input file " << filestart << " must be complex, can't continue." << endl;
           exit( 0 );
        }

        nx =  myFile.nx();
        ny =  myFile.ny();

        nx = nx/2;
        wave0.resize( nx, ny );

        //  save starting pix for later
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
                wave0.re(ix,iy) = myFile(ix,iy);
                wave0.im(ix,iy) = myFile(ix+nx,iy);
        }

        //  save parameters to verify successive images are same size etc.
        for( i=0; i<NPARAM; i++) sparam[i] = myFile.getParam( i );

        ax = sparam[pDX] * nx;
        by = sparam[pDY] * ny;
        v0     = sparam[pENERGY];
        nslic0 = (int) sparam[pNSLICES];
        cout << "Starting pix range " << sparam[pRMIN] << " to " << sparam[pRMAX] 
              << " real\n" << "           " << sparam[pIMIN] << " to " 
              << sparam[pIMAX] << " imag" << endl;
        cout << "Beam voltage = " << v0 << " kV" << endl;
        cout << "Old crystal tilt x,y = " << 1000.*sparam[pXCTILT] <<  ", " 
            << 1000.*sparam[pYCTILT] << " mrad" << endl;

    } else nslic0 = 0;     /* end if( lstart...) */

/*  calculate relativistic factor and electron wavelength */

    mm0 = 1.0F + v0/511.0F;
    wavlen = (float) wavelength( v0 );
    cout << "electron wavelength = " << wavlen << " Angstroms" << endl;

/*  read in specimen coordinates and scattering factors */

    natom = ReadXYZcoord( filein.c_str(), ncellx, ncelly, ncellz,
        &ax, &by, &cz, &Znum, &x, &y, &z, &occ, &wobble,
        description );

	if (!natom)
	{
		cout << "!!!Zero atomic coordinates read. Input any key to exit...!!!";
		cin >> chaa;
		exit(-1);
	}
	
    cout << natom << " atomic coordinates read in"  << endl;
    cout << description << endl;

    cout <<"Size in pixels Nx, Ny= " << nx << " x " << ny << " = " << nx*ny 
        << " beams" << endl;
    cout <<"Lattice constant a,b = " << ax << ", " << by << endl;


	//@@@@@ start temporary code
	// rotate the sample if necessary
	if (angle)
	{
		for (size_t k = 0; k < natom; k++)
		{
			x[k] = xc + (x[k] - xc) * cos(angle) + (z[k] - zc) * sin(angle);
			z[k] = zc + (-x[k] + xc) * sin(angle) + (z[k] - zc) * cos(angle);
		}
	}
	//@@@@@ end temporary code
	
	/*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
	//@@@@@ start temporary code
	// force square dimensions in (xz) plane, assuming that the x-size is larger or equal to the z-size
	if (zmin < xmin || zmax > xmax)
	{
		cout << "!!!Error: zmin < xmin or zmax > xmax in the XYZ file!!! Input any character to exit...";
		cin >> chaa;
		exit(-1);
	}
	zmin = xmin;
	zmax = xmax;
	//@@@@@ end temporary code
    cout << "Total specimen range is\n" 
        << xmin << " to " << xmax << " in x\n"
            << ymin << " to " << ymax << " in y\n"
        << zmin << " to " << zmax << " in z" << endl;
    if( lwobble == 1 )
        cout << "Range of thermal rms displacements (300K) = "
            << wmin << " to " << wmax << endl;

    // ---------  setup calculation -----
    //   set calculation flags
    aslice.lbeams = lbeams;
    aslice.lcross = lcross;
    aslice.lpartl = lpartl;
    aslice.lstart = lstart; 
    aslice.lwobble = lwobble;

    //   set calculation parameters (some already set above)
    param[ pAX ] = ax;			// supercell size
    param[ pBY ] = by;
    param[ pNX ] = (float) nx;
    param[ pNY ] = (float) ny;
    param[pENERGY]   =  v0;
    param[ pDELTAZ ] = (float) deltaz;	// slice thickness
    param[ pOAPERT ] = aobj;
    param[ pXCTILT ] = ctiltx;		// crystal tilt
    param[ pYCTILT ] = ctilty;
    param[pCAPERT] = acmax;		// condencer angles
    param[pCAPERTMIN] = acmin;
    param[ pTEMPER ] = (float) fabs( temperature );
    param[ pNWOBBLE ] = (float) nwobble;	//  number config. to average
    param[pWAVEL] = wavlen;			//  probably recal. autoslice::calculate()

    param[pMODE] = 6;  // save mode = autoslic

    if ( lpartl == 1 ) {
            param[pDEFOCUS] = df0;
            param[pOAPERT] = aobj;
            param[pDDF] = sigmaf;
            param[pCAPERT] = acmax;

            rx = 1.0F/ax;
            ry = 1.0F/by;
            cout << "Illumination angle sampling (in mrad) = "
                << 1000.*rx*wavlen << ", " << 1000.*ry*wavlen << "\n" << endl;
    }

    // ------- iterate the multislice algorithm proper -----------

    aslice.calculate( pix, wave0, depthpix, param, multiMode, natom, &iseed,
                Znum, x,y,z,occ,wobble, beams, hbeam, kbeam, nbout, ycross, dfdelt, nmode);
 
    if( lpartl == 1 ) {         //    with partial coherence
        nillum = aslice.nillum;
        cout << "Total number of illumination angle = "
                << nillum << endl;
        ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );  // recal same value in class
        cout << "Total number of defocus values = " << 2*ndf+1 << endl;
    }
   
    else if( lbeams ==1 ) {
            fp1.open( filebeam.c_str() );
            if( fp1.bad() ) {
                cout << "can't open file " << filebeam << endl;
                exit(0);
            }
            fp1 << " (h,k) = " ;
            for(ib=0; ib<nbout; ib++)
                fp1 << " (" << hbeam[ib] << "," << kbeam[ib] << ")";
            fp1 << endl;
            nzbeams = beams.ny();
            fp1 << "nslice, (real,imag) (real,imag) ...\n" << endl;
            for( islice=0; islice<nzbeams; islice++) {
                fp1 << setw(5) << islice;
                for( ib=0; ib<nbout; ib++) //  setprecision(4)
                    fp1 << "  " << setw(10) << beams.re(ib,islice)		//????? "%10.6f %10.6f",
                        << "  " << setw(10) << beams.im(ib,islice); 
                fp1 << endl;
            }
        fp1.close();

    } // end else 

/*  ------------------------------------------------------
    output results and find min and max to echo
    remember that complex pix are stored in the file in FORTRAN
        order for compatibility */

    pix.findRange( rmin, rmax, aimin, aimax );

    param[pRMAX]  = rmax;
    param[pIMAX]  = aimax;
    param[pRMIN]  = rmin;
    param[pIMIN]  = aimin;

    param[pDX] = dx = (float) ( ax/((float)nx) );
    param[pDY] = dy = (float) ( by/((float)ny) );

    //param[pNSLICES] = 0.0F;  /* ??? */
	//@@@@@ start temporary code
	//GRD file output
	xar::XArray2D<double> inten(ny, nx);
	IXAHWave2D* ph2new = CreateWavehead2D();
	ph2new->SetData(wavlen * 1.e-4, ymin * 1.e-4, ymax * 1.e-4, xmin * 1.e-4, xmax * 1.e-4);
	inten.SetHeadPtr(ph2new);
	for (ix = 0; ix < nx; ix++)
		for (iy = 0; iy < ny; iy++)
			inten[iy][ix] = pix.re(ix, iy) * pix.re(ix, iy) + pix.im(ix, iy) * pix.im(ix, iy);
	xar::XArData::WriteFileGRD(inten, fileout.c_str(), xar::eGRDBIN);

 /*   for( ix=0; ix<NPARAM; ix++ ) myFile.setParam( ix, param[ix] );

    if ( lpartl == 1 ) {
        myFile.resize( nx, ny );
        myFile.setnpix( 1 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
            myFile(ix,iy) = pix.re(ix,iy);
    } else {
        myFile.resize( 2*nx, ny );
        myFile.setnpix( 2 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            myFile(ix,iy)    = pix.re(ix,iy);
            myFile(ix+nx,iy) = pix.im(ix,iy);
        }
    }

    i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );
    if( i != 1 ) cout << "autoslice cannot write TIF file " << fileout << endl;
    cout << "pix range " << rmin << " to " << rmax << " real,\n" <<
            "          " << aimin << " to " << aimax << " imag" << endl;
*/
	//@@@@@ end temporary code

    /* ----- output depth cross section if requested ------- */
    if( lcross == 1 ){
        depthpix.findRange( rmin, rmax, aimin, aimax );
        myFile.setParam( pRMAX, rmax );
        myFile.setParam( pIMAX, 0.0F );
        myFile.setParam( pRMIN, rmin );
        myFile.setParam( pIMIN, 0.0F );
        myFile.setParam( pDY, dy = (float) ( deltaz ) );

        nzout = depthpix.ny();
        myFile.resize( nx, nzout );
        myFile.setnpix( 1 );
        for( ix=0; ix<nx; ix++) for( iz=0; iz<nzout; iz++) {
            myFile(ix,iz) = depthpix.re(ix,iz);
        }
        i = myFile.write( filecross.c_str(), rmin, rmax, aimin, aimax, dx, dy );

        if( i != 1 ) cout << "autoslice cannot write TIF file "
                << filecross << endl;
        cout << "depth pix range " << rmin << " to " << rmax << " real" << endl;
    }

    cout << "Total CPU time = " << cputim()-timer << " sec." << endl;
#ifdef USE_OPENMP
    cout << "wall time = " << walltim() - walltimer << " sec." << endl;
#endif
	//@@@@@ start temporary code
	//char a;
	//cout << "\nPress any key to exit ...";
	//cin >> a;
	//@@@@@ end temporary code
    return 0;

} /* end main() */



