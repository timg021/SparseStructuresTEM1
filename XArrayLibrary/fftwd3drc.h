//Header fftwd3drc.h
//
//
//	HEADER FILE TITLE:
//
//	C++ wrapper for FFTW double-precision (double) FFT routines of real-valued 3D arrays (forward "real-to-complex" FFT and inverse "complex-to-real" FFT)
//
//	COPYRIGHT:
//
//		TEG 2019
//
//
/*!
	\file		XArray.h
	\brief		C++ wrapper for FFTW double-precision (double) FFT routines of real-valued 3D arrays
	\par		Description:
		This is C++ wrapper for FFTW double-precision (double) FFT routines of real-valued 3D arrays (forward "real-to-complex" FFT and inverse "complex-to-real" FFT).

*/
#if !defined FFTWD3DRC_H
#define FFTWD3DRC_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//

#include "fftw3.h"
#include "XArray3D.h"


class Fftwd3drc
{
	// Enumerators
	// Structures
	// Constructors
public:
	//! Default constructor
	//Fftwd3drc() { if (!IsEmpty()) Cleanup(); uiflag = FFTW_ESTIMATE;  }
	//! Constructor allocating internal storage and creating "FFT plans"
	Fftwd3drc(int nz1, int ny1, int nx1, bool bMeasure = false)
	{ 
		if (!IsEmpty()) Cleanup();
		uiflag = bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE;
		nz = nz1; ny = ny1; nx = nx1;
		pin = (double*)fftw_malloc(sizeof(double) * GetNr());
		pout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * GetNc());
		aplan = fftw_plan_dft_r2c_3d(nz, ny, nx, pin, pout, uiflag);
		bplan = fftw_plan_dft_c2r_3d(nz, ny, nx, pout, pin, uiflag);
	}
private:
	//! Copy constructor (declared private to prohibit copying)
	Fftwd3drc(const Fftwd3drc& other) {  }
public:
	//! Destructor, destroys the plans and deallocates internal storage
	~Fftwd3drc() { if (!IsEmpty()) Cleanup(); }

	// Operators
private:
	//! Assignment operator (declared private to prohibit copying)
	Fftwd3drc& operator=(const Fftwd3drc& other) { return *this; }

	// Attributes
public:

	// Operations
public:
	// Member functions
	// Get full size of the real array
	inline int GetNr() { return nz * ny * nx; }
	// Get the last dimension of the complex array
	inline int GetNx2() { return int(nx / 2 + 1); }
	// Get full size of the complex array
	inline int GetNc() { return nz * ny * GetNx2(); }
	// Checks if the object is empty
	inline bool IsEmpty() { return (nz == 0 && ny == 0 && nx == 0 && pin == 0 && pout == 0 && aplan == 0 && bplan == 0); }
	// Empties the object: may be called if one wants to release the memory before the destructor is called
	void Cleanup();
	// Get/Set arrays
	inline double* GetReal() { return pin; }
	void GetRealXArray3D(xar::XArray3D<double>& aaa);
	void SetRealXArray3D(xar::XArray3D<double> aaa);
	inline fftw_complex* GetComplex() { return pout; }
	void GetComplexXArray3D(xar::XArray3D<xar::dcomplex>& aaa);
	void SetComplexXArray3D(xar::XArray3D<xar::dcomplex> aaa);
	// Forward FFT
	void ForwardFFT() { fftw_execute(aplan); }
	// Inverse FFT
	void InverseFFT() { fftw_execute(bplan); }
	// Print internal arrays
	void PrintRealArray(const char* message);
	void PrintComplexArray(const char* message);

private:
	// Member variables
	int nz; // z dimension (= Dim1 for XArray3D, = n0 for fftw)
	int ny; // y dimension (= Dim2 for XArray3D, = n1 for fftw)
	int nx; // x dimension (= Dim3 for XArray3D, = n2 for fftw)
	unsigned int uiflag; // FFTW_MEASURE or FFTW_ESTIMATE

	double* pin; // pointer to the real array
	fftw_complex* pout; // pointer to the complex array

	fftw_plan aplan; // forward plan
	fftw_plan bplan; // inverse plan
};


#endif	// FFTWD3DRC_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//