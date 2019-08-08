//Header fftwf3drc.h
//
//
//	HEADER FILE TITLE:
//
//	C++ wrapper for FFTW single-precision (float) FFT routines of real-valued 3D arrays (forward "real-to-complex" FFT and inverse "complex-to-real" FFT)
//
//	COPYRIGHT:
//
//		TEG 2019
//
//
/*!
	\file		XArray.h
	\brief		C++ wrapper for FFTW single-precision (float) FFT routines of real-valued 3D arrays
	\par		Description:
		This is C++ wrapper for FFTW single-precision (float) FFT routines of real-valued 3D arrays (forward "real-to-complex" FFT and inverse "complex-to-real" FFT).

*/
#if !defined FFTWF3DRC_H
#define FFTWF3DRC_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//

#include "fftw3.h"
#include "XArray3D.h"


class Fftwf3drc
{
	// Enumerators
	// Structures
	// Constructors
public:
	//! Default constructor
	Fftwf3drc() { if (!IsEmpty()) Cleanup(); uiflag = FFTW_ESTIMATE;  }
	//! Constructor allocating internal storage and creating "FFT plans"
	Fftwf3drc(int nx1, int ny1, int nz1, bool bMeasure = false)
	{ 
		if (!IsEmpty()) Cleanup();
		uiflag = bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE;
		nx = nx1; ny = ny1; nz = nz1;
		pin = (float*)fftwf_malloc(sizeof(float) * GetNr());
		pout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * GetNc());
		aplan = fftwf_plan_dft_r2c_3d(nx, ny, nz, pin, pout, uiflag);
		bplan = fftwf_plan_dft_c2r_3d(nx, ny, nz, pout, pin, uiflag);
	}
private:
	//! Copy constructor (declared private to prohibit copying)
	Fftwf3drc(const Fftwf3drc& other) {  }
public:
	//! Destructor, destroys the plans and deallocates internal storage
	~Fftwf3drc() { if (!IsEmpty()) Cleanup(); }

	// Operators
private:
	//! Assignment operator (declared private to prohibit copying)
	Fftwf3drc& operator=(const Fftwf3drc& other) { return *this; }

	// Attributes
public:

	// Operations
public:
	// Member functions
	// Get full size of the real array
	inline int GetNr() { return nx * ny * nz; }
	// Get the last dimension of the complex array
	inline int GetNz2() { return int(nz / 2 + 1); }
	// Get full size of the complex array
	inline int GetNc() { return nx * ny * GetNz2(); }
	// Checks if the object is empty
	inline bool IsEmpty() { return (nx == 0 && ny == 0 && nz == 0 && pin == 0 && pout == 0 && aplan == 0 && bplan == 0); }
	// Empties the object: may be called if one wants to release the memory before the destructor is called
	void Cleanup();
	// Get/Set arrays
	inline float* GetReal() { return pin; }
	void GetRealXArray3D(xar::XArray3D<float>& aaa);
	void SetRealXArray3D(xar::XArray3D<float> aaa);
	inline fftwf_complex* GetComplex() { return pout; }
	void GetComplexXArray3D(xar::XArray3D<xar::fcomplex>& aaa);
	void SetComplexXArray3D(xar::XArray3D<xar::fcomplex> aaa);
	// Forward FFT
	void ForwardFFT() { fftwf_execute(aplan); }
	// Inverse FFT
	void InverseFFT() { fftwf_execute(bplan); }
	// Print internal arrays
	void PrintRealArray(const char* message);
	void PrintComplexArray(const char* message);

private:
	// Member variables
	int nx; // x dimension
	int ny; // y dimension
	int nz; // z dimension
	unsigned int uiflag; // FFTW_MEASURE or FFTW_ESTIMATE

	float* pin; // pointer to the real array
	fftwf_complex* pout; // pointer to the complex array

	fftwf_plan aplan; // forward plan
	fftwf_plan bplan; // inverse plan
};


#endif	// FFTWF3DRC_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//