//Header XA_spln3.h
//
//
//	HEADER FILE TITLE:
//
//		Three-dimensional spline interpolations
//
/*!
	\file		XA_spln3.h
	\brief		Three dimensional spline interpolations and related operations
	\par		Description:
		This class implements simple spline interpolations of XArray3D<T> objects and related operations.
*/
#pragma once
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_head3.h"
#include "XA_fft2.h"

namespace xar
{
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//  
//---------------------------------------------------------------------------
//	MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//Class XArray3DSpln
//
//	Three-dimensional spline interpolations
//
/*!
	\brief		Three-dimensional spline interpolations and related operations
	\par		Description:
				This class template defines a 'wrapper' around the XArray3D<T> object
				on which it operates; it contains functions implementing simple
				spline interpolations of the 'wrapped' XArray3D<T> object.
	\remarks	The wrapped member 3D array is declared const, and functions of this class place result in a separate 3D array.
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
*/

	template <class T> class XArray3DSpln
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray3DSpln(const XArray3D<T>& rXAr3D);
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XArray3DSpln(const XArray3DSpln<T>& rCopy) : m_rXArray3D(rCopy.m_rXArray3D) { Initialize(); }
	public:
		//! Destructor 
		~XArray3DSpln() {}
	
	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XArray3DSpln<T>& rCopy);

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XArray3DSpln<T> objects for types T other than float or double
		//! Returns the xar::_eValueType corresponding to T
		static _eValueType GetValuetype(void);
		//! Returns a reference to the non-modifiable 'wrapped' XArray3D<T> object
		const XArray3D<T>& GetBaseObject() const { return m_rXArray3D; }
		//! Returns a reference to the 'wrapped' XArray3D<T> object
		const XArray3D<T>& GetBaseObject() { return m_rXArray3D; }

	// Operations
	public:
		//! Rotates 3D array around a given point with respect to two axes 
		void Rotate(XArray3D<T>& xaResult3D, double angleY, double angleX, double zc = -1.1e-11, double yc = -1.1e-11, double xc = -1.1e-11, T Backgr = T(0)) const;
		//! Rotates 3D array around a given point with respect to two axes 
		void Rotate1(XArray3D<T>& xaResult3D, double angleY, double angleX, double zc = -1.1e-11, double yc = -1.1e-11, double xc = -1.1e-11, T Backgr = T(0)) const;
		//! Calculates multislice propagation of a plane electron wave through the 3D distribution of the scaled electrostatic potential defined by the wrapped object rXAr3D
		void Multislice_eV(XArray2D<std::complex<T> >& xaCamp2D, double angleY, double angleX, double sliceTh, double q2max = -1.0) const;

	// Overridables

	// Implementation
	private:
	// Member variables	
		//! Reference to the 'wrapped' XArray3D<T> object that contains array of data at equidistant x-points 
		const XArray3D<T>& m_rXArray3D;
		//! Array dimensions
		index_t m_nz, m_ny, m_nx;
		//! Wavehead3D data or its default equivalent
		double m_wl, m_zlo, m_zhi, m_ylo, m_yhi, m_xlo, m_xhi;
		//! Auxiliary data
		double m_zst, m_yst, m_xst;

	// Helper functions
		void Initialize();
	};

}  //namespace xar closed


//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//! Returns the xar::_eValueType corresponding to T=char
template<> inline xar::_eValueType xar::XArray3DSpln<char>::GetValuetype() { return eXAChar; }
//! Returns the xar::_eValueType corresponding to T=short
template<> inline xar::_eValueType xar::XArray3DSpln<short>::GetValuetype() { return eXAShort; }
//! Returns the xar::_eValueType corresponding to T=long
template<> inline xar::_eValueType xar::XArray3DSpln<long>::GetValuetype() { return eXALong; }
//! Returns the xar::_eValueType corresponding to T=float
template<> inline xar::_eValueType xar::XArray3DSpln<float>::GetValuetype() { return eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> inline xar::_eValueType xar::XArray3DSpln<double>::GetValuetype() { return eXADouble; }

// NOTE!!! Complex-valued arrays would require specialized templates for some functions of this class

//! Constructor
template <class T> xar::XArray3DSpln<T>::XArray3DSpln(const XArray3D<T>& rXAr3D) : m_rXArray3D(rXAr3D)
{
	GetValuetype();
	Initialize();
}


//! Initialize member variables (has to be called again after m_rXArray3D changes)
template <class T> void xar::XArray3DSpln<T>::Initialize()
{
	m_nz = m_rXArray3D.GetDim1();
	m_ny = m_rXArray3D.GetDim2();
	m_nx = m_rXArray3D.GetDim3();
	
	const IXAHWave3D* pHead = GetIXAHWave3D(m_rXArray3D);
	if (pHead)
	{
		m_wl = pHead->GetWl();
		m_zlo = pHead->GetZlo();
		m_zhi = pHead->GetZhi();
		m_ylo = pHead->GetYlo();
		m_yhi = pHead->GetYhi();
		m_xlo = pHead->GetXlo();
		m_xhi = pHead->GetXhi();
	}
	else
	{
		m_wl = 0.0001;
		m_zlo = 0;
		m_zhi = (double)m_nz;
		m_ylo = 0;
		m_yhi = (double)m_ny;
		m_xlo = 0;
		m_xhi = (double)m_nx;
	}
	m_zst = GetZStep(m_rXArray3D);
	m_yst = GetYStep(m_rXArray3D);
	m_xst = GetXStep(m_rXArray3D);
}


//! Assignment  (declared protected to prohibit copying)
template <class T> void xar::XArray3DSpln<T>::operator=(const XArray3DSpln<T>& rCopy)
{ 
	if (this == &rCopy) 
		return; 
	else
		throw std::logic_error("logic_error in XArray3DSpln<T>::operator= (array reference cannot be changed after construction)");
}

		
//---------------------------------------------------------------------------
//Function XArray3DSpln<T>::Rotate
//
//	Rotates 3D array around two axes
//
/*!
	\brief		Rotates 3D array around two axes
	\param		xaResult3D Resultant rotated 3D array
	\param		angleY rotation angles (in radians) in xz plane (i.e. around y axis)
	\param		angleX rotation angles (in radians) in yz' plane (i.e. around x' axis)
	\param		zc Z-coordinate of the centre of rotation
	\param		yc Y-coordinate of the centre of rotation
	\param		xc X-coordinate of the centre of rotation
	\param		Backgr Value for filling "background" areas around the rotated array
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function rotates 3D array around a specified point by the specified angles	using trilinear extrapolation

	//!!! NOTE it seems that the function Rotate1 below works much better (more accurately)
		
*/	
template <class T> void xar::XArray3DSpln<T>::Rotate(XArray3D<T>& xaResult3D, double angleY, double angleX, double zc, double yc, double xc, T Backgr) const
{
	// check and enforce default coordinates of the centre of rotation
	if (zc == -1.1e-11) zc = 0.5 * (m_zlo + m_zhi);
	if (yc == -1.1e-11) yc = 0.5 * (m_ylo + m_yhi);
	if (xc == -1.1e-11) xc = 0.5 * (m_xlo + m_xhi);

	// we don't expect xaResult3D to be necesserily properly "formatted/dimensioned" before this function is called
	if (xaResult3D.GetDim1() != m_nz || xaResult3D.GetDim2() != m_ny || xaResult3D.GetDim3() != m_nx) xaResult3D.Resize(m_nz, m_ny, m_nx, Backgr);
	else xaResult3D.Fill(Backgr);

	// create new head (note that the new head is created even if the array did not have a had before the rotation)
	IXAHWave3D* pHead = CreateWavehead3D();
	pHead->SetData(m_wl, m_zlo, m_zhi, m_ylo, m_yhi, m_xlo, m_xhi);
	xaResult3D.SetHeadPtr(pHead);

	// calculate the coordinate illumination angle parameters
	double xxx;
	double sinangleY(sin(angleY)), cosangleY(cos(angleY)), sinangleX(sin(angleX)), cosangleX(cos(angleX));
	vector<double> x_sinangleY(m_nx), x_cosangleY(m_nx);
	for (index_t i = 0; i < m_nx; i++)
	{
		xxx = m_xlo + m_xst * i - xc;
		x_sinangleY[i] = (-xxx) * sinangleY;
		x_cosangleY[i] = xxx * cosangleY;
	}
	double yyy;
	vector<double> y_sinangleX(m_ny), y_cosangleX(m_ny);
	for (index_t j = 0; j < m_ny; j++)
	{
		yyy = m_ylo + m_yst * j - yc;
		y_sinangleX[j] = (-yyy) * sinangleX;
		y_cosangleX[j] = yyy * cosangleX;
	}
	double zzz;
	vector<double> z_sinangleY(m_nz), z_cosangleY(m_nz);
	for (index_t n = 0; n < m_nz; n++)
	{
		zzz = m_zlo + m_zst * n - zc;
		z_sinangleY[n] = zzz * sinangleY;
		z_cosangleY[n] = zzz * cosangleY;
	}

	int ii, jj, nn, nx1 = int(m_nx) - 1, ny1 = int(m_ny) - 1, nz1 = int(m_nz) - 1;
	double zz, zz_sinangleX, zz_cosangleX, dK, dx0, dx1, dy0, dy1, dz0, dz1, dz0K, dz1K;
	bool bii, bjj, bnn;

	for (int n = 0; n < m_nz; n++)
	{
		for (int i = 0; i < m_nx; i++)
		{
			//rotation around Y axis
			xxx = xc + x_cosangleY[i] + z_sinangleY[n]; // x coordinate with respect to the rotated 3D sample
			dx1 = (xxx - m_xlo) / m_xst; 
			ii = int(floor(dx1)); if (ii < 0 || ii >= m_nx) continue;
			bii = (ii < nx1);
			dx1 -= ii; dx0 = 1.0 - dx1;
			zz = x_sinangleY[i] + z_cosangleY[n];
			zz_sinangleX = yc + zz * sinangleX;
			zz_cosangleX = zc + zz * cosangleX;
			
			for (int j = 0; j < m_ny; j++)
			{
				// rotation around X' axis
				yyy = y_cosangleX[j] + zz_sinangleX; // y coordinate with respect to the rotated 3D sample
				zzz = y_sinangleX[j] + zz_cosangleX;  // z coordinate with respect to the rotated 3D sample
				dy1 = (yyy - m_ylo) / m_yst; 
				jj = int(floor(dy1)); if (jj < 0 || jj >= m_ny) continue;
				bjj = (jj < ny1);
				dy1 -= jj; dy0 = 1.0 - dy1;
				dz1 = (zzz - m_zlo) / m_zst; 
				nn = int(floor(dz1)); if (nn < 0 || nn >= m_nz) continue;
				bnn = (nn < nz1);
				dz1 -= nn; dz0 = 1.0 - dz1;
				dK = m_rXArray3D[n][j][i];
				dz0K = dz0 * dK;
				dz1K = dz1 * dK;
				xaResult3D[nn][jj][ii] += T(dx0 * dy0 * dz0K);
				if (bii)
				{
					xaResult3D[nn][jj][ii + 1] += T(dx1 * dy0 * dz0K);
					if (bjj)
					{
						xaResult3D[nn][jj + 1][ii + 1] += T(dx1 * dy1 * dz0K);
						if (bnn) xaResult3D[nn + 1][jj + 1][ii + 1] += T(dx1 * dy1 * dz1K);
					}
					if (bnn) xaResult3D[nn + 1][jj][ii + 1] += T(dx1 * dy0 * dz1K);
				}
				if (bjj)
				{
					xaResult3D[nn][jj + 1][ii] += T(dx0 * dy1 * dz0K);
					if (bnn) xaResult3D[nn + 1][jj + 1][ii] += T(dx0 * dy1 * dz1K);
				}
				if (bnn) xaResult3D[nn + 1][jj][ii] += T(dx0 * dy0 * dz1K);
			}
		}
	}
}

//---------------------------------------------------------------------------
//Function XArray3DSpln<T>::Rotate1
//
//	Rotates 3D array around two axes
//
/*!
	\brief		Rotates 3D array around two axes
	\param		xaResult3D Resultant rotated 3D array
	\param		angleY rotation angles (in radians) in xz plane (i.e. around y axis)
	\param		angleX rotation angles (in radians) in yz' plane (i.e. around x' axis)
	\param		zc Z-coordinate of the centre of rotation
	\param		yc Y-coordinate of the centre of rotation
	\param		xc X-coordinate of the centre of rotation
	\param		Backgr Value for filling "background" areas around the rotated array
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function rotates 3D array around a specified point by the specified angles	using trilinear interpolation

	//NOTE: !!! it seems that this function works much better (more accurately) than the preceding function Rotate

*/
template <class T> void xar::XArray3DSpln<T>::Rotate1(XArray3D<T>& xaResult3D, double angleY, double angleX, double zc, double yc, double xc, T Backgr) const
{
	// check and enforce default coordinates of the centre of rotation
	if (zc == -1.1e-11) zc = 0.5 * (m_zlo + m_zhi);
	if (yc == -1.1e-11) yc = 0.5 * (m_ylo + m_yhi);
	if (xc == -1.1e-11) xc = 0.5 * (m_xlo + m_xhi);

	// we don't expect xaResult3D to be necesserily properly "formatted/dimensioned" before this function is called
	if (xaResult3D.GetDim1() != m_nz || xaResult3D.GetDim2() != m_ny || xaResult3D.GetDim3() != m_nx) xaResult3D.Resize(m_nz, m_ny, m_nx, Backgr);
	else xaResult3D.Fill(Backgr);

	// create new head (note that the new head is created even if the array did not have a had before the rotation)
	IXAHWave3D* pHead = CreateWavehead3D();
	pHead->SetData(m_wl, m_zlo, m_zhi, m_ylo, m_yhi, m_xlo, m_xhi);
	xaResult3D.SetHeadPtr(pHead);

	// calculate the coordinate illumination angle parameters
	double xxx;
	double sinangleY(sin(angleY)), cosangleY(cos(angleY)), sinangleX(sin(angleX)), cosangleX(cos(angleX));
	vector<double> x_sinangleY(m_nx), x_cosangleY(m_nx);
	for (index_t i = 0; i < m_nx; i++)
	{
		xxx = m_xlo + m_xst * i - xc;
		x_sinangleY[i] = xxx * sinangleY;
		x_cosangleY[i] = xxx * cosangleY;
	}
	double yyy;
	vector<double> y_sinangleX(m_ny), y_cosangleX(m_ny);
	for (index_t j = 0; j < m_ny; j++)
	{
		yyy = m_ylo + m_yst * j - yc;
		y_sinangleX[j] = yyy * sinangleX;
		y_cosangleX[j] = yyy * cosangleX;
	}
	double zzz;
	vector<double> z_sinangleX(m_nz), z_cosangleX(m_nz);
	for (index_t n = 0; n < m_nz; n++)
	{
		zzz = m_zlo + m_zst * n - zc;
		z_sinangleX[n] = zzz * sinangleX;
		z_cosangleX[n] = zzz * cosangleX;
	}

	int ii, jj, nn, nx2 = int(m_nx) - 2, ny2 = int(m_ny) - 2, nz2 = int(m_nz) - 2;
	double zz, zz_sinangleY, zz_cosangleY, dx0, dx1, dy0, dy1, dz0, dz1;

	for (int n = 0; n < m_nz; n++)
	{
		for (index_t j = 0; j < m_ny; j++)
		{
			// inverse rotation around X' axis
			yyy = yc + y_cosangleX[j] - z_sinangleX[n]; // y coordinate with respect to the rotated 3D sample
			dy1 = (yyy - m_ylo) / m_yst;
			jj = (int)(floor(dy1)); if (jj < 0 || jj > ny2) continue;
			dy1 -= jj; dy0 = 1.0 - dy1;
			zz = y_sinangleX[j] + z_cosangleX[n];
			zz_sinangleY = xc - zz * sinangleY;
			zz_cosangleY = zc + zz * cosangleY;

			for (index_t i = 0; i < m_nx; i++)
			{
				// inverse rotation around Y axis
				xxx = x_cosangleY[i] + zz_sinangleY; // x coordinate with respect to the rotated 3D sample
				zzz = x_sinangleY[i] + zz_cosangleY; // z coordinate with respect to the rotated 3D sample
				dx1 = (xxx - m_xlo) / m_xst;
				ii = (int)floor(dx1); if (ii < 0 || ii > nx2) continue;
				dx1 -= ii; dx0 = 1.0 - dx1;
				dz1 = (zzz - m_zlo) / m_zst;
				nn = (int)floor(dz1); if (nn < 0 || nn > nz2) continue;
				dz1 -= nn; dz0 = 1.0 - dz1;
				xaResult3D[n][j][i] = m_rXArray3D[nn][jj][ii] * T(dx0 * dy0 * dz0) +
				m_rXArray3D[nn][jj][ii + 1] * T(dx1 * dy0 * dz0) +
				m_rXArray3D[nn][jj + 1][ii + 1] * T(dx1 * dy1 * dz0) +
				m_rXArray3D[nn + 1][jj + 1][ii + 1] * T(dx1 * dy1 * dz1) +
				m_rXArray3D[nn + 1][jj][ii + 1] * T(dx1 * dy0 * dz1) +
				m_rXArray3D[nn][jj + 1][ii] * T(dx0 * dy1 * dz0) +
				m_rXArray3D[nn + 1][jj + 1][ii] * T(dx0 * dy1 * dz1) +
				m_rXArray3D[nn + 1][jj][ii] * T(dx0 * dy0 * dz1);
			}
		}
	}
}


//---------------------------------------------------------------------------
//Function XArray3DSpln<T>::Multislice_eV
//
// Calculates multislice propagation of a plane electron wave through the 3D distribution of the scaled electrostatic potential defined by the wrapped object rXAr3D
//
/*!
	\brief		Calculates multislice propagation of a plane electron wave through the 3D distribution of electrostatic potential defined by the wrapped object rXAr3D
	\param		xaCamp2D Resultant transmitted complex wave amplitude in the exit plane
	\param		angleY rotation angles (in radians) in xz plane (i.e. around y axis)
	\param		angleX rotation angles (in radians) in yz' plane (i.e. around x' axis)
	\param      sliceTh slice thickness in the units of the wrapped object
	\param		q2max Defines the optional bandwidth limit (q2max <= 0.0 is interepreted as infinite aperture)
	\param		rXar3D is supposed to contain the 3D distribution of -2*Pi/lambda*delta(x,y,z) = -Pi/(lambda*E)*V(x,y,z)
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function calculates multislice propagation of a plane electron wave through the 3D distribution of the scaled electrostatic potential defined by the wrapped object rXAr3D

*/
template <class T> void xar::XArray3DSpln<T>::Multislice_eV(XArray2D<std::complex<T> >& xaCamp2D, double angleY, double angleX, double sliceTh, double q2max) const
{
	if (sliceTh <= 0 || sliceTh < m_zst) throw std::invalid_argument("invalid slice thickness in XArray3DSpln<T>::Multislice3DSpln (sliceTh is not positive or smaller than z-step in the wrapped 3D array)");

	if (xaCamp2D.GetDim1() != m_ny || xaCamp2D.GetDim2() != m_nx) 
		throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (x or y dimension is different from x or y dimension of the 3D wrapped array with electrostatic potential)");

	index_t i = 2;
	while (i < m_ny) i *= 2;
	if (i != m_ny) throw std::invalid_argument("invalid dim1 in XArray3DSpln<T>::Multislice3DSpln (dim1 is not a power of 2)");
	index_t j = 2;
	while (j < m_nx) j *= 2;
	if (j != m_nx) throw std::invalid_argument("invalid dim2 in XArray3DSpln<T>::Multislice3DSpln (dim2 is not a power of 2)");

	IXAHWave2D* ph2 = GetIXAHWave2D(xaCamp2D);
	if (!ph2) throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (no Wavehead2D)");
	ph2->Validate();

	//if (ph2->GetWl() != m_wl) throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (wavelength is different from that in the wrapped 3D array)");
	double wl = ph2->GetWl(); // we will use this wavelength even if it does not coincide with m_wl
	if (ph2->GetYlo() != m_ylo) throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (m_ylo is different from that in the wrapped 3D array)");
	if (ph2->GetYhi() != m_yhi) throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (m_yhi is different from that in the wrapped 3D array)");
	if (ph2->GetXlo() != m_xlo) throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (m_xlo is different from that in the wrapped 3D array)");
	if (ph2->GetXhi() != m_xhi) throw std::invalid_argument("invalid argument xaCamp2D in XArray3DSpln<T>::Multislice3DSpln (m_xhi is different from that in the wrapped 3D array)");

	// calculate the coordinate illumination angle parameters
	double xxx;
	double sinangleY(sin(angleY)), cosangleY(cos(angleY)), sinangleX(sin(angleX)), cosangleX(cos(angleX));
	double xc = (m_xlo + m_xhi) / 2.0, yc = (m_ylo + m_yhi) / 2.0, zc = (m_zlo + m_zhi) / 2.0;
	vector<double> x_sinangleY(m_nx), x_cosangleY(m_nx);
	for (index_t i = 0; i < m_nx; i++)
	{
		xxx = m_xlo + m_xst * i - xc;
		x_sinangleY[i] = xxx * sinangleY;
		x_cosangleY[i] = xxx * cosangleY;
	}
	double yyy;
	vector<double> y_sinangleX(m_ny), y_cosangleX(m_ny);
	for (index_t j = 0; j < m_ny; j++)
	{
		yyy = m_ylo + m_yst * j - yc;
		y_sinangleX[j] = yyy * sinangleX;
		y_cosangleX[j] = yyy * cosangleX;
	}
	double zzz;
	vector<double> z_sinangleX(m_nz), z_cosangleX(m_nz);
	for (index_t n = 0; n < m_nz; n++)
	{
		zzz = m_zlo + m_zst * n - zc;
		z_sinangleX[n] = zzz * sinangleX;
		z_cosangleX[n] = zzz * cosangleX;
	}

	//bool bEdge;
	int nn0(0), nn1(0), nnstep = int(sliceTh / m_zst + 0.5);
	int ii, jj, nn, nx2 = int(m_nx) - 2, ny2 = int(m_ny) - 2, nz2 = int(m_nz) - 2;
	double zz, zz_sinangleY, zz_cosangleY, dx0, dx1, dy0, dy1, dz0, dz1, sliced2 = sliceTh / 2.0;
	T phinterp;

	xar::XArray2D<T> phi(m_ny, m_nx);
	xar::XArray2DFFT<T> xafft2(xaCamp2D);

	//xafft2.Fresnel(sliced2, q2max);
	while (true) // cycle over the slices
	{
		nn0 = nn1; nn1 += nnstep;
		if (nn1 > m_nz) 
			if (nn0 >= m_nz) break;
			else
			{
				nn1 = int(m_nz);
				sliceTh *= double(nn1 - nn0) / double(nnstep);
			}
				
		// integrate the exp(iV) within the current slice and multiply the complex amplitude by the result
		phi.Fill(T(0));
		for (int n = nn0; n < nn1; n++)
		{
			//if (n == nn0 || n == nn1 - 1) bEdge = true; else bEdge = false;
			for (index_t j = 0; j < m_ny; j++)
			{
				// inverse rotation around X' axis
				yyy = yc  + y_cosangleX[j] - z_sinangleX[n]; // y coordinate with respect to the rotated 3D sample
				dy1 = (yyy - m_ylo) / m_yst;
				jj = (int)(floor(dy1)); if (jj < 0 || jj > ny2) continue;
				dy1 -= jj; dy0 = 1.0 - dy1;
				zz = y_sinangleX[j] + z_cosangleX[n];
				zz_sinangleY = xc - zz * sinangleY;
				zz_cosangleY = zc + zz * cosangleY;

				for (index_t i = 0; i < m_nx; i++)
				{
					// inverse rotation around Y axis
					xxx = x_cosangleY[i] + zz_sinangleY; // x coordinate with respect to the rotated 3D sample
					zzz = x_sinangleY[i] + zz_cosangleY; // z coordinate with respect to the rotated 3D sample
					dx1 = (xxx - m_xlo) / m_xst;
					ii = (int)floor(dx1); if (ii < 0 || ii > nx2) continue;
					dx1 -= ii; dx0 = 1.0 - dx1;
					dz1 = (zzz - m_zlo) / m_zst;
					nn = (int)floor(dz1); if (nn < 0 || nn > nz2) continue;
					dz1 -= nn; dz0 = 1.0 - dz1;
					phinterp = m_rXArray3D[nn][jj][ii] * T(dx0 * dy0 * dz0) +
						m_rXArray3D[nn][jj][ii + 1] * T(dx1 * dy0 * dz0) +
						m_rXArray3D[nn][jj + 1][ii + 1] * T(dx1 * dy1 * dz0) +
						m_rXArray3D[nn + 1][jj + 1][ii + 1] * T(dx1 * dy1 * dz1) +
						m_rXArray3D[nn + 1][jj][ii + 1] * T(dx1 * dy0 * dz1) +
						m_rXArray3D[nn][jj + 1][ii] * T(dx0 * dy1 * dz0) +
						m_rXArray3D[nn + 1][jj + 1][ii] * T(dx0 * dy1 * dz1) +
						m_rXArray3D[nn + 1][jj][ii] * T(dx0 * dy0 * dz1);
					//if (bEdge) phi[j][i] += T(0.5) * phinterp; else phi[j][i] += phinterp;
					phi[j][i] += phinterp;
				}
			}
		}
		// multiply the complex amplitude by exp(-i*2*Pi/lambda*slice_projected_phase)
		MultiplyExpiFi(xaCamp2D, phi);
		// Fresnel propagate slice thickness forward
		xafft2.Fresnel(sliceTh, q2max);
	}
	xafft2.Fresnel(m_zst, q2max); // we want to propagate to zhi, rather than zhi - zst, which is what zlo + zst * (m_nz - 1) is equal to
}


//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//

// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
//template class xar::XArray3DSpln<char>;
//template class xar::XArray3DSpln<short>;
//template class xar::XArray3DSpln<long>;
template class xar::XArray3DSpln<float>;
template class xar::XArray3DSpln<double>;

//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//