//Header XArray.h
//
//
//	HEADER FILE TITLE:
//
//		Common base class template for XArray classes
//
//	COPYRIGHT:
//
//		Commonwealth Scientific and Industrial Research Organisation, 2001-2007
//					All Rights Reserved
//
//
/*!
	\file		XArray.h
	\brief		Common base class template for XArray classes
	\par		Description:
		This is a common base class template for XArrayND<T> classes with different dimensions N.
		Also defined here are some basic XArray<complex<T> > - related non-member functions.

*/
#if !defined XARRAY_H
#define XARRAY_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include <cmath>
#include <limits>
#include <algorithm>

#include "XArrayB.h" // XArrayBase class
#include "IXAHead.h" // base IXAHead class

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
//Class XArray<T>
//
//	Common base class template for XArray classes
//
/*!
	\brief		Common base class template for XArray classes
	\par		Description:
				This class is designed as an abstract implementation class which augments XArrayBase<T> base class
				functionality by providing additional useful dimensionality-neutral operations for XArray classes.
				Many useful functions are inherited from the base class, XArrayBase<T>.
				XArray class is also responsible for (has ownership of) the associated head.
	\remarks	This class is supposed to be 'abstract' and cannot be instantiated (only its
				dimension-specific descendants, like XArray1D<T>, can be instantiated)
	\remarks	Resizing functions AcceptMemBuffer, ReleaseMemBuffer, Resize, Truncate, Swap, etc. are defined in derived classes
				as they are dimension-specific
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
	\warning	The absence of appropriate specializations of the function GetValuetype()
				will prevent instantiation of XArray<T> object for new types T
*/
	template <class T> class XArray : public XArrayBase<T>
	{
	// Enumerators
	// Structures
	// Constructors
		//! Constructors are protected to prevent instantiation of objects of this class
		// GetValuetype() is called in constructors only to prevent the compillation 
		//		of XArray<T>-derived objects with unsupported T types
	protected:
		//! Default constructor
		XArray() : m_pHead(0) {  GetValuetype(); }
		//! Constructor with a predefined size
		explicit XArray(index_t NumPoints, T tVal = T()) : XArrayBase<T>(NumPoints, tVal), m_pHead(0) { GetValuetype(); }
		 //! Promotion from XArrayBase
		explicit XArray(const XArrayBase<T>& rXArrayBase) : XArrayBase<T>(rXArrayBase), m_pHead(0) { GetValuetype(); }
		//! Construction from a raw memory buffer
		XArray(T* ptBufBegin, T* ptBufEnd) : XArrayBase<T>(ptBufBegin, ptBufEnd), m_pHead(0) { GetValuetype(); } 
		 //! Copy constructor
		XArray(const XArray<T>& rXArray) : XArrayBase<T>(rXArray), m_pHead(rXArray.m_pHead ? rXArray.m_pHead->Clone() : 0) {}
		 //! Destructor, note that XArray has ownership of its head
		~XArray() { delete m_pHead; }

	// Operators
	public:
		//! Makes this array a (deep) copy of the rXArray
		void operator=(const XArray<T>& rXArray);
		//! Makes this array a (deep) copy of the rXArrayBase and sets the head pointer to 0
		void operator=(const XArrayBase<T>& rXArrayBase);
		//! Adds a given value to each array element
		void operator+=(T tVal);
		//! Subtracts a given value from each array element
		void operator-=(T tVal);
		//! Multiplies each array element by a given value
		void operator*=(T tVal);
		//! Divides each array element by a given value  (checks for division by zero)
		void operator/=(T tVal);
		//! Raises each array element to a given power
		void operator^=(T tVal);
		// The following operators will throw an exception if SameHead returns false
		// (this is why m_pHead must be a member of XArray, and not of derived classes)
		//
		//! Performs elementwise addition of the two arrays
		void operator+=(const XArray<T>& rXArray);
		//! Performs elementwise subtraction of the two arrays
		void operator-=(const XArray<T>& rXArray);
		//! Performs elementwise multiplication of the two arrays
		void operator*=(const XArray<T>& rXArray);
		//! Performs elementwise division of the two arrays (checks for division by zero)
		void operator/=(const XArray<T>& rXArray);

	// Attributes
	public:
		//! Returns the xar::_eValueType corresponding to T
		// The absense of relevant specialization will prevent instantiation
		// of XArray<T> classes for unsupported T-types (this is the intended behavior)
		static _eValueType GetValuetype(void);
		//! Returns constant pointer to the associated head (zero by default)
		const IXAHead* GetHeadPtr(void) const { return m_pHead; }
		//! Returns a modifiable pointer to the associated head (zero by default)
		IXAHead* GetHeadPtr(void) { return m_pHead; }
		//! Sets a new head pointer (deletes the old head!)
		void SetHeadPtr(IXAHead* pHead) { if (pHead) pHead->Validate(); delete m_pHead; m_pHead = pHead; }
		//! Calculates various standard norms (metrics) of the XArray data
		double Norm(_eNormID eMode) const; 
		//! Calculates chi-square distance from another XArray image
		double Chi2(const XArray<T>& rXArray, double dblRelStDevAver, bool bUsePoissonStat) const;

	// Operations
	public:
		//! If m_pHead==0, then returns rXArray.m_pHead==0, else calls an overridable function m_pHead->Equivalent(rXArray.m_pHead)
		bool SameHead(const XArray<T>& rXArray) const;
		//! Provides slow (but checked and polymorphic) read access to an element
		double GetAt(index_t index) const;
		//! Provides slow (but checked and polymorphic) write access to an element
		void SetAt(index_t index, double dblVal);
		//! Provides slow (but checked and polymorphic) read access to an element
		dcomplex GetCmplAt(index_t index) const;
		//! Provides slow (but checked and polymorphic) write access to an element
		void SetCmplAt(index_t index, dcomplex cxdVal);
		//! Fills the array with a given fixed value
		void Fill(T tVal) { std::fill(XArrayBase<T>::begin(), XArrayBase<T>::end(), tVal); }
		//! Replaces each array member with its modulus
		void Abs(void);
		//! Replaces each array member with its square modulus
		void Abs2(void);
		//! Replaces each array member with its exponent
		void Exp(void);
		//! Replaces each array member with its logarithm
		void Log(void);
		//! Replaces each array member with its sine
		void Sin(void); 
		//! Replaces each array member with its cosine
		void Cos(void);
		//! Replaces each array member with its tangent
		void Tan(void);
		//! Replaces each array member with its arcsine
		void Asin(void);
		//! Replaces each array member with its arccosine
		void Acos(void);
		//! Replaces each array member with its arctangent
		void Atan(void);

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
		 //! This pointer is an interface to a loosely coupled head (it is 0 by default)
		IXAHead* m_pHead;

	// Member functions

	};

}  //namespace xar closed

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
namespace xar
{
	//
	//NOTE: the absense of relevant specialization will prevent instantiation
	//		of XArray<T> classes for unsupported T-types (this is the intended behavior)
	//
	inline _eValueType XArray<char>::GetValuetype() { return eXAChar; }	
	inline _eValueType XArray<short>::GetValuetype() { return eXAShort; }	
	inline _eValueType XArray<long>::GetValuetype() { return eXALong; }	
	inline _eValueType XArray<float>::GetValuetype() { return eXAFloat; }
	inline _eValueType XArray<double>::GetValuetype() { return eXADouble; }
	inline _eValueType XArray<fcomplex>::GetValuetype() { return eXAFComplex; }
	inline _eValueType XArray<dcomplex>::GetValuetype() { return eXADComplex; }


	template <class T> bool XArray<T>::SameHead(const XArray<T>& xa) const
	{
		if (m_pHead)
		{
			return m_pHead->Equivalent(xa.m_pHead);
		}
		else
		{
			return xa.m_pHead == 0; 
		}
	}


	//***** XArray member operators

	template <class T> void XArray<T>::operator=(const XArray<T>& xa)
	{
		if (this == &xa) return;
		if (xa.m_pHead) xa.m_pHead->Validate();
		XArrayBase<T>::operator=(xa);
		delete m_pHead;
		m_pHead = xa.m_pHead ? xa.m_pHead->Clone() : 0;
	}


	template <class T> void XArray<T>::operator=(const XArrayBase<T>& va)
	{
		XArrayBase<T>::operator=(va);
		delete m_pHead; m_pHead = 0;
	}


	template <class T> void XArray<T>::operator+=(const XArray<T>& xa)
	{
		if (xa.size() != XArrayBase<T>::size()) throw std::range_error("range_error lhs.size != rhs.size in XArray<T>::operator+=");
		if (!SameHead(xa))
			throw std::invalid_argument("invalid_argument 'rhs' in XArray<T>::operator+= (different head)"); 	
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) += xa[i];
	}
	
	
	template <class T> void XArray<T>::operator-=(const XArray<T>& xa)
	{
		if (xa.size() != XArrayBase<T>::size()) throw std::range_error("range_error lhs.size != rhs.size in XArray<T>::operator-=");
		if (!SameHead(xa))
			throw std::invalid_argument("invalid_argument 'rhs' in XArray<T>::operator-= (different head)"); 	
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) -= xa[i];
	}
	
	
	template <class T> void XArray<T>::operator*=(const XArray<T>& xa)
	{
		if (xa.size() != XArrayBase<T>::size()) throw std::range_error("range_error lhs.size != rhs.size in XArray<T>::operator*=");
		if (!SameHead(xa))
			throw std::invalid_argument("invalid_argument 'rhs' in XArray<T>::operator*= (different head)"); 	
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) *= xa[i];
	}

	
	template <class T> void XArray<T>::operator/=(const XArray<T>& xa)
	// This operator can leave *this in an incorrect state, but the alternatives are too expensive
	{
		if (xa.size() != XArrayBase<T>::size()) throw std::range_error("range_error lhs.size != rhs.size in XArray<T>::operator/=");
		if (!SameHead(xa))
			throw std::invalid_argument("invalid_argument 'rhs' in XArray<T>::operator/= (different head)"); 	
		for (index_t i=0; i<XArrayBase<T>::size(); i++) 
		{
			if (xa[i] == T(0)) throw std::invalid_argument("invalid_argument 'rhs' in XArray<T>::operator/= (division by zero)"); 	
			XArrayBase<T>::operator[](i) /= xa[i];
		}
	}

	
	template <class T> void XArray<T>::operator+=(T tVal)
	{
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) += tVal;
	}
	

	template <class T> void XArray<T>::operator-=(T tVal)
	{
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) -= tVal;
	}
	

	template <class T> void XArray<T>::operator*=(T tVal)
	{
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) *= tVal;
	}


	template <class T> void XArray<T>::operator/=(T tVal)
	{
		if (tVal == T(0)) throw std::invalid_argument("invalid_argument 'rhs' in XArray<T>::operator/= (division by zero)"); 	
		else 
		{
			T tAval = T(1) / tVal;
			for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) *= tAval;
		}
	}
	

	template <class T> void XArray<T>::operator^=(T tVal)
	// This operator is specialized separately for complex T
	{
		double dblVal(tVal);
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = T(pow((double)XArrayBase<T>::operator[](i), dblVal));
	}
	

	//***** member functions

	template <class T> inline void XArray<T>::Abs() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = T(fabs(XArrayBase<T>::operator[](i)));
	}


	template <class T> inline void XArray<T>::Abs2() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = XArrayBase<T>::operator[](i) * XArrayBase<T>::operator[](i);
	}


	template <class T> inline void XArray<T>::Exp() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)exp((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Log() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)log((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Sin() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)sin((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Cos() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)cos((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Tan() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)tan((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Asin() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)asin((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Acos() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)acos((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline void XArray<T>::Atan() 
	// This function is specialized separately for complex T
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++) XArrayBase<T>::operator[](i) = (T)atan((double)XArrayBase<T>::operator[](i));
	}


	template <class T> inline double XArray<T>::GetAt(index_t index) const 
	// This function is specialized separately for complex T
	{ 
		return double(XArrayBase<T>::at(index));
	}
	

	template <class T> inline void XArray<T>::SetAt(index_t index, double val)
	// This function is specialized separately for complex T
	{ 
		XArrayBase<T>::at(index) = T(val);
	}
	

	template <class T> inline dcomplex XArray<T>::GetCmplAt(index_t index) const
	// This function is specialized separately for complex T
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<T>::GetCmplAt (must be complex)"); 	
	}
	

	template <class T> inline void XArray<T>::SetCmplAt(index_t index, dcomplex val)
	// This function is specialized separately for complex T
	{
		throw std::invalid_argument("invalid_argument '*this' in XArray<T>::SetCmplAt (must be complex)"); 
	}

	//---------------------------------------------------------------------------
	//Function XArray<T>::Norm
	//
	//	Calculates various standard norms (metrics) of the XArray data
	//
	/*!
		\brief		Calculates various standard norms (metrics) of the XArray data
		\param		eMode	determines the norm (metric) to be calculated
		\exception	std::invalid_argument is thrown if the eMode is unknown or if the requested norm cannot be evalueated for this XArray
		\return		\a anorm the calculated value of the norm
		\par		Description:
			This function can calculate multiple different norms (metrics) of this XArray object (e.g. maximum and minimum 
			value of the elements, RMS metric, etc). The available norms are enumerated by the _eNormID enumerator defined
			in the XA_ini.h module
		\par		Example:
\verbatim
double dblRMS_metric = myXArray2D.Norm(eNormL2); 
\endverbatim
	*/
	template <class T> double XArray<T>::Norm(_eNormID eMode) const
	{
		index_t np = XArrayBase<T>::size();
		index_t i;
		double anorm=0.0, abra=0.0;

		switch (eMode)
		{
		case eNormTypeMin: // T-type min
			if (std::numeric_limits<T>::is_integer) abra = std::numeric_limits<T>::min();
			else abra = -std::numeric_limits<T>::max();
			break;
		case eNormTypeMax: // T-type max
			abra = std::numeric_limits<T>::max();
			break;
		case eNormIndexOfMin: // index-of-min
			anorm = XArrayBase<T>::front(); abra = 0;
			for (i=0; i<np; i++) 
				if (!(XArrayBase<T>::operator[](i) >= anorm)) 
					if (!(XArrayBase<T>::operator[](i) < anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (bad data in array)"); 
					else { anorm = double(XArrayBase<T>::operator[](i)); abra = double(i); }
			anorm = abra;
			break;
		case eNormIndexOfMax: // index-of-max
			anorm = XArrayBase<T>::front(); abra = 0;
			for (i=0; i<np; i++) 
				if (!(XArrayBase<T>::operator[](i) <= anorm)) 
					if (!(XArrayBase<T>::operator[](i) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (bad data in array)"); 
					else { anorm = double(XArrayBase<T>::operator[](i)); abra = double(i); }
			anorm = abra;
			break;
		case eNormYMin: // y-min
			throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (no header)"); 
			break;
		case eNormYMax: // y-max
			throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (no header)"); 
			break;
		case eNormXMin: // x-min
			throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (no header)"); 
			break;
		case eNormXMax: // x-max
			throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (no header)"); 
			break;
		case eNormMin: // min
			anorm = XArrayBase<T>::front();
			for (i=1; i<np; i++) 
				if (!(XArrayBase<T>::operator[](i) >= anorm))
					if (!(XArrayBase<T>::operator[](i) < anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (bad data in array)"); 
					else anorm = XArrayBase<T>::operator[](i);
			break;
		case eNormMax: // max
			anorm = XArrayBase<T>::front();
			for (i=1; i<np; i++) 
				if (!(XArrayBase<T>::operator[](i) <= anorm))
					if (!(XArrayBase<T>::operator[](i) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (bad data in array)"); 
					else anorm = XArrayBase<T>::operator[](i);
			break;
		case eNormC0: // C0-Norm
			anorm = fabs(XArrayBase<T>::front());
			for (i=0; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) <= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormL1: // l1-Norm (not normalized)
			for (i=0; i<np; i++) anorm += fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormL2: // l2-Norm (not normalized)
			for (i=0; i<np; i++) anorm += (double)XArrayBase<T>::operator[](i) * (double)XArrayBase<T>::operator[](i);
			anorm = sqrt(anorm);
			break;
		case eNormAver: // average
			for (i=0; i<np; i++) anorm += XArrayBase<T>::operator[](i);
			anorm /= np;
			break;
		case eNormStdDev: // st.deviation
			for (i=0; i<np; i++) abra += XArrayBase<T>::operator[](i);
			abra /= np;
			for (i=0; i<np; i++) anorm += (XArrayBase<T>::operator[](i) - abra) * (XArrayBase<T>::operator[](i) - abra);
			anorm = sqrt(anorm / np);
			break;
		case eNormL2N: // l2 (normalized)
			for (i=0; i<np; i++) anorm += (double)XArrayBase<T>::operator[](i) * (double)XArrayBase<T>::operator[](i);
			anorm = sqrt(anorm / np);
			break;
		case eNormL1N: // l1 (normalized)
			for (i=0; i<np; i++) anorm += fabs(XArrayBase<T>::operator[](i));
			anorm /= np;
			break;
		case eNormEntropy: // entropy
			{
				double alog2 = 1.0 / log(2.0);
				double totI = 0;
				for (i=0; i<np; i++) 
				{
					if (XArrayBase<T>::operator[](i) < 0) throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Norm(eNormEntropy) (negative values)"); 
					totI += XArrayBase<T>::operator[](i);
				}
				if (totI != 0)
				{
					double mlogtotI = -log(totI) * alog2;
					for (i=0; i<np; i++) 
					{
						if (XArrayBase<T>::operator[](i) != 0) anorm += XArrayBase<T>::operator[](i) * (log((double)XArrayBase<T>::operator[](i)) * alog2 + mlogtotI);
					}
					anorm /= -totI;
				}
			}
			break;
		default:
			throw std::invalid_argument("invalid_argument 'eMode' in XArray<T>::Norm"); 
		}
		return anorm;
	}


	//---------------------------------------------------------------------------
	//Function XArray<T>::Chi2
	//
	//	Calculates chi-square distance from another XArray image
	//
	/*!
		\brief		Calculates chi-square distance from another XArray image
		\param		rXArray another image to be compared to this one
		\param		dblRelStDevAver	(sigma) relative standard deviation corresponding to average intensity
		\param		bUsePoissonStat if true, Poisson (photon counting) statistics is used, else Gaussian additive noise implied
		\exception	std::invalid_argument is thrown if 
		\return		\a chi2 the calculated value of the chi-square distance
		\par		Description:
			This function calculates the chi-square distance between two images
			according to the formula, (1/sigma^2) Sum{ (I1[i]-I0[i])^2 / <I0> / I0[i] }.
			This formula corresponds to the Poisson noise (almost). If bUsePoissonStat=false,
			a modified formula corresponding to additive noise with constant 'sigma' is used,
			(1/sigma^2/<I0>^2) Sum{ (I1[i]-I0[i])^2 }.
		\par		Example:
\verbatim
double dblChi2 = myXArray2D.Chi2(otherXArray); 
\endverbatim
	*/
	template <class T> double XArray<T>::Chi2(const XArray<T>& rXArray, double dblRelStDevAver, bool bUsePoissonStat) const
	{
		if (rXArray.size() != XArrayBase<T>::size()) 
			throw std::invalid_argument("invalid_argument 'rXArray' in XArray<T>::Chi2 (different size)"); 				
		if (!SameHead(rXArray))
			throw std::invalid_argument("invalid_argument 'rXArray' in XArray<T>::Chi2 (different head)"); 	
		if (dblRelStDevAver <= 0)
			throw std::invalid_argument("invalid_argument 'dblRelStDevAver' in XArray<T>::Chi2 (must be positive)"); 	
		double chi2 = 0.0, temp;
		double dblAver = Norm(eNormAver);
		if (bUsePoissonStat) // Poisson photon statistics
		{
			if (Norm(eNormMin) <= 0)
				throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Chi2 (not all values are positive)"); 	
			for (index_t i=0; i<XArrayBase<T>::size(); i++) 
			{
				temp = XArrayBase<T>::operator[](i) - rXArray[i];
				chi2 += temp * temp / XArrayBase<T>::operator[](i);
			}
			chi2 /= dblRelStDevAver * dblRelStDevAver * dblAver;
		}
		else  // additive Gaussian noise with constant sigma
		{
			if (dblAver <= 0)
				throw std::invalid_argument("invalid_argument '*this' in XArray<T>::Chi2 (average value is not positive)"); 	
			for (index_t i=0; i<XArrayBase<T>::size(); i++) 
			{
				temp = XArrayBase<T>::operator[](i) - rXArray[i];
				chi2 += temp * temp;
			}
			chi2 /= dblRelStDevAver * dblRelStDevAver * dblAver * dblAver;
		}
		return chi2;
	}


	//
	//***** complex-specific  specializations
	// It appears that these functions have to be 'inline' to be considered by the compiler
	//

	inline double XArray<fcomplex>::GetAt(index_t index) const
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::GetAt (this cannot be complex)");
	}


	inline void XArray<fcomplex>::SetAt(index_t index, double val)
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::SetAt (this cannot be complex)");
	}


	inline double XArray<dcomplex>::GetAt(index_t index) const 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::GetAt (this cannot be complex)");
	}


	inline void XArray<dcomplex>::SetAt(index_t index, double val)
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::SetAt (this cannot be complex)");
	}


	inline dcomplex XArray<fcomplex>::GetCmplAt(index_t index) const
	{
		return dcomplex(at(index));
	}


	inline void XArray<fcomplex>::SetCmplAt(index_t index, dcomplex val)
	{
		at(index) = fcomplex(val); 
	}


	inline dcomplex XArray<dcomplex>::GetCmplAt(index_t index) const
	{
		return at(index); 
	}


	inline void XArray<dcomplex>::SetCmplAt(index_t index, dcomplex val)
	{
		at(index) = val; 
	}

	//
	// The following specializations are necessary because the standard global math functions 
	// do not work for std::complex, while std::functions do not work for non-complex types
	//

	inline void XArray<fcomplex>::operator^=(fcomplex cxfVal)
	{
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::pow(XArrayBase<T>::operator[](i), cxfVal);
	}


	inline void XArray<dcomplex>::operator^=(dcomplex cxdVal)
	{
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::pow(XArrayBase<T>::operator[](i), cxdVal);
	}


	inline void XArray<fcomplex>::Abs() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Abs (member function undefined)"); 	
	}
	

	inline void XArray<fcomplex>::Abs2() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Abs2 (member function undefined)"); 	
	}


	inline void XArray<fcomplex>::Log() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::log(XArrayBase<T>::operator[](i));
	}


	inline void XArray<fcomplex>::Exp() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::exp(XArrayBase<T>::operator[](i));
	}


	inline void XArray<fcomplex>::Sin() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::sin(XArrayBase<T>::operator[](i));
	}


	inline void XArray<fcomplex>::Cos() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::cos(XArrayBase<T>::operator[](i));
	}


	inline void XArray<fcomplex>::Tan() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Tan (member function undefined)"); 		
	}


	inline void XArray<fcomplex>::Asin() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Asin (member function undefined)"); 	
	}


	inline void XArray<fcomplex>::Acos() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Acos (member function undefined)"); 	
	}


	inline void XArray<fcomplex>::Atan() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Atan (member function undefined)"); 	
	}


	inline void XArray<dcomplex>::Abs() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Abs (member function undefined)"); 	
	}
	

	inline void XArray<dcomplex>::Abs2() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Abs2 (member function undefined)"); 	
	}


	inline void XArray<dcomplex>::Log() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::log(XArrayBase<T>::operator[](i));
	}


	inline void XArray<dcomplex>::Exp() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::exp(XArrayBase<T>::operator[](i));
	}


	inline void XArray<dcomplex>::Sin() 
	{
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::sin(XArrayBase<T>::operator[](i));
	}


	inline void XArray<dcomplex>::Cos() 
	{ 
		for (index_t i=0; i<XArrayBase<T>::size(); i++)  XArrayBase<T>::operator[](i) = std::cos(XArrayBase<T>::operator[](i));
	}


	inline void XArray<dcomplex>::Tan() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Tan (member function undefined)"); 		
	}


	inline void XArray<dcomplex>::Asin() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Asin (member function undefined)"); 	
	}


	inline void XArray<dcomplex>::Acos() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Acos (member function undefined)"); 	
	}


	inline void XArray<dcomplex>::Atan() 
	{ 
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Atan (member function undefined)"); 	
	}


	inline double XArray<fcomplex>::Norm(_eNormID eMode) const
	// NOTE: this function must be defined inline in the Array1D.h to be considered by the compiler alongside 
	// with the default template version of the Norm function
	{
		index_t  np = XArrayBase<T>::size();
		index_t  i;
		double anorm=0.0, abra=0.0;

		switch (eMode)
		{
		case eNormYMin: //y-min
			throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (no header)"); 
			break;
		case eNormYMax: //y-max
			throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (no header)"); 
			break;
		case eNormXMin: //x-min
			throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (no header)"); 
			break;
		case eNormXMax: //x-max
			throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (no header)"); 
			break;
		case eNormMin: //min
			anorm = fabs(XArrayBase<T>::front());
			for (i=1; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) >= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) < anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormMax: //max
			anorm = fabs(XArrayBase<T>::front());
			for (i=1; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) <= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormC0: //C0-Norm
			anorm = fabs(XArrayBase<T>::front());
			for (i=0; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) <= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormL1: //l1-Norm (not normalized)
			for (i=0; i<np; i++) anorm += fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormL2: //l2-Norm (not normalized)
			for (i=0; i<np; i++) anorm += XArrayBase<T>::operator[](i).real() * XArrayBase<T>::operator[](i).real() +  XArrayBase<T>::operator[](i).imag() * XArrayBase<T>::operator[](i).imag();
			anorm = sqrt(anorm);
			break;
		case eNormL2N: //l2 (normalized)
			for (i=0; i<np; i++) anorm += XArrayBase<T>::operator[](i).real() * XArrayBase<T>::operator[](i).real() +  XArrayBase<T>::operator[](i).imag() * XArrayBase<T>::operator[](i).imag();
			anorm = sqrt(anorm / np);
			break;
		case eNormL1N: //l1 (normalized)
			for (i=0; i<np; i++) anorm += fabs(XArrayBase<T>::operator[](i));
			anorm /= np;
			break;
		default:
			throw std::invalid_argument("invalid_argument 'eMode' in XArray<fcomplex>::Norm"); 
		}
		return anorm;
	}


	inline double XArray<dcomplex>::Norm(_eNormID eMode) const
	// NOTE: this function must be defined inline in the Array1D.h to be considered by the compiler alongside 
	// with the default template version of the Norm function
	{
		index_t  np = XArrayBase<T>::size();
		index_t i;
		double anorm = 0.0, abra = 0.0;

		switch (eMode)
		{
		case eNormYMin: //y-min
			throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (no header)"); 
			break;
		case eNormYMax: //y-max
			throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (no header)"); 
			break;
		case eNormXMin: //x-min
			throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (no header)"); 
			break;
		case eNormXMax: //x-max
			throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (no header)"); 
			break;
		case eNormMin: //min
			anorm = fabs(XArrayBase<T>::front());
			for (i=1; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) >= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) < anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormMax: //max
			anorm = fabs(XArrayBase<T>::front());
			for (i=1; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) <= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormC0: //C0-Norm
			anorm = fabs(XArrayBase<T>::front());
			for (i=0; i<np; i++) 
				if (!(fabs(XArrayBase<T>::operator[](i)) <= anorm)) 
					if (!(fabs(XArrayBase<T>::operator[](i)) > anorm)) throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Norm (bad data in array)"); 
					else anorm = fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormL1: //l1-Norm (not normalized)
			for (i=0; i<np; i++) anorm += fabs(XArrayBase<T>::operator[](i));
			break;
		case eNormL2: //l2-Norm (not normalized)
			for (i=0; i<np; i++) anorm += XArrayBase<T>::operator[](i).real() * XArrayBase<T>::operator[](i).real() +  XArrayBase<T>::operator[](i).imag() * XArrayBase<T>::operator[](i).imag();
			anorm = sqrt(anorm);
			break;
		case eNormL2N: //l2 (normalized)
			for (i=0; i<np; i++) anorm += XArrayBase<T>::operator[](i).real() * XArrayBase<T>::operator[](i).real() +  XArrayBase<T>::operator[](i).imag() * XArrayBase<T>::operator[](i).imag();
			anorm = sqrt(anorm / np);
			break;
		case eNormL1N: //l1 (normalized)
			for (i=0; i<np; i++) anorm += fabs(XArrayBase<T>::operator[](i));
			anorm /= np;
			break;
		default:
			throw std::invalid_argument("invalid_argument 'eMode' in XArray<dcomplex>::Norm"); 
		}
		return anorm;
	}

	inline double XArray<fcomplex>::Chi2(const XArray<fcomplex>& rXArray, double dblRelStDevAver, bool bUsePoissonStat) const
	{
		throw std::invalid_argument("invalid_argument '*this' in XArray<fcomplex>::Chi2 (member function undefined)"); 	
	}

	inline double XArray<dcomplex>::Chi2(const XArray<dcomplex>& rXArray, double dblRelStDevAver, bool bUsePoissonStat) const
	{
		throw std::invalid_argument("invalid_argument '*this' in XArray<dcomplex>::Chi2 (member function undefined)"); 	
	}

}
//---------------------------------------------------------------------------
//	RELATED IN-LINE NON-MEMBER DEFINITIONS
//
namespace xar
{

	//---------------------------------------------------------------------------
	//Function XArray<T>::MakeComplex
	//
	//	Makes a complex XArray object C = A + ib or C = A * exp(ib) from a real XArray object A and a scalar b
	//
	/*!
		\brief		Makes a complex XArray object C = A + ib or C = A * exp(ib) from a real XArray object A and a scalar b
		\param		A	Real XArray object representing real part or modulus of a complex object to be constructed
		\param		b	Real scalar value representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray object constructed by this function
		\param		bMakePolar if \b true, then C = A * exp(ib) is constructed, else C = A + ib is constructed
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions called from 
					inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray object C = A + ib or C = A * exp(ib) from a real XArray object and a scalar.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter A
		\par		Example:
\verbatim
XArray1D<float> A(20, 1.0);
XArray1D<fcomplex> C; 
MakeComplex(A, 0.0f, C, false);
\endverbatim
	*/
	template <class T> void MakeComplex(const XArray<T>& A, T b, XArray< std::complex<T> >& C, bool bMakePolar)
	{
		if (A.GetHeadPtr()) A.GetHeadPtr()->Validate();
		C.resize(A.XArrayBase<T>::size());
		if (bMakePolar) for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) C[i] = std::polar<T>(A[i], b);
		else for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) C[i] = std::complex<T>(A[i], b);
		C.SetHeadPtr(A.GetHeadPtr() ? A.GetHeadPtr()->Clone() : 0);
	}


	//---------------------------------------------------------------------------
	//Function XArray<T>::MakeComplex
	//
	//	Makes a complex XArray object C = a + iB or C = a * exp(iB) from a real XArray object B and a scalar a
	//
	/*!
		\brief		Makes a complex XArray object C = a + iB or C = a * exp(iB) from a real XArray object B and a scalar a
		\param		a	Real scalar value representing real part or modulus of a complex object to be constructed
		\param		B	Real XArray object representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray object constructed by this function
		\param		bMakePolar if \b true, then C = a * exp(iB) is constructed, else C = a + iB is constructed
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions called from 
					inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray object C = a + iB or C = a * exp(iB) from a real XArray object and a scalar.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter B
		\par		Example:
\verbatim
XArray1D<float> B(20, 1.0);
XArray1D<fcomplex> C;
MakeComplex(1.0f, B, C, false);
\endverbatim
	*/
	template <class T> void MakeComplex(T a, const XArray<T>& B, XArray< std::complex<T> >& C, bool bMakePolar)
	{
		if (B.GetHeadPtr()) B.GetHeadPtr()->Validate();
		C.resize(B.XArrayBase<T>::size());
		if (bMakePolar) for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) C[i] = std::polar<T>(a, B[i]);
		else for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) C[i] = std::complex<T>(a, B[i]);
		C.SetHeadPtr(B.GetHeadPtr() ? B.GetHeadPtr()->Clone() : 0);
	}
	

	//---------------------------------------------------------------------------
	//Function XArray<T>::MakeComplex
	//
	//	Makes a complex XArray object C = A + iB or C = A * exp(iB) from 2 real XArray objects, A and B
	//
	/*!
		\brief		Makes a complex XArray object C = A + iB or C = A * exp(iB) from 2 real XArray objects, A and B
		\param		A	Real XArray object representing real part or modulus of a complex object to be constructed
		\param		B	Real XArray object representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray object constructed by this function
		\param		bMakePolar if \b true, then C = A * exp(iB) is constructed, else C = A + iB is constructed
		\exception	std::invalid_argument is thrown if A and B have differented sizes or heads
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions called 
					from inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray object C = A + iB or C = A * exp(iB) from 2 real XArray objects.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter A
		\par		Example:
\verbatim
XArray1D<float> A(20, 1.0f), B(20, 2.0f);
XArray1D<fcomplex> C; 
MakeComplex(A, B, C, false);
\endverbatim
	*/
	template <class T> void MakeComplex(const XArray<T>& A, const XArray<T>& B, XArray< std::complex<T> >& C, bool bMakePolar)
	{
		if (A.size() != B.XArrayBase<T>::size()) throw std::invalid_argument("invalid_argument 'A and B' in MakeComplex(A, B, C) (different sizes)"); 	
		if (A.GetHeadPtr()) A.GetHeadPtr()->Validate();
		if (!A.SameHead(B))
			throw std::invalid_argument("invalid_argument 'A and B' in MakeComplex(A, B, C) (different heads)"); 	

		C.resize(A.XArrayBase<T>::size());
		if (bMakePolar) for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) C[i] = std::polar<T>(A[i], B[i]);
		else for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) C[i] = std::complex<T>(A[i], B[i]);
		C.SetHeadPtr(A.GetHeadPtr() ? A.GetHeadPtr()->Clone() : 0);
	}
	

	//---------------------------------------------------------------------------
	//Function XArray<T>::Re
	//
	//	Calculates the real part of a complex XArray object
	//
	/*!
		\brief		Calculates the real part of a complex XArray object
		\param		C	Input complex XArray object
		\param		A	Output real XArray object to be made equal to the real part of C
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the real part of every element of a complex XArray object and 
			copies the result into a real XArray object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0));
XArray1D<double> A; 
Re(C, A);
\endverbatim
	*/
	template <class T> void Re(const XArray< std::complex<T> >& C, XArray<T>& A)
	{ 
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();
		A.resize(C.XArrayBase<T>::size()); 
		for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) A[i] = std::real(C[i]);
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);
	}
	
	
	//---------------------------------------------------------------------------
	//Function XArray<T>::Im
	//
	//	Calculates the imaginary part of a complex XArray object
	//
	/*!
		\brief		Calculates the imaginary part of a complex XArray object
		\param		C	Input complex XArray object
		\param		A	Output real XArray object to be made equal to the imaginary part of C
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the imaginary part of every element of a complex XArray object
			and copies the result into a real XArray object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0));
XArray1D<double> A; 
Im(C, A);
\endverbatim
	*/
	template <class T> void Im(const XArray< std::complex<T> >& C, XArray<T>& A)
	{ 
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();
		A.resize(C.XArrayBase<T>::size()); 
		for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) A[i] = std::imag(C[i]);
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);
	}
	
	
	//---------------------------------------------------------------------------
	//Function XArray<T>::Abs
	//
	//	Calculates the modulus of a complex XArray object
	//
	/*!
		\brief		Calculates the modulus of a complex XArray object
		\param		C	Input complex XArray object
		\param		A	Output real XArray object to be made equal to the modulus of C
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the modulus of every element of a complex XArray object
			and copies the result into a real XArray object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0));
XArray1D<double> A;
Abs(C, A);
\endverbatim
	*/	
	template <class T> void Abs(const XArray< std::complex<T> >& C, XArray<T>& A)
	{ 
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();
		A.resize(C.XArrayBase<T>::size());
		for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) A[i] = T(::sqrt(C[i].real() * C[i].real() + C[i].imag() * C[i].imag()));
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);	
	}
	
	
	//---------------------------------------------------------------------------
	//Function XArray<T>::Arg
	//
	//	Calculates the argument of a complex XArray object
	//
	/*!
		\brief		Calculates the argument of a complex XArray object
		\param		C	Input complex XArray object
		\param		A	Output real XArray object to be made equal to the argument of C
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the argument of every element of a complex XArray object and copies
			the result into a real XArray object. The head for the output object A is copied from the input object C
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0));
XArray1D<double> A;
Arg(C, A);
\endverbatim
	*/	
	template <class T> void Arg(const XArray< std::complex<T> >& C, XArray<T>& A)
	{ 
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();
		A.resize(C.XArrayBase<T>::size()); 
		for (index_t i = 0; i < C.XArrayBase<T>::size(); i++) A[i] = T(::atan2(C[i].imag(), C[i].real()));
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);
	}
	
	
	//---------------------------------------------------------------------------
	//Function XArray<T>::Abs2
	//
	//	Calculates the square modulus of a complex XArray object
	//
	/*!
		\brief		Calculates the square modulus of a complex XArray object
		\param		C	Input complex XArray object
		\param		A	Output real XArray object to be made equal to the square modulus of C
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the square modulus of every element of a complex XArray object
			and copies the result into a real XArray object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0));
XArray1D<double> A;
Abs2(C, A);
\endverbatim
	*/	
	template <class T> void Abs2(const XArray< std::complex<T> >& C, XArray<T>& A)
	{ 
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();
		A.resize(C.XArrayBase< std::complex<T> >::size());
		for (index_t i = 0; i < C.XArrayBase< std::complex<T> >::size(); i++) A[i] = C[i].real() * C[i].real() + C[i].imag() * C[i].imag();
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);
	}
	
	
	//---------------------------------------------------------------------------
	//Function XArray<T>::Conjg
	//
	//	Conjugates a complex XArray object
	//
	/*!
		\brief		Conjugates a complex XArray object
		\param		C	Input and output complex XArray object
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function conjugates every element of a complex XArray object. The head does not change.
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0));
Conjg(C);
\endverbatim
	*/	
	template <class T> void Conjg(XArray< std::complex<T> >& C) 
	{ 
		T* pFront = reinterpret_cast<T*>(&C[0]);
		pFront--;
		for(index_t i = 0; i < C.XArrayBase<T>::size(); i++) 
		{
			pFront += 2;
			*pFront = -(*pFront);
			//C[i] = std::complex<T>(C[i].real(), -C[i].imag());
		}
	}
	
	
	//---------------------------------------------------------------------------
	//Function XArray<T>::CArg
	//
	//	Calculates the 1D-continuous phase of a complex XArray object
	//
	/*!
		\brief		Calculates the 1D-continuous phase of a complex XArray object
		\param		C	Input complex XArray object
		\param		A	Output real XArray object to be made equal to the phase of C
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the argument of every element of a complex XArray object using the 
			value of the argument of the preceding element to calculate the appropriate 2pi-multiple shift,
			and copies the result into a real XArray object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray1D<dcomplex> C(20, dcomplex(1.0, 2.0)); 
XArray1D<double> A; 
CArg(C, A);
\endverbatim
	*/	
	template <class T> void CArg(const XArray< std::complex<T> >& C, XArray<T>& A) 
	{
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();	
		A[0] = carg(C[0], T(0));
		for (index_t i = 1; i < C.XArrayBase< std::complex<T> >::size(); i++) A[i] = carg(C[i], A[i-1]);
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);
	}


} //namespace xar closed

//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings are generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray<char>;
	template class xar::XArray<short>;
	template class xar::XArray<long>;
	template class xar::XArray<float>;
	template class xar::XArray<double>;
	template class xar::XArray<xar::fcomplex>;
	template class xar::XArray<xar::dcomplex>;
#endif

//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XARRAY_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
