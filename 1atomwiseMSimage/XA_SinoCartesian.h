//Header XA_SinoCartesian.h
//
//
//	HEADER FILE TITLE:
//
//		Converts sinograms into Cartesian grid and back
//
/*!
	\file		XA_SinoCartesian.h
	\brief		Converts sinograms into Cartesian grid and back
	\par		Description:
		Converts sinograms into Cartesian grid and back
*/
#pragma once
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_head2.h"

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
//Class XArraySinoCart
//
//	Two-dimensional geometrical transformations
//
/*!
	\brief		Converts sinograms into Cartesian grid and back
	\par		Description:
				This class template defines a 'wrapper' around the XArray2D<T> object
				on which it operates; it contains functions implementing convers sinograms into Cartesian grid and back
	\remarks    If a Wavehead2D is attached to the XArray2D<T> object, it is transformed too
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
*/
	template <class T> class XArraySinoCart
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArraySinoCart(XArray2D<T>& rXAr2D) : m_rXArray2D(rXAr2D) {}
	protected:
		//! Copy constructor (declared protected to prohibit explicit copying)
		XArraySinoCart(const XArraySinoCart<T>& rCopy) : m_rXArray2D(rCopy.m_rXArray2D) {}
	public:
		//! Destructor 
		~XArraySinoCart() {}
	
	// Operators
	protected:
		//! Assignment (declared protected to prohibit explicit copying)
		void operator=(const XArraySinoCart<T>& rCopy);

	// Attributes
	public:
		//! Returns a reference to the non-modifiable 'wrapped' XArray2D<T> object
		const XArray2D<T>& GetBaseObject() const { return m_rXArray2D; }
		//! Returns a reference to the 'wrapped' XArray2D<T> object
		XArray2D<T>& GetBaseObject() { return m_rXArray2D; }

	// Operations
	public:
		//!	Converts the array from Cartesian to polar coordinates
		void Cart2Polar(bool& bOdd);
		//!	Converts the array from polar to Cartesian coordinates
		void Polar2Cart(bool bOdd);

	// Overridables

	// Member variables	
		//! Reference to the XArray2D<T> object that is being operated upon
		XArray2D<T>& m_rXArray2D;
	};


	//---------------------------------------------------------------------------
	//	TEMPLATE MEMBER DEFINITIONS
	//

	//! Assignment  (declared protected to prohibit explicit copying)
	template <class T> void XArraySinoCart<T>::operator=(const XArraySinoCart<T>& rCopy)
	{ 
		if (this == &rCopy) 
			return; 
		else 
			m_rXArray2D = rCopy.m_rXArray2D;
	}

	
	//---------------------------------------------------------------------------
	//Function XArraySinoCart<T>::Cart2Polar
	//
	//	In-place conversion of the array from Cartesian to polar coordinates
	//
	/*!
		\brief		Performs in-place conversion of the array from Cartesian to polar coordinates
		\param		bOdd (out) is an auxiliary parameter that specifies whether the size of the original array if odd (bOdd = true) or even (bOdd = false)
		\par		Description:
			This function performes in-place conversion of the array from the Cartesian to polar coordinates			
	*/	
	template<class T> void XArraySinoCart<T>::Cart2Polar( bool& bOdd )
	{
		index_t nx( m_rXArray2D.GetDim2() );

		if( m_rXArray2D.GetDim1() != nx )
			throw std::invalid_argument( "Invalid argument 'm_rXArray2D' in XArraySinoCart<T>::Cart2Polar() (MUST be square)" );

		IXAHWave2D* ph2( GetIXAHWave2D( m_rXArray2D ) );

		if( ph2 != 0 && ph2->GetXStep( nx ) != ph2->GetYStep( nx ) )
			throw std::invalid_argument( "Invalid argument 'm_rXArray2D' in XArraySinoCart<T>::Cart2Polar() (MUST be square)" );

		bOdd = !(( nx % 2 ) == 0);

		index_t nr( bOdd ? ( nx - 1 ) / 2 + 1 : nx / 2 );
		index_t nphi( static_cast<index_t>( xar::tPI * nr + 0.5 ) );
		double dphi( xar::tPI / nphi );

		xar::XArray2D<T> XArPolar( nphi, nr );
		double num2( 0.5 * ( nx - 1 ) );

		for( index_t i = 0; i < nphi; i++ )
		{
			double phi( i * dphi - xar::PI2 ); //@@@

			for( index_t j = 0; j < nr; j++ )
			{
				double r( static_cast<double>( j ) );
				double x( r * cos( phi ) + num2 );
				double y( r * sin( phi ) + num2 );

				index_t ix( static_cast<index_t>( x ) );
				index_t iy( static_cast<index_t>( y ) );
				index_t ix1( ix + 1 );
				index_t iy1( iy + 1 );

				if( ix1 < nx && iy1 < nx )
				{
					double dx( x - static_cast<double>( ix ) );
					double dy( y - static_cast<double>( iy ) );
					double a00( ( 1. - dy ) * ( 1. - dx ) );
					double a01( ( 1. - dy ) * dx );
					double a10( dy * ( 1. - dx ) );
					double a11( dy * dx );
					XArPolar[i][j] = static_cast<T>( m_rXArray2D[iy][ix] * a00 + m_rXArray2D[iy][ix1] * a01 + m_rXArray2D[iy1][ix] * a10 + m_rXArray2D[iy1][ix1] * a11);
				}
			}
		}

		// if m_rXArray2D has a header, create a header
		if( ph2 != 0 )
		{
			IXAHWave2D* phNew( CreateWavehead2D() );
			phNew->SetData( ph2->GetWl(), 0, dphi * ( nphi - 1 ), 0, ph2->GetXStep( nx ) * ( nr - 1 ) );
			XArPolar.SetHeadPtr( phNew );
		}
		
		m_rXArray2D.Swap( XArPolar );
	}


	//---------------------------------------------------------------------------
	//Function XArraySinoCart<T>::Polar2Cart
	//
	//	In-place conversion of the array from polar to Cartesian coordinates
	//
	/*!
		\brief		Performs in-place conversion of the array from polar to Cartesian coordinates
		\param		bOdd (in) is an auxiliary parameter that specifies whether the size of the array (in the Cartesian coordinates) if odd (bOdd = true) or even (bOdd = false)
		\par		Description:
			This function performes in-place conversion of the array from the polar to the Cartesian coordinates			
	*/	
	template<class T> void XArraySinoCart<T>::Polar2Cart( bool bOdd )
	{
		index_t nphi( m_rXArray2D.GetDim1() );
		index_t nr( m_rXArray2D.GetDim2() );

		index_t nx( bOdd ? 2 * nr - 1 : 2 * nr );
		double invdphi( nphi / xar::tPI );

		xar::XArray2D<T> XArCart( nx, nx ); // a square array
		double num2( 0.5 * ( nx - 1 ) );

		for( index_t i = 0; i < nx; i++ )
		{
			double y( static_cast<double>( i ) - num2 );

			for( index_t j = 0; j < nx; j++ )
			{
				double x( static_cast<double>( j ) - num2 );
				double ro( sqrt( x * x + y * y ) );
				index_t iro( static_cast<index_t>( ro ) );//index_t(bOdd ? ro : ro - 0.5);
				index_t iro1( iro + 1 );

				if( iro1 < nr )
				{
					double phi;

					if( ro == 0 )
						phi = 0;
					else
					{
						phi = acos( x / ro );

						if( y < 0 )
							phi = xar::tPI - phi;

						phi += xar::PI2; //@@@

						if( phi >= xar::tPI )
							phi -= xar::tPI;

						//phi -= int(phi / xar::tPI) * xar::tPI;
						//if (phi < 0) phi += xar::tPI;
						phi *= invdphi;
					}

					double d_ro( ro - static_cast<double>( iro ) );
					index_t iphi( static_cast<index_t>( phi ) );
					index_t iphi1( iphi + 1 );

					if( iphi1 >= nphi )
						iphi1 -= nphi;

					double d_phi( phi - static_cast<double>( iphi ) );
					double a00( ( 1. - d_phi ) * ( 1. - d_ro ) );
					double a01( ( 1. - d_phi ) * d_ro );
					double a10( d_phi * ( 1. - d_ro ) );
					double a11( d_phi * d_ro );

					XArCart[i][j] = static_cast<T>( m_rXArray2D[iphi][iro] * a00 + m_rXArray2D[iphi][iro1] * a01 + m_rXArray2D[iphi1][iro] * a10 + m_rXArray2D[iphi1][iro1] * a11 );
				}
			}
		}

		IXAHWave2D* ph2( GetIXAHWave2D( m_rXArray2D ) );

		// if m_rXArray2D has a header, create a header
		if( ph2 != 0 )
		{
			IXAHWave2D* phNew( CreateWavehead2D() );
			phNew->SetData( ph2->GetWl(), 0, ph2->GetXStep( nr ) * ( nx - 1 ), 0, ph2->GetXStep( nr ) * (nx - 1) );
			XArCart.SetHeadPtr( phNew );
		}
		
		m_rXArray2D.Swap( XArCart );
	}

	}  //namespace xar closed
//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings are generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.


// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArraySinoCart<char>;
	template class xar::XArraySinoCart<short>;
	template class xar::XArraySinoCart<long>;
	template class xar::XArraySinoCart<float>;
	template class xar::XArraySinoCart<double>;
	template class xar::XArraySinoCart<xar::fcomplex>;
	template class xar::XArraySinoCart<xar::dcomplex>;
#endif

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
