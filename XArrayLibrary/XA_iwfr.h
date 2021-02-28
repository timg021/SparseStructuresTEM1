//Header XA_iwfr.h
//
//
//	HEADER FILE TITLE:
//
//		Phase retrieval algorithms based on the iwfr (Iterative Wave Function Reconstruction) method
//
/*!
	\file		XA_iwfr.h
	\brief		Phase retrieval algorithms based on the iwfr (Iterative Wave Function Reconstruction) method
	\par		Description:
		This header contains a class that provides phase retrieval services 
		based on the iwfr (Iterative Wave Function Reconstruction) method for XArray2D<T> objects
*/

#pragma once

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_spln2.h"
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
//Class XA_IWFR<T>
//
//	Phase retrieval algorithms based on the 1st iwfr and Rytov approximations
//
/*!
	\brief		Phase retrieval algorithms based on the iwfr (Iterative Wave Function Reconstruction) method
	\par		Description:
				This class template provides iwfr (Iterative Wave Function Reconstruction) 
				services for XArray2D<T> objects
	\remarks	An object of this class represents an interface exposing several functions
				that provide iwfr (Iterative Wave Function Reconstruction) phase retrieval services				
	\remarks	This class can only be created for T=float or T=double
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XA_IWFR
	{
	// Typedefs
	public:
		typedef T type;
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XA_IWFR() {}
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XA_IWFR(const XA_IWFR<T>& rCopy) { GetValuetype(); }
	public:
		//! Destructor
		~XA_IWFR() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XA_IWFR<T>& rCopy) { if (this == &rCopy) return; }

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XA_IWFR<T> objects for unsupported types T
		//! Returns the xar::_eValueType corresponding to T
		static xar::_eValueType GetValuetype(void);

	// Operations
	public:
		//! Retrieves complex wave function from a vector of defocused images
		void Iwfr(vector<XArray2D<T> >& vint0, XArray2D<std::complex<T> >& campOut, vector<double> vdefocusdist, double zmiddle, double k2maxo, double Cs3, double Cs5, int kmax, double epsilon, bool bVerboseOutput);
		void CTFL2(vector<XArray2D<T> >& vint0, XArray2D<T>& fiOut, vector<double> vdefocusdist, double q2max, double Cs3, double Cs5, double& alpha);
		void PhaseB7(const XArray2D<T>& intIn, XArray2D<T>& fiOut);
		void MinLogAmp(const XArray2D<T>& intIn, XArray2D<T>& fiOut);

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
	// Member functions

	};

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//! Returns the xar::_eValueType corresponding to T=float
template<> inline xar::_eValueType XA_IWFR<float>::GetValuetype() { return xar::eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> inline xar::_eValueType XA_IWFR<double>::GetValuetype() { return xar::eXADouble; }

//! Retrieves phase from several defocused images using Iterative Wave Function Reconstruction method
// vint0 - input vector of defocused images
// campOut - output complex amplitude in the plane zmiddle
// vdefocusdist - vector of defocus distances
// zmiddle - z-coordinate of the defocus plane where the output is produced
// k2maxo - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// kmax - maximum number of iterations
// epsilon - minimum L2 error (interrupts the iterations)
// bVerboseOutput - if TRUE, then verbose diagnostic output is printed
template <class T> void XA_IWFR<T>::Iwfr(vector< XArray2D<T> >& vint0, XArray2D<std::complex<T> >& campOut, vector<double> vdefocusdist, double zmiddle, double k2maxo, double Cs3, double Cs5, int kmax, double epsilon, bool bVerboseOutput)
{
	int ndefocus = (int)vint0.size();
	if (ndefocus < 1) throw std::exception("Error: input image vector empty in IWFR()");

	bool bAbort(false);
	double dtemp; // auxilliary variable
	double ssej(0.0), ssejm1(0.0); // current and previous average reconstruction errors
	vector<double> verr(ndefocus); // reconstruction errors in individual defocus planes
	vector<double> vint0_L1(ndefocus); // L1 norms of the initial defocused intensities
	vector<XArray2D<T> > vint(ndefocus); // iterated defocused intensities
	vector<XArray2D<std::complex<T> > > vcamp(ndefocus); // defocused complex amplitudes

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(vint0[0].GetHeadPtr()->Clone());
	IXAHWave2D* ph2 = GetIXAHWave2D(vint0[0]);
	ph2->Validate();
	index_t ny = vint0[0].GetDim1();
	index_t nx = vint0[0].GetDim2();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);

	for (int k = 0; k < kmax; k++)
	{
		// apply initial or newly reconstructed phases (the second case is equal to restoring the original moduli)
		// and propagate each defocused amplitude to the "middle" plane z = zmiddle
		#pragma omp parallel for
		for (int n = 0; n < ndefocus; n++)
		{
			try
			{
				if (k == 0) // create initial defocused complex amplitude
				{
					vint0_L1[n] = vint0[n].Norm(eNormL1);
					if (bVerboseOutput) printf("\nL1 norm of input defocused intensity no. %d = %g", n, vint0_L1[n]);
					if (vint0_L1[n] == 0) throw std::exception("Error: input intensity file is empty in IWFR()");
					vint0[n] ^= 0.5; // intensity --> real amplitude
					MakeComplex(vint0[n], 0.0, vcamp[n], true); // apply initial zero phases
				}
				else // apply phases obtained on the previous iteration
					for (index_t j = 0; j < vcamp[n].GetDim1(); j++)
						for (index_t i = 0; i < vcamp[n].GetDim2(); i++)
						{
							dtemp = abs(vcamp[n][j][i]);
							if (dtemp) vcamp[n][j][i] *= vint0[n][j][i] / dtemp;
							else vcamp[n][j][i] = std::polar(vint0[n][j][i], 0.0);
						}
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(zmiddle - vdefocusdist[n], false, k2maxo, Cs3, Cs5); // propagate to z = 0
			}
			catch (std::exception& E)
			{
				printf("\n\n!!!Exception: %s\n", E.what());
				bAbort = true;
			}
		}
		if (bAbort) throw std::runtime_error("at least one thread has thrown an exception in IWFR().");

		// average complex amplitudes in the middle plane
		campOut = vcamp[0];
		for (index_t n = 1; n < ndefocus; n++) campOut += vcamp[n];
		campOut /= double(ndefocus);

		// propagate the averaged complex amplitude from the middle plane to the individual defocus planes
		#pragma omp parallel for shared(campOut)
		for (int n = 0; n < ndefocus; n++)
		{
			try
			{
				vcamp[n] = campOut;
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(vdefocusdist[n] - zmiddle, false, k2maxo, Cs3, Cs5); // propagate to z = z[n]
				Abs(vcamp[n], vint[n]);
				vint[n] -= vint0[n];
				verr[n] = pow(vint[n].Norm(eNormL2), 2.0) / vint0_L1[n];
			}
			catch (std::exception& E)
			{
				printf("\n\n!!!Exception in IWFR(): %s\n", E.what());
				bAbort = true;
			}
		}
		if (bAbort) throw std::runtime_error("at least one thread has thrown an exception in IWFR().");

		// calculate the current reconstruction error and
		// if the difference with the previous error is smaller than the defined minimum
		// or, if the error started to increase, interrupt the iterations
		ssej = 0.0;
		for (index_t n = 0; n < ndefocus; n++) ssej += verr[n];
		ssej /= double(ndefocus);

		if (bVerboseOutput)
		{
			if (k == 0) printf("\nIteration number %d; SSE_aver error at 0th iteration = %g\n", k, ssej);
			else printf("\nIteration number %d; SSE_aver error difference with previous iteration = %g\n", k, ssejm1 - ssej);

			for (index_t n = 0; n < ndefocus; n++) printf("SSE(%zd) = %g ", n, verr[n]);
		}

		if (k > 0 && (ssejm1 - ssej) < epsilon) break;
		else ssejm1 = ssej;
	}
}

//! Retrieves phase from several defocused images using L2 minimization of the CTF function
// vint0 - input vector of defocused images
// fiOut - output phase in the plane vdefocusdist[0]
// vdefocusdist - vector of defocus distances
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// alpha - Tikhonov regularization parameter, on exit it is replaced by 0.1 times the minimal non-zero CTF^4 (as in the denominator of eq.(A6) of paper D. PAGANIN et al, J.Micros. 214 (2004) 51-61)
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
// NOTE: it could be useful later to implement a modified version of this algorithm, where the phase is not simply retrieved at defocus[0], which is highly dependent on the image[0], 
// but instead it would be retrieved at each defocus[n], and then the complex amplitude would be averaged, e.g. at the middle plane, as in IWFR
template <class T> void XA_IWFR<T>::CTFL2(vector<XArray2D<T> >& vint0, XArray2D<T>& fiOut, vector<double> vdefocusdist, double q2max, double Cs3, double Cs5, double& alpha)
{
	bool bAper(q2max > 0);
	int ndefocus = (int)vint0.size();
	if (ndefocus < 2) throw std::exception("Error: at least two defocused images are required in CTFL2()");

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(vint0[0].GetHeadPtr()->Clone());
	IXAHWave2D* ph2 = GetIXAHWave2D(vint0[0]);
	ph2->Validate();
	index_t ny = vint0[0].GetDim1();
	index_t nx = vint0[0].GetDim2();
	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	// convert input defocused images to contrast functions and complexify for FFT (should change the code to use real FFT at a later time)
	vector<XArray2D<std::complex<T> > > vcamp(ndefocus);
	for (int n = 0; n < ndefocus; n++)
	{
		vint0[n] -= vint0[n].Norm(eNormAver); // image --> contrast function * (-1)
		MakeComplex(vint0[n], T(0), vcamp[n], false);
	}

	index_t i = 2;
	while (i < ny) i *= 2;
	if (i != ny) throw std::invalid_argument("invalid_argument in CTFL2() (m_dim1 is not a power of 2)");
	index_t j = 2;
	while (j < nx) j *= 2;
	if (j != nx) throw std::invalid_argument("invalid_argument in CTFL2() (m_dim2 is not a power of 2)");

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::CTFL2 (evanescent waves present)");

	vector<double> dblDistance(ndefocus);
	dblDistance[0] = 0.0;
	for (int n = 1; n < ndefocus; n++)
	{
		dblDistance[n] = vdefocusdist[n] - vdefocusdist[0];
		if (dblDistance[n] == 0)
			throw std::runtime_error("runtime_error in XA_IWFR<T>::CTFL2 (duplicate propagation distance)");
	}

	//********* Fourier transforming initial amplitude u_n(i,j)
	vector<T*> u(ndefocus);
	OouraFft<T> fft;
	for (int n = 0; n < ndefocus; n++)
	{
		u[n] = reinterpret_cast<T*>(&(vcamp[n].front()));
		fft.Complex2D((std::complex<T> *) u[n], ny, nx, OouraFft<T>::eDirFwd);
	}

	// prepare output phase array
	fiOut.Resize(ny, nx, T(0));
	fiOut.SetHeadPtr(vint0[0].GetHeadPtr()->Clone());
	XArray2D<std::complex<T> > vcampOut;
	MakeComplex(fiOut, T(0), vcampOut, false);
	T* fiC;
	fiC = reinterpret_cast<T*>(&(vcampOut.front()));
	fft.Complex2D((std::complex<T> *) fiC, ny, nx, OouraFft<T>::eDirFwd);

	//********* Dividing F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi2 = 1.0 / xap2;
	double deta2 = 1.0 / yap2;

	vector<double> fac2(ndefocus);
	double fac2min = PI * wl * abs(dblDistance[1]);
	for (int n = 1; n < ndefocus; n++)
	{
		fac2[n] = PI * wl * dblDistance[n];
		if (abs(fac2[n]) < fac2min) fac2min = abs(fac2[n]);
	}
	double alphaNew = 0.1 * pow(fac2min * std::min(dcsi2, deta2), 4);
	if (alpha <= 0) alpha = alphaNew;

	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	index_t k, kj;
	double eta2, q2;
	double dtemp, sintemp, sumsin2, costemp, sumuk, sumuk1, Cstemp(0);

	for (long i = -long(nyd2); i < 0; i++)
	{
		kj = nxy2 + nx2 * i + nx2;
		eta2 = deta2 * i * i;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3); 
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
		kj = nxy2 + nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
	}
	for (long i = 0; i < long(nyd2); i++)
	{
		kj = nx2 * i + nx2;
		eta2 = deta2 * i * i;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
		kj = nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
	}
	alpha = alphaNew; // pass the information about the smallest possible non-zero squared squared CTF value back to the calling program in order to optionally adjust alpha next time

	//********* inverse Fourier transforming
	fft.Complex2D((std::complex<T> *) fiC, ny, nx, OouraFft<T>::eDirInv);
	T fact = T(1.0) / nxy;
	for (k = 0; k < nxy2; k++)	fiC[k] *= fact;
	Re(vcampOut, fiOut);
}


//! Retrieves phase from a single defocused image as phase = -0.5 * log(Intensity)
// vint0 - input defocused image
// fiOut - output phase in the same plane
template <class T> void XA_IWFR<T>::MinLogAmp(const XArray2D<T>& intIn, XArray2D<T>& fiOut)
{
	const double* pInt = &(intIn[0][0]);
	fiOut = intIn;
	double* pPha = &(fiOut[0][0]);
	for (index_t i = 0; i < fiOut.size(); i++)
		//pPha[i] = -0.5 * log(pInt[i]);
		pPha[i] = 0.5 * (1.0 - pInt[i]);
}


//! Retrieves phase from a single defocused image according to eq. (B7) in the our 2nd Ultramicroscopy paper
// vint0 - input defocused image
// fiOut - output phase in the same plane
template <class T> void XA_IWFR<T>::PhaseB7(const XArray2D<T>& intIn, XArray2D<T>& fiOut)
{
		const double* pInt = &(intIn[0][0]);
		fiOut = intIn;
		double* pPha = &(fiOut[0][0]);
		double Kmax2(0.0), dtemp;
		for (index_t i = 0; i < intIn.size(); i++)
		{
			dtemp = 1.0 - pInt[i];
			dtemp *= dtemp;
			pPha[i] = dtemp;
			if (dtemp > Kmax2) Kmax2 = dtemp;
		}
		for (index_t i = 0; i < fiOut.size(); i++)
			pPha[i] = -0.5 * sqrt(Kmax2 - pPha[i]);
}


} // namespace xar closed

//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//

// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class XA::XA_IWFR<float>;
	template class XA::XA_IWFR<double>;
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
