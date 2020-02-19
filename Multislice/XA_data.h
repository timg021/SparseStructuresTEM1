//Header XA_data.h
//
//
//	HEADER FILE TITLE:
//
//		Data file I/O facilities
//
/*!
	\file		XA_data.h
	\brief		Data file I/O facilities
	\par		Description:
		This header contains functions that provide file I/O for XArray objects and
		DAT, GRD and GRC data files
*/
#if !defined XA_DATA_H
#define XA_DATA_H

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include <limits>
#include "XA_head1.h"
#include "XA_head2.h"
#include "XA_file.h"


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
//	ENUMERATED DataA TYPES
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
//Class XArData
//
//	Data file I/O facilities for XArrays	
//
/*!
	\brief		Data file I/O facilities for XArrays
	\par		Description:
				This class provides data file I/O facilities for XArrayND<T> objects
				and DAT, GRD, and GRC files
	\remarks	An object of this class represents a simple wrapper around several template 
				functions that provide data file I/O for XArrays				
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	class XArData
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArData() {}
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XArData(const XArData& rCopy) {}
	public:
		//! Destructor
		~XArData() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XArData& rCopy) {}

	// Attributes
	public:
		
	// Operations
	public:

	//Function WriteFileDAT
	//
	//	Writes real XArray1D objects into ASC Data files
	//
	/*!
		\brief		Writes real XArray1D objects into an ASC Data files
		\param		rXAr1D	Reference to an XArray1D<T> object to be saved
		\param		pchFilename	Full name of the file to be written into
		\exception  std::invalid_argument is thrown if XArray1D object contains complex values
		\exception  std::runtime_error is thrown if any of the file write operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function writes a real 1D XArray object into an ASC Data file.
			The head of the object, except wavelength, is transformed into the argument (first) column.
			The array dataa is always stored as "float" values.
			The Data file consists of 2 columns of "float" values:
			(1) the column of arguments, x
			(2) the column of corresponding array values
	*/	
	template <class T> static void WriteFileDAT(const XArray1D<T>& rXAr1D, const char* pchFilename)
	{
		double dblXlo, dblStep;

		const IXAHWave1D* pHead = GetIXAHWave1D(rXAr1D);
		if (pHead)
		{
			dblXlo = pHead->GetXlo();
			dblStep = GetStep(rXAr1D);
		}
		else
		{
			dblXlo = 1;
			dblStep = 1;
		}

		// open file
		FilePtr fp(pchFilename, "wt");

		// write array
		for (index_t i = 0; i < rXAr1D.size(); i++)
		{
			if (fprintf(fp, "%14.6e\t%14.6e\n", dblXlo + dblStep * i, rXAr1D[i]) <= 0) 
				throw std::runtime_error("runtime_error in XArData::WriteFileDAT (error writing array)");
		}
	}

	//! Specialization for T = fcomplex (only throws an exception)
	static void WriteFileDAT(const XArray1D<fcomplex>& rXAr1D, const char* pchFilename) 
	{
			throw std::invalid_argument("invalid_argument 'rXAr1D' in XArData::WriteFileDAT (fcomplex dataa)");
	}

	//! Specialization for T = dcomplex (only throws an exception)
	static void WriteFileDAT(const XArray1D<dcomplex>& rXAr1D, const char* pchFilename)
	{
			throw std::invalid_argument("invalid_argument 'rXAr1D' in XArData::WriteFileDAT (dcomplex dataa)");
	}


	//---------------------------------------------------------------------------
	//Function ReadFileDAT
	//
	//	Reads real XArray1D objects from an ASC Data file
	//
	/*!
		\brief		Reads real XArray1D objects from an ASC Data file
		\param		rXAr1D	Reference to an XArray1D<T> object to be read in
		\param		pchFilename	Full name of the file to be read from
		\param		dblWavelength	Wavelength to be used in the Wavehead1D
		\param		bCheck	Determines if to perform the check for the equidistance of x-dataa
		\exception  std::invalid_argument is thrown if XArray1D object has a complex value type
		\exception  std::runtime_error is thrown if any of the file read operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function reads a 1D XArray object from an ASC Data file.
			The head of the object, except wavelength, is read from the file.
			The Data file consists of 2 columns of "float" values:
			(1) the column of arguments, x
			(2) the column of corresponding array values
	*/	
	template <class T> static void ReadFileDAT(XArray1D<T>& rXAr1D, const char* pchFilename, double dblWavelength, bool bCheck = true)
	{
		char buf[129], buf1[129], buf2[129];

		// open file for reading as ASCII text
		FilePtr fp(pchFilename, "rt");

		// first, count the number of rows in the file
		index_t i, iSize = 0;
		for (;;) 
		{ 
			if (fgets(buf, 128, fp) <= 0)
				break;
			else
			{
				iSize++;
			}
		}
		// now, read the file dataa into rXAr1D
		rewind (fp);
		rXAr1D.Resize(iSize);
		XArray1D<double> temp(iSize); // temp buffer for x-dataa
		// read array
		for (i = 0; i < iSize; i++)
		{
			fgets(buf, 128, fp);
			sscanf(buf, "%s %s", buf1, buf2);
			temp[i] = strtod(buf1, 0);
			rXAr1D[i] = T(strtod(buf2, 0));
		}
		if (bCheck)
		{
			// Check whether the x-dataa set is equidistant
			double dblStep = (temp[iSize - 1] - temp[0]) / (iSize - 1);
			for (i = 1; i < iSize; i++)
			{
				if (fabs(temp[i] - temp[i - 1] - dblStep) > dblStep * 0.001)
					throw std::invalid_argument("invalid_argument 'rXAr1D' in XArData::ReadFileDAT (not equidistant)");
			}
		}
		IXAHWave1D* ph1 = CreateWavehead1D();
		ph1->SetData(dblWavelength, temp[0], temp[iSize - 1]);
		rXAr1D.SetHeadPtr(ph1);
	}

	//! Specialization for T = fcomplex (only throws an exception)
	static void ReadFileDAT(XArray1D<fcomplex>& rXAr1D, const char* pchFilename, double dblWavelength, bool bCheck)
	{
			throw std::invalid_argument("invalid_argument 'rXAr1D' in XArData::ReadFileDAT (fcomplex object)");
	}

	//! Specialization for T = dcomplex (only throws an exception)
	static void ReadFileDAT(XArray1D<dcomplex>& rXAr1D, const char* pchFilename, double dblWavelength, bool bCheck)
	{
			throw std::invalid_argument("invalid_argument 'rXAr1D' in XArData::ReadFileDAT (dcomplex object)");
	}


	//---------------------------------------------------------------------------
	//Function WriteFileGRD
	//
	//	Writes real XArray2D object into an ASC or BIN GRD file
	//
	/*!
		\brief		Writes real XArray2D object into an ASC or BIN GRD file
		\param		rXAr2D	Reference to an XArray2D<T> object to be saved
		\param		pchFilename	Full name of the file to be written into
		\param		filetype eGRDBIN or eGRDASC for BIN and ASC files respectively
		\exception  std::invalid_argument is thrown if XArray2D object contains complex values
		\exception  std::invalid_argument is thrown if XArray2D object does not have a Wavehead2D
		\exception  std::invalid_argument is thrown if 'filetype' parameter is not eGRDBIN or eGRDASC
		\exception  std::runtime_error is thrown if any of the file write operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function writes a real 2D XArray object into a BIN or ASC GRD file (GRD file format
			is used	e.g. in SURFER program by Golden Software). The head of the object, except wavelength, 
			is saved in the file. The array data is always stored as "float" values. The GRD file
			contains:
			(1) string "DSAA" or "DSBB" (in ASC and BIN GRD files, respectively)
			(2) nx ny  - array dimensions as "short"s separated by a white space
			(3) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a white space)
			(4) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a white space)
			(5) zlo zhi - mnimum and maximum array value as "double"s (in ASC format, separated by a white space)
			(6) u[i][j] - array values as "float"s (in ASC format, separated by a white space), 
							j index changes most rapidly
	*/	
	template <class T> static void WriteFileGRD(const XArray2D<T>& rXAr2D, const char* pchFilename, _eFileType filetype)
	{
		if (!rXAr2D.GetHeadPtr())
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRD (no head)");
		const IXAHWave2D* ph2 = dynamic_cast<const IXAHWave2D*>(rXAr2D.GetHeadPtr());
		if (!ph2)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRD (the head is not a Wavehead2D)");
		ph2->Validate();

		if (rXAr2D.GetValuetype() == eXAFComplex || rXAr2D.GetValuetype() == eXADComplex)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRD (complex data)");

		// the following statements will throw exception if there is bad data in the array
		double dblDataMax = rXAr2D.Norm(eNormMax);
		double dblDataMin = rXAr2D.Norm(eNormMin);

		index_t ny = rXAr2D.GetDim1();
		index_t nx = rXAr2D.GetDim2();
		double wl = ph2->GetWl(); 
		double ylo = ph2->GetYlo();
		double yhi = ph2->GetYhi();
		double xlo = ph2->GetXlo();
		double xhi = ph2->GetXhi();

		index_t i, j;
		short stemp;
		//float ftemp;

		// open file
		FilePtr fp;
		int fd;
		switch (filetype)
		{
		case eGRDBIN: 
			//fp.Open(pchFilename, "wb"); 
			fd = _open(pchFilename, _O_WRONLY | _O_BINARY | _O_CREAT | _O_TRUNC | _O_SEQUENTIAL, _S_IWRITE); 
			if (fd == -1) 
			{
				throw std::runtime_error("cannot open file for writing in XArData::WriteFileGRD");
			}
			break;
		case eGRDASC: 
			fp.Open(pchFilename, "wt"); 
			break;
		default: throw std::invalid_argument("invalid_argument 'filetype' in XArData::WriteFileGRD");
		}

		// write GRD head
		switch (filetype)
		{
		case eGRDBIN:
			if (_write(fd, "DSBB", 4) != 4)
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing DSBB)");
			}
			stemp = (short)(nx);
			if (_write(fd, &stemp, unsigned(sizeof(short))) != int(sizeof(short)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing nx)");
			}
			stemp = (short)(ny);
			if (_write(fd, &stemp, unsigned(sizeof(short))) != int(sizeof(short)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing ny)");
			}
			if (_write(fd, &xlo, unsigned(sizeof(double))) != int(sizeof(double)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing xlo)");
			}
			if (_write(fd, &xhi, unsigned(sizeof(double))) != int(sizeof(double)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing xhi)");
			}
			if (_write(fd, &ylo, unsigned(sizeof(double))) != int(sizeof(double)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing ylo)");
			}
			if (_write(fd, &yhi, unsigned(sizeof(double))) != int(sizeof(double)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing yhi)");
			}
			if (_write(fd, &dblDataMin, unsigned(sizeof(double))) != int(sizeof(double)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing zlo)");
			}
			if (_write(fd, &dblDataMax, unsigned(sizeof(double))) != int(sizeof(double)))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing zhi)");
			}
			break;
		case eGRDASC:
			if (fprintf(fp, "%s\n", "DSAA") <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing DSAA)");
			if (fprintf(fp, "%zi   %zi\n", nx, ny) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing nx or ny)");
			if (fprintf(fp, "%.15g   %.15g\n", xlo, xhi) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing xlo or xhi)");
			if (fprintf(fp, "%.15g   %.15g\n", ylo, yhi) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing ylo or yhi)");
			if (fprintf(fp, "%.15g   %.15g\n", dblDataMin, dblDataMax) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing ylo or yhi)");
			break;
		}

		// write array
		switch (filetype)
		{
		case eGRDBIN:
			{
				float* arr = 0;
				if (rXAr2D.GetValuetype() != eXAFloat) // a buffer for type conversion is required
				{
					arr = new float[nx]; // the size of this buffer has to be small enough for network transmission
					if (arr == 0) throw std::runtime_error("runtime_error in XArData::WriteFileGRD (failed to allocate memory buffer)");
				}
				for(index_t i = 0; i < ny; i++)
				{
					if (rXAr2D.GetValuetype() != eXAFloat)
						for(index_t j = 0; j < nx; j++)
							//*arr++ = rXAr2D[i][j]; // does not work!!!
							arr[j] = float(rXAr2D[i][j]);
					else
						arr = const_cast<float*>((const float*)(rXAr2D[i]));
					if (_write(fd, arr, unsigned(nx * sizeof(float))) != int(nx * sizeof(float)))
					{
						if (rXAr2D.GetValuetype() != eXAFloat) delete [] arr;
						_close(fd);
						throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing array into file)");
					}
				}
				if (rXAr2D.GetValuetype() != eXAFloat) delete [] arr;
				_close(fd);
			}
			break;
		case eGRDASC:
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
					if (fprintf(fp, "%.6g\n", (float)rXAr2D[i][j]) <= 0) 
						throw std::runtime_error("runtime_error in XArData::WriteFileGRD (error writing array)");
			break;
		}

	}


	//---------------------------------------------------------------------------
	//Function ReadFileGRD
	//
	//	Reads an XArray2D object from an ASC or BIN GRD file
	//
	/*!
		\brief		Reads an XArray2D object from an ASC or BIN GRD file
		\param		rXAr2D	Reference to an XArray2D<T> object to be read in
		\param		pchFilename	Full name of the file to be read from
		\param		dblWavelength	Wavelength to be used in the Wavehead2D
		\param		bForceRead	Determines if the read data may be forced into an XArray2D object with integral value type
		\exception  std::invalid_argument is thrown if bForceRead is false and XArray2D object has integral value type
		\exception  std::runtime_error is thrown if any of the file read operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function reads a 2D XArray object from a BIN or ASC GRD file (GRD file format is used
			e.g. in SURFER program by Golden Software). The head of the object, except wavelength, 
			is read from the file. As file data is always of "float" type, the loss of data is very likely,
			if the GRD file is read into an XArray2D<T> object with integral type T (unless all numbers
			in the GRD file are in fact of integral type smaller or equal in size to T).
			The GRD file contains:
			(1) strings "DSAA" or "DSBB" (in ASC and BIN GRD files, respectively)
			(2) nx ny  - array dimensions as "short"s (in ASC format, separated by a white space)
			(3) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a white space)
			(4) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a white space)
			(5) zlo zhi - mnimum and maximum array value as "double"s (in ASC format, separated by a white space)
			(6) u[i][j] - array values as "float"s (in ASC format, separated by a white space), 
							j index changes most rapidly
	*/	
	template <class T> static void ReadFileGRD(XArray2D<T>& rXAr2D, const char* pchFilename, double dblWavelength, bool bForceRead = false)
	{
		if (std::numeric_limits<T>::is_integer && !bForceRead)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::ReadFileGRD (integer XArray2D)");

		char strtemp[5], buf[129], buf1[129];
		short stemp;
		index_t i, j, nx, ny;
		double dtemp;
		double wl, xlo, xhi, ylo, yhi;

		strcpy(strtemp, "ZZZZ");
		wl = dblWavelength;
		
		// open file
		try 
		{ 
			// open file for reading as BIN
			//FilePtr fp(pchFilename, "rb");
			int fd = _open(pchFilename, _O_RDONLY | _O_BINARY | _O_SEQUENTIAL, _S_IREAD); 
			if (fd == -1) 
			{
				throw std::runtime_error("cannot open file for reading in XArData::ReadFileGRD");
			}
			
			// read attribute
			//if (fread(strtemp, sizeof(char), 4, fp) != 4 || strcmp(strtemp, "DSBB"))
			strtemp[4] = '\0';
			if (_read(fd, strtemp, 4) != 4 || strcmp(strtemp, "DSBB"))
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading file attribute)");

			// read head
			//if (fread(&stemp, sizeof(short), 1, fp) - 1)
			if (_read(fd, &stemp, sizeof(short)) != sizeof(short))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading nx)");
			}
			nx = stemp;
			//if (fread(&stemp, sizeof(short), 1, fp) - 1)
			if (_read(fd, &stemp, sizeof(short)) != sizeof(short))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading ny)");
			}
			ny = stemp;
			
			rXAr2D.Resize(ny, nx);	

			//if (fread(&xlo, sizeof(double), 1, fp) - 1)
			if (_read(fd, &xlo, sizeof(double)) != sizeof(double))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading xlo)");
			}
			//if (fread(&xhi, sizeof(double), 1, fp) - 1)
			if (_read(fd, &xhi, sizeof(double)) != sizeof(double))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading xhi)");
			}
			//if (fread(&ylo, sizeof(double), 1, fp) - 1)
			if (_read(fd, &ylo, sizeof(double)) != sizeof(double))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading ylo)");
			}
			//if (fread(&yhi, sizeof(double), 1, fp) - 1)
			if (_read(fd, &yhi, sizeof(double)) != sizeof(double))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading yhi)");
			}
			//if (fread(&dtemp, sizeof(double), 1, fp) - 1)
			if (_read(fd, &dtemp, sizeof(double)) != sizeof(double))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading zlo)");
			}
			//if (fread(&dtemp, sizeof(double), 1, fp) - 1)
			if (_read(fd, &dtemp, sizeof(double)) != sizeof(double))
			{
				_close(fd);
				throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading zhi)");
			}

			// read array
			//for (i=0; i<ny; i++) 
				//for (j=0; j<nx; j++)
				//{ 
				//	fread(&ftemp, sizeof(float), 1, fp);
				//	rXAr2D[i][j] = T(ftemp);
				//}
			float* arr(0);
			if (rXAr2D.GetValuetype() != eXAFloat) // a buffer for type conversion is required
			{
				arr = new float[nx]; // the size of this buffer has to be small enough for network transmission
				if (arr == 0) throw std::runtime_error("runtime_error in XArData::ReadFileGRD (failed to allocate memory buffer)");
			}
			for (index_t i = 0; i < ny; i++)
			{
				if (rXAr2D.GetValuetype() == eXAFloat) arr = (float*)rXAr2D[i];
				if (_read(fd, arr, unsigned(nx * sizeof(float))) != int(nx * sizeof(float)))
				{
					if (rXAr2D.GetValuetype() != eXAFloat) delete [] arr;
					_close(fd);
					throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading array from file)");
				}
				if (rXAr2D.GetValuetype() != eXAFloat)
					for (index_t j = 0; j < nx; j++) 
						rXAr2D[i][j] = T(arr[j]);
			}
			if (rXAr2D.GetValuetype() != eXAFloat) delete [] arr;
			_close(fd);
		}
		catch (std::runtime_error&) // if BIN reading failed, try ASC
		{ 
			try
			{
				// open file for reading as ASC
				FilePtr fp(pchFilename, "rt");

				// read attribute
				if (fscanf(fp, "%4c\n", strtemp) <= 0 || strcmp(strtemp, "DSAA"))
					throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading file attribute)");
				
				// read head
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading nx or ny)");
				nx = strtol(buf, 0, 10); ny = strtol(buf1, 0, 10);
				rXAr2D.Resize(ny, nx);		
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading xlo or xhi)");
				xlo = strtod(buf, 0); xhi = strtod(buf1, 0);
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading ylo or yhi)");
				ylo = strtod(buf, 0); yhi = strtod(buf1, 0);
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading zlo or zhi)");
			
				// read array
				for (i=0; i<ny; i++) 
					for (j=0; j<nx; j++)
					{ 
						if (fscanf(fp, "%s", buf) <= 0)
							throw std::runtime_error("runtime_error in XArData::ReadFileGRD (error reading array)");
						rXAr2D[i][j] = T(strtod(buf, 0));
					}
			}
			catch (std::runtime_error&) { throw; } // if ASC reading failed - abort
		}

		// create and attach a new head
		IXAHWave2D* ph2 = CreateWavehead2D();
		ph2->SetData(wl, ylo, yhi, xlo, xhi);
		rXAr2D.SetHeadPtr(ph2);
	}


	//! Reports the type of a GRD file
	static _eFileType ReadFileGRDType(const char* pchFilename)
	{
		char strtemp[5];

		strcpy(strtemp, "ZZZZ");

		_eFileType filetype;
		
		//open file
		try 
		{ 
			//open file for reading as BIN
			FilePtr fp(pchFilename, "rb");
			filetype = eGRDBIN;
			
			//read attribute
			if (fread(strtemp, sizeof(char), 4, fp) != 4 || strcmp(strtemp, "DSBB"))
				throw std::runtime_error("runtime_error in XArData::ReadFileGRDType (error reading file attribute)");
		}
		catch (std::runtime_error&) //if BIN reading failed, try ASC
		{ 
			try
			{
				//open file for reading as ASC
				FilePtr fp(pchFilename, "rt");
				filetype = eGRDASC;

				//read attribute
				if (fscanf(fp, "%4c\n", strtemp) <= 0 || strcmp(strtemp, "DSAA"))
					throw std::runtime_error("runtime_error in XArData::ReadFileGRDType (error reading file attribute)");
			}
			catch (std::runtime_error&) { throw; }//if ASC reading failed - abort
		}

		return filetype;
	}


	//---------------------------------------------------------------------------
	//Function WriteFileGRC
	//
	//	Writes dcomplex XArray2D object into an ASC or BIN GRC file
	//
	/*!
		\brief		Writes dcomplex XArray2D object into an ASC or BIN GRC file
		\param		rXAr2D	Reference to an XArray2D object to be saved
		\param		pchFilename	Full name of the file to be written into
		\param		filetype eGRCBIN or eGRCASC for BIN and ASC files respectively
		\exception  std::invalid_argument is thrown if XArray2D object does not have a Wavehead2D
		\exception  std::invalid_argument is thrown if 'filetype' parameter is not eGRCBIN or eGRCASC
		\exception  std::runtime_error is thrown if any of the file write operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function writes a dcomplex 2D XArray object into a BIN or ASC GRC file. The head
			of the object is saved in the file. The array data is always stored as pairs of "float" values. 
			The GRC file contains:
			(1) -5 - this number identifies GRC format
			(2)	wl 0.0 - wavelength and 0.0 (reserved value) as "double"s (in ASC format, separated by a white space)
			(3) nx ny  - array dimensions as "short"s (in ASC format, separated by a white space)
			(4) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a white space)
			(5) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a white space)
			(6) Re(u[i][j]) Im(u[i][j]) - array values as pairs of "float" numbers (in ASC format,
										separated by a white space),  j index changes most rapidly
	*/	
	static void WriteFileGRC(const XArray2D<dcomplex>& rXAr2D, const char* pchFilename, _eFileType filetype)
	{
		if (!rXAr2D.GetHeadPtr())
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRC (no head)");
		const IXAHWave2D* ph2 = dynamic_cast<const IXAHWave2D*>(rXAr2D.GetHeadPtr());
		if (!ph2)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRC (the head is not a Wavehead2D)");
		ph2->Validate();

		index_t ny = rXAr2D.GetDim1();
		index_t nx = rXAr2D.GetDim2();
		double wl = ph2->GetWl(); 
		double xlo = ph2->GetXlo();
		double xhi = ph2->GetXhi();
		double ylo = ph2->GetYlo();
		double yhi = ph2->GetYhi();

		index_t i, j;
		short stemp;
		float ftemp;
		double dtemp;

		//open file
		FilePtr fp;
		switch (filetype)
		{
		case eGRCBIN: fp.Open(pchFilename, "wb"); break;
		case eGRCASC: fp.Open(pchFilename, "wt"); break;
		default: throw std::invalid_argument("invalid_argument 'filetype' in XArData::WriteFileGRC");
		}

		//write GRC head
		switch (filetype)
		{
		case eGRCBIN:
			stemp = -5;
			if (fwrite(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing file attribute)");
			if (fwrite(&wl, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing wl)");
			dtemp = 0.0;
			if (fwrite(&dtemp, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing theta)");
			stemp = (short)(nx);
			if (fwrite(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing nx)");
			stemp = (short)(ny);
			if (fwrite(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing ny)");
			if (fwrite(&xlo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing xlo)");
			if (fwrite(&xhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing xhi)");
			if (fwrite(&ylo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing zlo)");
			if (fwrite(&yhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing zhi)");
			break;
		case eGRCASC:
			if (fprintf(fp, "%i\n", -5) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing file attribute)");
			dtemp = 0.0;
			if (fprintf(fp, "%.15g   %.15g\n", wl, dtemp) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing wl or theta)");
			if (fprintf(fp, "%zi   %zi\n", nx, ny) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing nx or ny)");
			if (fprintf(fp, "%.15g   %.15g\n", xlo, xhi) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing xlo or xhi)");
			if (fprintf(fp, "%.15g   %.15g\n", ylo, yhi) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing ylo or yhi)");
			break;
		}

		//write array
		switch (filetype)
		{
		case eGRCBIN:
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
				{ 
					ftemp=(float)std::real(rXAr2D[i][j]);
					fwrite(&ftemp, sizeof(float), 1, fp);
					ftemp=(float)std::imag(rXAr2D[i][j]);
					fwrite(&ftemp, sizeof(float), 1, fp);
				}
				break;
		case eGRCASC:
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
					if (fprintf(fp, "%.6g   %.6g\n", (float)std::real(rXAr2D[i][j]), (float)std::imag(rXAr2D[i][j])) <= 0)
						throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing array)");
			break;
		}
	}

	//---------------------------------------------------------------------------
	//Function WriteFileGRC
	//
	//	Writes fcomplex XArray2D object into an ASC or BIN GRC file
	//
	/*!
		\brief		Writes fcomplex XArray2D object into an ASC or BIN GRC file
		\param		rXAr2D	Reference to an XArray2D object to be saved
		\param		pchFilename	Full name of the file to be written into
		\param		filetype eGRCBIN or eGRCASC for BIN and ASC files respectively
		\exception  std::invalid_argument is thrown if XArray2D object does not have a Wavehead2D
		\exception  std::invalid_argument is thrown if 'filetype' parameter is not eGRCBIN or eGRCASC
		\exception  std::runtime_error is thrown if any of the file write operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function writes an fcomplex 2D XArray object into a BIN or ASC GRC file. The head
			of the object is saved in the file. The array data is always stored as pairs of "float" values. 
			The GRC file contains:
			(1) -5 - this number identifies GRC format
			(2)	wl 0.0 - wavelength and 0.0 (reserved value) as "double"s (in ASC format, separated by a white space)
			(3) nx ny  - array dimensions as "short"s (in ASC format, separated by a white space)
			(4) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a white space)
			(5) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a white space)
			(6) Re(u[i][j]) Im(u[i][j]) - array values as pairs of "float" numbers (in ASC format,
										separated by a white space), j index changes most rapidly
	*/	
	static void WriteFileGRC(const XArray2D<fcomplex>& rXAr2D, const char* pchFilename, _eFileType filetype)
	{
		if (!rXAr2D.GetHeadPtr())
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRC (no head)");
		const IXAHWave2D* ph2 = dynamic_cast<const IXAHWave2D*>(rXAr2D.GetHeadPtr());
		if (!ph2)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in XArData::WriteFileGRC (the head is not a Wavehead2D)");
		ph2->Validate();

		index_t ny = rXAr2D.GetDim1();
		index_t nx = rXAr2D.GetDim2();
		double wl = ph2->GetWl(); 
		double xlo = ph2->GetXlo();
		double xhi = ph2->GetXhi();
		double ylo = ph2->GetYlo();
		double yhi = ph2->GetYhi();

		index_t i, j;
		short stemp;
		float ftemp;
		double dtemp;

		// open file
		FilePtr fp;
		switch (filetype)
		{
		case eGRCBIN: fp.Open(pchFilename, "wb"); break;
		case eGRCASC: fp.Open(pchFilename, "wt"); break;
		default: throw std::invalid_argument("invalid_argument 'filetype' in XArData::WriteFileGRC");
		}

		// write GRC head
		switch (filetype)
		{
		case eGRCBIN:
			stemp = -5;
			if (fwrite(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing file attribute)");
			if (fwrite(&wl, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing wl)");
			dtemp = 0.0;
			if (fwrite(&dtemp, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing theta)");
			stemp = (short)(nx);
			if (fwrite(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing nx)");
			stemp = (short)(ny);
			if (fwrite(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing ny)");
			if (fwrite(&xlo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing xlo)");
			if (fwrite(&xhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing xhi)");
			if (fwrite(&ylo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing zlo)");
			if (fwrite(&yhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing zhi)");
			break;
		case eGRCASC:
			if (fprintf(fp, "%i\n", -5) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing file attribute)");
			dtemp = 0.0;
			if (fprintf(fp, "%.15g   %.15g\n", wl, dtemp) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing wl or theta)");
			if (fprintf(fp, "%zi   %zi\n", nx, ny) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing nx or ny)");
			if (fprintf(fp, "%.15g   %.15g\n", xlo, xhi) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing xlo or xhi)");
			if (fprintf(fp, "%.15g   %.15g\n", ylo, yhi) <= 0)
				throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing ylo or yhi)");
			break;
		}

		// write array
		switch (filetype)
		{
		case eGRCBIN:
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
				{ 
					ftemp=(float)std::real(rXAr2D[i][j]);
					fwrite(&ftemp, sizeof(float), 1, fp);
					ftemp=(float)std::imag(rXAr2D[i][j]);
					fwrite(&ftemp, sizeof(float), 1, fp);
				}
				break;
		case eGRCASC:
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
					if (fprintf(fp, "%.6g   %.6g\n", std::real(rXAr2D[i][j]), std::imag(rXAr2D[i][j])) <= 0)
						throw std::runtime_error("runtime_error in XArData::WriteFileGRC (error writing array)");
			break;
		}
	}

	//---------------------------------------------------------------------------
	//Function ReadFileGRC
	//
	//	Reads dcomplex XArray2D object from an ASC or BIN GRC file
	//
	/*!
		\brief		Reads dcomplex XArray2D object from an ASC or BIN GRC file
		\param		rXAr2D	Reference to an XArray2D object to be read in
		\param		pchFilename	Full name of the file to be read from
		\exception  std::runtime_error is thrown if any of the file read operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function reads a dcomplex 2D XArray object from a BIN or ASC GRC file. The Wavehead2D
			head of the object is read from the file.
			The GRC file contains:
			(1) -5 - this number identifies GRC format
			(2)	wl 0.0 - wavelength and 0.0 (reserved value) as "double"s (in ASC format, separated by a white space)
			(3) nx ny  - array dimensions as "short"s (in ASC format, separated by a white space)
			(4) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a white space)
			(5) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a white space)
			(6) Re(u[i][j]) Im(u[i][j]) - array values as pairs of "float" numbers (in ASC format,
									separated by a white space), j index changes most rapidly
	*/	
	static void ReadFileGRC(XArray2D<dcomplex>& rXAr2D, const char* pchFilename)
	{
		char buf[129], buf1[129];
		short stemp;
		index_t i, j, nx, ny;
		float ftemp, ftemp1;
		double dtemp;
		double wl, xlo, xhi, ylo, yhi;

		// open file
		try 
		{ 
			// open file for reading as BIN
			FilePtr fp(pchFilename, "rb");
			
			// read attribute
			if (fread(&stemp, sizeof(short), 1, fp) != 1 || stemp != -5)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file attribute)");

			// read head
			if (fread(&wl, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading wavelength)");
			if (fread(&dtemp, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading theta)");
			if (fread(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading nx)");
			nx = stemp;
			if (fread(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading ny)");
			ny = stemp;
			rXAr2D.Resize(ny, nx);		
			if (fread(&xlo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading xlo)");
			if (fread(&xhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading xhi)");
			if (fread(&ylo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading ylo)");
			if (fread(&yhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading yhi)");

			// read array
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
				{ 
					fread(&ftemp, sizeof(float), 1, fp);
					fread(&ftemp1, sizeof(float), 1, fp);
					rXAr2D[i][j] = dcomplex(ftemp, ftemp1);
				}
		}
		catch (std::runtime_error&) // if BIN reading failed, try ASC
		{ 
			try
			{
				// open file for reading as ASC
				FilePtr fp(pchFilename, "rt");

				// read attribute
				if (fscanf(fp, "%s", buf) <= 0 || strtol(buf, 0, 10) != -5)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file attribute)");
				
				// read head
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file wavelength or theta)");
				wl = strtod(buf, 0);
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file nx or ny)");
				nx = strtol(buf, 0, 10); ny = strtol(buf1, 0, 10);
				rXAr2D.Resize(ny, nx);		
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file xlo or xhi)");
				xlo = strtod(buf, 0); xhi = strtod(buf1, 0);
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file ylo or yhi)");
				ylo = strtod(buf, 0); yhi = strtod(buf1, 0);
			
				// read array
				for (i=0; i<ny; i++) 
					for (j=0; j<nx; j++)
					{ 
						if (fscanf(fp, "%s %s", buf, buf1) <= 0)
							throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading array)");
						rXAr2D[i][j] = dcomplex(strtod(buf, 0), strtod(buf1, 0));
					}
			}
			catch (std::runtime_error&) { throw; }//if ASC reading failed - abort
		}

		// create and attach a new head
		IXAHWave2D* ph2 = CreateWavehead2D();
		ph2->SetData(wl, ylo, yhi, xlo, xhi);
		rXAr2D.SetHeadPtr(ph2);
	}

	//---------------------------------------------------------------------------
	//Function ReadFileGRC
	//
	//	Reads fcomplex XArray2D object from an ASC or BIN GRC file
	//
	/*!
		\brief		Reads fcomplex XArray2D object from an ASC or BIN GRC file
		\param		rXAr2D	Reference to an XArray2D object to be read in
		\param		pchFilename	Full name of the file to be read from
		\exception  std::runtime_error is thrown if any of the file read operations fail
		\exception	std::exception and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function reads an fcomplex 2D XArray object from a BIN or ASC GRC file. The Wavehead2D
			head of the object is read from the file.
			The GRC file contains:
			(1) -5 - this number identifies GRC format
			(2)	wl 0.0 - wavelength and 0.0 (reserved value) as "double"s (in ASC format, separated by a white space)
			(3) nx ny  - array dimensions as "short"s (in ASC format, separated by a white space)
			(4) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a white space)
			(5) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a white space)
			(6) Re(u[i][j]) Im(u[i][j]) - array values as pairs of "float" numbers (in ASC format,
							separated by a white space), j index changes most rapidly
	*/	
	static void ReadFileGRC(XArray2D<fcomplex>& rXAr2D, const char* pchFilename)
	{
		char buf[129], buf1[129];
		short stemp;
		index_t i, j, nx, ny;
		float ftemp, ftemp1;
		double dtemp;
		double wl, xlo, xhi, ylo, yhi;

		// open file
		try 
		{ 
			// open file for reading as BIN
			FilePtr fp(pchFilename, "rb");
			
			// read attribute
			if (fread(&stemp, sizeof(short), 1, fp) != 1 || stemp != -5)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file attribute)");

			// read head
			if (fread(&wl, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading wavelength)");
			if (fread(&dtemp, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading theta)");
			if (fread(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading nx)");
			nx = stemp;
			if (fread(&stemp, sizeof(short), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading ny)");
			ny = stemp;
			rXAr2D.Resize(ny, nx);		
			if (fread(&xlo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading xlo)");
			if (fread(&xhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading xhi)");
			if (fread(&ylo, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading ylo)");
			if (fread(&yhi, sizeof(double), 1, fp) - 1)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading yhi)");

			// read array
			for (i=0; i<ny; i++) 
				for (j=0; j<nx; j++)
				{ 
					fread(&ftemp, sizeof(float), 1, fp);
					fread(&ftemp1, sizeof(float), 1, fp);
					rXAr2D[i][j] = fcomplex(ftemp, ftemp1);
				}
		}
		catch (std::runtime_error&) // if BIN reading failed, try ASC
		{ 
			try
			{
				// open file for reading as ASC
				FilePtr fp(pchFilename, "rt");

				// read attribute
				if (fscanf(fp, "%s", buf) <= 0 || strtol(buf, 0, 10) != -5)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file attribute)");
				
				// read head
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file wavelength or theta)");
				wl = strtod(buf, 0);
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file nx or ny)");
				nx = strtol(buf, 0, 10); ny = strtol(buf1, 0, 10);
				rXAr2D.Resize(ny, nx);		
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file xlo or xhi)");
				xlo = strtod(buf, 0); xhi = strtod(buf1, 0);
				if (fscanf(fp, "%s %s", buf, buf1) <= 0)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading file ylo or yhi)");
				ylo = strtod(buf, 0); yhi = strtod(buf1, 0);
			
				// read array
				for (i=0; i<ny; i++) 
					for (j=0; j<nx; j++)
					{ 
						if (fscanf(fp, "%s %s", buf, buf1) <= 0)
							throw std::runtime_error("runtime_error in XArData::ReadFileGRC (error reading array)");
						rXAr2D[i][j] = fcomplex(float(strtod(buf, 0)), float(strtod(buf1, 0)));
					}
			}
			catch (std::runtime_error&) { throw; } // if ASC reading failed - abort
		}

		// create and attach a new head
		IXAHWave2D* ph2 = CreateWavehead2D();
		ph2->SetData(wl, ylo, yhi, xlo, xhi);
		rXAr2D.SetHeadPtr(ph2);
	}


	//! Reads a complex Wave2D object type from an ASC or BIN GRC file. 
	static _eFileType ReadFileGRCType(const char* pchFilename)
	{
		short stemp;
		char buf[129];

		_eFileType filetype;
		
		// open file
		try 
		{ 
			// open file for reading as BIN

			FilePtr fp(pchFilename, "rb");
			filetype = eGRCBIN;
			
			// read attribute
			if (fread(&stemp, sizeof(short), 1, fp) != 1 || stemp != -5)
				throw std::runtime_error("runtime_error in XArData::ReadFileGRCType (error reading file attribute)");
		}
		catch (std::runtime_error&) // if BIN reading failed, try ASC
		{ 
			try
			{
				// open file for reading as ASC
				FilePtr fp(pchFilename, "rt");
				filetype = eGRCASC;

				// read attribute
				if (fscanf(fp, "%s", buf) <= 0 || strtol(buf, 0, 10) != -5)
					throw std::runtime_error("runtime_error in XArData::ReadFileGRCType (error reading file attribute)");
			}
			catch (std::runtime_error&) { throw; } // if ASC reading failed - abort
		}

		return filetype;
	}

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
	// Member functions

	};


} // namespace xar closed

//---------------------------------------------------------------------------
//	INLINE AND TEMPLATE MEMBER DEFINITIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//


#endif  // XA_DATA_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
