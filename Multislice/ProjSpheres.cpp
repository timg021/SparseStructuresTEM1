#include "ProjSpheres.h"

namespace xar
{

	XArray3D<double>& ProjectSphereN1(XArray2D<dcomplex>& camp, vector<double> R, vector<double> xr, vector<double> yr)
	{
		//Calculates a 3D array of lenghts of projection 'chords' through a set of solid spheres indexed by k for rays (along z) indexed by (i,j).
		//	This function should be called first in multi - slice calculations, with the output array ro then used for each slice.
		//	Parameters :
		//	camp - input complex amplitude array - only used here for extracting the dimensions
		//	R - 1D array of radiuses of the spheres
		//	xr - 1D array of x - coordinates of the centres of the spheres
		//	yr - 1D array of y - coordinates of the centres of the spheres
		//	returns - 3D array of projection chord lengths


		size_t nk = R.size();
		if (nk != xr.size() || nk != yr.size())
			throw std::invalid_argument("some of input arrays R, xr, yr have different sizes in ProjectSphereN1()");

		size_t nx = camp.GetDim2();
		size_t ny = camp.GetDim1();
		IXAHWave2D* ph2 = GetIXAHWave2D(camp);
		if (!ph2)
			throw std::invalid_argument("input array camp does not have a Wavehead2D header in ProjectSphereN1()");
		ph2->Validate();

		double xst = ph2->GetXStep(nx);
		double yst = ph2->GetYStep(ny);
		double xlo = ph2->GetXlo();
		double ylo = ph2->GetYlo();

		vector<double> R2(R.size());
		for (size_t i = 0; i < R.size(); i++) R2[i] = R[i] * R[i];

		XArray3D<double>& ro(*new XArray3D<double>(ny, nx, nk)); // ((nx, ny, nk))
		double x, y, xk, yk, r2;

		for (size_t i = 0; i < ny; i++)
		{
			y = ylo + yst * i;
			for (size_t j = 0; j < nx; j++)
			{
				x = xlo + xst * j;
				for (size_t k = 0; k < nk; k++)
				{
					xk = x - xr[k];
					yk = y - yr[k];
					r2 = R2[k] - xk * xk - yk * yk;
					if (r2 > 0) ro[i][j][k] = sqrt(r2);
				}
			}
		}

		return ro;
	}


	void ProjectSphereN2(XArray2D<dcomplex>& camp, XArray3D<double>& ro, vector<dcomplex> nc, vector<double> zr, double zlo, double zhi)
	{
		//Calculates projection (within a z-slice) of a complex amplitude through a set of solid spheres, given a 3D array of projection 'chords'
		//This function should be called for each slice in multi - slice calculations, after a single initial call to ProjectSphereN1()
		//Parameters :
		//camp - 2D complex amplitude array
		//ro - 3D array of projection chord lengths as produced by ProjectionSphere1A() function
		//nc - 1D array of complex refraction indexes, 1 - delta + i * beta, of the uniform spheres
		//zr - 1D array of z - coordinates of the centres of the spheres
		//z0 - minimal coordinate of the z - slice
		//z1 - maximal coordinate of the z - slice

		size_t nx = camp.GetDim2();
		size_t ny = camp.GetDim1();
		size_t nk = nc.size();

		if (ro.GetDim1() != ny || ro.GetDim2() != nx || ro.GetDim3() != nk)
			throw std::invalid_argument("input array ro has incorrect dimensions in ProjectSphereN2()");

		IXAHWave2D* ph2 = GetIXAHWave2D(camp);
		if (!ph2)
			throw std::invalid_argument("input array camp does not have a Wavehead2D header in ProjectSphereN2()");
		ph2->Validate();

		double roijk, zk, z0ro, z1ro;
		dcomplex ff = dcomplex(0.0, 1.0) * tPI / ph2->GetWl(), cc;
		for (size_t i = 0; i < nc.size(); i++) nc[i] = (nc[i] - 1.0) * ff; // - i * 2 * pi / wl * delta - 2 * pi / wl * beta

		for (size_t i = 0; i < ny; i++)
			for (size_t j = 0; j < nx; j++)
			{
				cc = dcomplex(0, 0);
				for (size_t k = 0; k < nk; k++)
				{
					roijk = ro[i][j][k];
					if (roijk)
					{
						zk = zr[k];
						if (zlo < zk + roijk && zhi > zk - roijk)
						{
							z0ro = max(zk - roijk, zlo); // left end of the integration interval
							z1ro = min(zk + roijk, zhi); // right end of the integration interval
							cc += nc[k] * (z1ro - z0ro); //accumulate all absorptionand phase shifts along this ray
							camp[i][j] *= exp(cc);
						}
					}
				}
			}
	}


	void MultisliceSphereN(XArray2D<dcomplex>& camp, vector<dcomplex> nc, vector<double> R, vector<double> xr, vector<double> yr,
		vector<double> zr, double zlo, double zhi, int nslices)
	{
		//Calculates multislice approximation for paraxial light propagation through a set of spheres
		//This function repeadetly calls a pair of functions for calculation of the projection approximationand for the Fresnel propagation through each slice
		//Parameters :
		//camp - 2D complex amplitude array
		//nc - 1D array of complex refraction indexes, 1 - delta + i * beta, of the uniform spheres
		//R - 1D array of radiuses of the spheres
		//xr - 1D array of x - coordinates of the centres of the spheres
		//yr - 1D array of y - coordinates of the centres of the spheres
		//zr - 1D array of z - coordinates of the centres of the spheres
		//zlo - minimal coordinate of the z - slices
		//zhi - maximal coordinate of the z - slices
		//nslices - number of z - slices

		//check the validity of function parameters
		if (!camp.GetHeadPtr())
			throw std::invalid_argument("camp does not have a header in MultisliceSphereN()");
		size_t nk = R.size();
		if (nk != nc.size() || nk != xr.size() || nk != yr.size() || nk != zr.size())
			throw std::invalid_argument("some of input arrays R, nc, xr, yr, zr have different sizes in MultisliceSphereN()");
		for (size_t i = 0; i < nc.size(); i++)
			if (1 - nc[i].real() < 0 || nc[i].imag() < 0)
				throw std::invalid_argument("some of nc[i] have negative delta or beta in MultisliceSphereN()");
		for (size_t i = 0; i < nc.size(); i++)
			if (R[i] <= 0)
				throw std::invalid_argument("some of R[i] are not positive in MultisliceSphereN()");
		if (nslices < 1)
			throw std::invalid_argument("number of slices is less than 1 in MultisliceSphereN()");
		if (zhi <= zlo)
			throw std::invalid_argument("zhi <= zlo in MultisliceSphereN()");
		for (size_t i = 0; i < nc.size(); i++)
			if (zlo > zr[i] - R[i] || zhi < zr[i] + R[i])
				throw std::invalid_argument("some spheres do not fit into [zlo, zhi] in MultisliceSphereN()");

		double sliceT = (zhi - zlo) / nslices; //slice thickness

		//calculate the array of projection chord lengths(required only if ProjectSphereN2 is called for each slice below)
		XArray3D<double> ro = ProjectSphereN1(camp, R, xr, yr);

		//iterate over all slices
		double z0, z1;
		XArray2DFFT<double> FS(camp);
		for (int m = 0; m < nslices; m++)
		{
			printf("\nslice = %d", m + 1);
			z0 = zlo + sliceT * m; //beginning of the slice
			z1 = z0 + sliceT; //end of the slice
			//calculate the projection of the complex amplitude through the next slice
			ProjectSphereN2(camp, ro, nc, zr, z0, z1);
			//calculate free - space propagation through the next slice
			FS.Fresnel(z1 - z0, true);
		}

	}


	void MultisliceSphereNF(int nx, int ny, vector<double> vHead, vector<dcomplex> nc, vector<double> R, vector<double> xr, vector<double> yr,
		vector<double> zr, double zlo, double zhi, int nslices, const char* outfilename)
	{
		//Same as MultisliceSphereN, but writes the output directly to a GRDC file"""

		XArray2D<dcomplex>& camp(*new XArray2D<dcomplex>(ny, nx, dcomplex(1, 0)));
		IXAHWave2D* ph2new = CreateWavehead2D();
		if (vHead.size() != 5)
			throw std::invalid_argument("incorrect size of vHead vector in MultisliceSphereNF()");
		ph2new->SetData(vHead[0], vHead[1], vHead[2], vHead[3], vHead[4]);
		camp.SetHeadPtr(ph2new);

		MultisliceSphereN(camp, nc, R, xr, yr, zr, zlo, zhi, nslices);

		XArray2D<double> inten;
		Abs2(camp, inten);
		XArData::WriteFileGRD(inten, outfilename, eGRDBIN);
		printf("\nOutput file = %s", outfilename);
	}


}