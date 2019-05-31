#include "XArray2D.h"
#include "XArray3D.h"
#include "IXAHWave.h"
#include "XA_head2.h"
#include "XA_data.h"

using namespace xar;


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
	/*
	nk = R.shape[0]
	if nk != nc.shape[0] or nk != xr.shape[0] or nk != yr.shape[0] or nk != zr.shape[0]:
	raise ValueError("some of input arrays R, nc, xr, yr, zr have different sizes in MultisliceSphereN()")

	nctemp0 = numpy.real(1 - nc)
	nctemp1 = numpy.imag(nc)
	if numpy.amin(R) <= 0 or numpy.amin(nctemp0) < 0 or numpy.amin(nctemp1) < 0 or nslices < 1 or zhi <= zlo or zlo > numpy.amin(zr - R) or zhi < numpy.amax(zr + R) :
		raise ValueError("incorrect input parameter(s) in MultisliceSphereN()")

		sliceT = (zhi - zlo) / nslices # slice thickness

		# calculate the array of projection chord lengths(required only if ProjectSphereN2 is called for each slice below)
		ro = self.ProjectSphereN1(camp, R, xr, yr)

		# iterate over all slices
		for m in range(0, nslices) :
			print("process = ", pID, "slice = ", m + 1)
			z0 = zlo + sliceT * m # beginning of the slice
			z1 = z0 + sliceT # end of the slice
			# calculate the projection of the complex amplitude through the next slice
			#camp = self.ProjectSphereN(camp, nc, R, xr, yr, zr, z0, z1);
	camp = self.ProjectSphereN2(camp, ro, nc, zr, z0, z1)
	# calculate free - space propagation through the next slice
	camp = self.NearFresnel(camp, z1 - z0, False)
	*/
}


void MultisliceSphereNF(int nx, int ny, vector<double> vHead, vector<dcomplex> nc, vector<double> R, vector<double> xr, vector<double> yr,
	vector<double> zr, double zlo, double zhi, int nslices, const char* outfilename)
{
	//Same as MultisliceSphereN, but writes the output directly to a GRDC file"""
		
	XArray2D<dcomplex>& camp(*new XArray2D<dcomplex>(ny, nx, dcomplex(1,0)));
	IXAHWave2D* ph2new = CreateWavehead2D();
	if (vHead.size() != 5)
		throw std::invalid_argument("incorrect size of vHead vector in MultisliceSphereNF()");
	ph2new->SetData(vHead[0], vHead[1], vHead[2], vHead[3], vHead[4]);
	camp.SetHeadPtr(ph2new);

	MultisliceSphereN(camp, nc, R, xr, yr, zr, zlo, zhi, nslices);

	XArray2D<double> inten;
	Abs2(camp, inten);
	XArData::WriteFileGRD(inten, outfilename, eGRDBIN);
}


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
	

	int nk = R.size();
	if (nk != xr.size() || nk != yr.size())
		throw std::invalid_argument("some of input arrays R, xr, yr have different sizes in ProjectSphereN1()");

	int nx = camp.GetDim2();
	int ny = camp.GetDim1();
	IXAHWave2D* ph2 = GetIXAHWave2D(camp);
	if (!ph2)
		throw std::invalid_argument("input array camp does not have a Wavehead2D header in ProjectSphereN1()");
	ph2->Validate();

	double xst = ph2->GetXStep(nx);
	double yst = ph2->GetYStep(ny);
	double xlo = ph2->GetXlo();
	double ylo = ph2->GetYlo();
	
	vector<double> R2(R.size());
	for (int i = 0; i < R.size(); i++) R2[i] = R[i] * R[i];
		
	XArray3D<double>& ro(*new XArray3D<double>(ny, nx, nk)); // ((nx, ny, nk))
	double x, y, xk, yk, r2;

	for (int i = 0; i < ny; i++)
	{
		y = ylo + yst * i;
		for (int j = 0; j < nx; j++)
		{
			x = xlo + xst * j;
			for (int k = 0; k < nk; k++)
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

	int nx = camp.GetDim2();
	int ny = camp.GetDim1();
	int nk = nc.size();

	if (ro.GetDim1() != nk || ro.GetDim2() != ny || ro.GetDim3() != nx)
		throw std::invalid_argument("input array ro has incorrect dimensions in ProjectSphereN2()");

	IXAHWave2D* ph2 = GetIXAHWave2D(camp);
	if (!ph2)
		throw std::invalid_argument("input array camp does not have a Wavehead2D header in ProjectSphereN2()");
	ph2->Validate();

/*
	dcomplex df = tPI *  camp.GetHeadPtr()->
		ifi = 1j * self.tpiwl * (nc - 1); #   - i * 2 * pi / wl * delta - 2 * pi / wl * beta

		for i in range(0, nx) :
			for j in range(0, ny) :
				cc = complex(0, 0)
				for k in range(0, nk) :
					roijk = ro[i, j, k]
					if roijk :
						zk = zr[k]
						if z0 < zk + roijk and z1 > zk - roijk :
							z0ro = max(zk - roijk, z0) # left end of the integration interval
							z1ro = min(zk + roijk, z1) # right end of the integration interval
							cc += ifi[k] * (z1ro - z0ro) # accumulate all absorption and phase shifts along this ray
							camp[i, j] *= exp(cc)

	*/
}