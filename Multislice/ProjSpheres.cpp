#include "XArray2D.h"
#include "XArray3D.h"

using namespace xar;

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
	
/*
		nk = R.shape[0]
		if nk != xr.shape[0] or nk != yr.shape[0]:
	raise ValueError("some of input arrays R, xr, yr have different sizes in ProjectSphereN1()")

		nx = camp.shape[0]
		ny = camp.shape[1]
		xst = (self.xhi - self.xlo) / nx
		yst = (self.yhi - self.ylo) / ny
		R2 = R * *2;
*/
	XArray3D<double>* pRo = new XArray3D<double>(1,1,1); // ((nx, ny, nk))
/*
		for i in range(0, nx) :
			x = self.xlo + xst * i
			for j in range(0, ny) :
				y = self.ylo + yst * j
				for k in range(0, nk) :
					xk = x - xr[k]
					yk = y - yr[k]
					r2 = R2[k] - xk * xk - yk * yk
					if r2 > 0 :
	ro[i, j, k] = sqrt(r2)
*/
	return pRo;
}

/*
def ProjectSphereN2(self, camp, ro, nc, zr, z0, z1) :
	"""Calculates projection (within a z-slice) of a complex amplitude through a set of solid spheres, given a 3D array of projection 'chords'
	This function should be called for each slice in multi - slice calculations, after a single initial call to ProjectSphereN1()
	Parameters :
	camp - 2D complex amplitude array
	ro - 3D array of projection chord lengths as produced by ProjectionSphere1A() function
	nc - 1D array of complex refraction indexes, 1 - delta + i * beta, of the uniform spheres
	zr - 1D array of z - coordinates of the centres of the spheres
	z0 - minimal coordinate of the z - slice
	z1 - maximal coordinate of the z - slice
	returns - input camp multipied by projection integrals of the complex refractive index(as produced by the set of spheres)
	"""
	nx = camp.shape[0]
	ny = camp.shape[1]
	nk = nc.shape[0]

	if (ro.shape != (nx, ny, nk)) :
		raise ValueError("input array ro has incorrect dimensions in ProjectSphereN2()")

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
							return camp
*/