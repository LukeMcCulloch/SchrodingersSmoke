// This file collects functions for calculus on a 3-torus discretized
// into a grid of resx * resy * resz many rectangular voxels of side
// length dx, dy, dz.
//
// Each cell (edge, square or voxel) of the grid is indexed by
// the index (ix,iy,iz) of the vertex where the function x+y+z,
// when restricted to the cell, takes on its minimum.
//
// We only deal with primal forms, so there are four kinds of objects:
// functions, 1-forms, 2-forms and 3-forms. For operations involving the
// Hodge star that would result in a dual form, we always apply another
// Hodge star to the result in order to remain in the realm of primal forms.

#include "math.h";

// For a function f compute the 1-form df
void DerivativeOfFunction( export float vx, vy, vz; const string f;
	const int ix, iy, iz, resx, resy, resz )
{
	int ixp = (ix+1)%resx;
	int iyp = (iy+1)%resy;
	int izp = (iz+1)%resz;
	float f0 = volumeindex(0,f,set(ix,iy,iz));
	
	vx = volumeindex(0,f,set(ixp,iy,iz)) - f0;
	vy = volumeindex(0,f,set(ix,iyp,iz)) - f0;
	vz = volumeindex(0,f,set(ix,iy,izp)) - f0;
}

// For a 1-form w compute the 2-form dw
void DerivativeOfOneForm( export float ux, uy, uz; const string wx, wy, wz;
	const int ix, iy, iz, resx, resy, resz )
{
	int ixp = (ix+1)%resx;
	int iyp = (iy+1)%resy;
	int izp = (iz+1)%resz;
	vector voxel = set(ix,iy,iz);
	float vx = volumeindex(0,wx,voxel);
	float vy = volumeindex(0,wy,voxel);
	float vz = volumeindex(0,wz,voxel);
	
	float cx  = vy - volumeindex(0,wy,set(ix,iy,izp));
	cx += volumeindex(0,wz,set(ix,iyp,iz)) - vz;
	
	float cy  = vz - volumeindex(0,wz,set(ixp,iy,iz));
	cy += volumeindex(0,wx,set(ix,iy,izp)) - vx;

	float cz  = vx - volumeindex(0,wx,set(ix,iyp,iz));
	cz += volumeindex(0,wy,set(ixp,iy,iz)) - vy;
	
	ux = cx;
	uy = cy;
	uz = cz;
}

// For a 2-form w compute the 3-form dw
void DerivativeOfTwoForm( export float f; const string wx, wy, wz;
	const int ix, iy, iz, resx, resy, resz )
{
	int ixp = (ix+1)%resx;
	int iyp = (iy+1)%resy;
	int izp = (iz+1)%resz;
	vector voxel = set(ix,iy,iz);
	
	f  = (volumeindex(0,wx,set(ixp,iy,iz)) - volumeindex(0,wx,voxel));
	f += (volumeindex(0,wy,set(ix,iyp,iz)) - volumeindex(0,wy,voxel));
	f += (volumeindex(0,wz,set(ix,iy,izp)) - volumeindex(0,wz,voxel));
}

// For a 1-form w compute the function *d*w
void Div( export float f; const string wx, wy, wz;
	const int ix, iy, iz, resx, resy, resz; const float dx, dy, dz )
{
	int ixm = (ix-1)%resx;
	int iym = (iy-1)%resy;
	int izm = (iz-1)%resz;
	vector voxel = set(ix,iy,iz);
	
	f  = (volumeindex(0,wx,voxel) - volumeindex(0,wx,set(ixm,iy,iz))) / (dx*dx);
	f += (volumeindex(0,wy,voxel) - volumeindex(0,wy,set(ix,iym,iz))) / (dy*dy);
	f += (volumeindex(0,wz,voxel) - volumeindex(0,wz,set(ix,iy,izm))) / (dz*dz);
}
/*
// For a 2-form w compute the 1-form *d*w
void Curl( export float cx, cy, cz; const string wx, wy, wz;
	const int ix, iy, iz, resx, resy, resz; const float dx, dy, dz )
{
	
}

// For a 3-form w compute the 2-form *d*w
void Grad( export float vx, vy, vz; const string w;
	const int ix, iy, iz, resx, resy, resz )
{
	
}
*/
// For 1-forms v,w compute the function *(v wedge *w). Note that
// v wedge *w first assigns a volume integral to each edge. We use half
// of this integral as a contribution to the volume integral over the
// dual cell around each of the two vertices adjacent to the edge.
// Finally, we use a Hodge star to convert the resulting dual 3-form to
// a primal function.
void InnerProductOfOneForms( export float f;  const string vx,vy,vz,wx, wy, wz;
	const int ix, iy, iz, resx, resy, resz; const float dx, dy, dz )
{
	int ixm = (ix-1) % resx;
	int iym = (iy-1) % resy;
	int izm = (iz-1) % resz;
	vector voxel = set(ix,iy,iz);
	
	f  = volumeindex(0,vx,voxel) * volumeindex(0,wx,voxel) / (dx*dx);
	f += volumeindex(0,vy,voxel) * volumeindex(0,wy,voxel) / (dy*dy);
	f += volumeindex(0,vz,voxel) * volumeindex(0,wz,voxel) / (dz*dz);
	voxel = set(ixm,iy,iz);
	f += volumeindex(0,vx,voxel) * volumeindex(0,wx,voxel) / (dx*dx);
	voxel = set(ix,iym,iz);
	f += volumeindex(0,vy,voxel) * volumeindex(0,wy,voxel) / (dy*dy);
	voxel = set(ix,iy,izm);
	f += volumeindex(0,vz,voxel) * volumeindex(0,wz,voxel) / (dz*dz);
}

// For 2-forms v,w compute the 3-form v wedge *w. Note that
// v wedge *w first assigns a volume integral to each face. We use half
// of this integral as a contribution to the volume integral over each
// of the two voxels adjacent to the face.

void InnerProductOfTwoForms( export float f; const string vx, vy, vz, wx, wy, wz;
	const int ix, iy, iz, resx, resy, resz; const float dx, dy, dz )
{
	int ixp = (ix+1) % resx;
	int iyp = (iy+1) % resy;
	int izp = (iz+1) % resz;
	vector voxel = set(ix,iy,iz);

	f  = volumeindex(0,vx,voxel) * volumeindex(0,wx,voxel) * dy*dz/dx;
	f += volumeindex(0,vy,voxel) * volumeindex(0,wy,voxel) * dz*dx/dy;
	f += volumeindex(0,vz,voxel) * volumeindex(0,wz,voxel) * dx*dy/dz;
	voxel = set(ixp,iy,iz);
	f += volumeindex(0,vx,voxel) * volumeindex(0,wx,voxel) * dy*dz/dx;
	voxel = set(ix,iyp,iz);
	f += volumeindex(0,vy,voxel) * volumeindex(0,wy,voxel) * dz*dx/dy;
	voxel = set(ix,iy,izp);
	f += volumeindex(0,vz,voxel) * volumeindex(0,wz,voxel) * dx*dy/dz;
	f *= .5;
}

// For a 1-form w compute the corresponding vector field w^sharp
// by averaging to vertices
void Sharp( export float ux, uy, uz; const string vx, vy, vz;
	const int ix, iy, iz, resx, resy, resz; const float dx, dy, dz )
{
	int ixm = (ix-1) % resx;
	int iym = (iy-1) % resy;
	int izm = (iz-1) % resz;
	vector voxel = set(ix,iy,iz);
	
	ux = .5*(volumeindex(0,vx,set(ixm,iy,iz)) + volumeindex(0,vx,voxel))/dx;
	uy = .5*(volumeindex(0,vy,set(ix,iym,iz)) + volumeindex(0,vy,voxel))/dy;
	uz = .5*(volumeindex(0,vz,set(ix,iy,izm)) + volumeindex(0,vz,voxel))/dz;
}

// For a 1-form w compute the corresponding vector field w^sharp
// as a staggered vector field living on edges
void StaggeredSharp( export float ux, uy, uz; const string vx, vy, vz;
	const int ix, iy, iz; const float dx, dy, dz )
{
	vector voxel = set(ix,iy,iz);
	
	ux = volumeindex(0,vx,voxel)/dx;
	uy = volumeindex(0,vy,voxel)/dy;
	uz = volumeindex(0,vz,voxel)/dz;
}
	

// Poisson solve in Fourier domain
void PoissonSolve( export float fre, fim; 
    const int ix, iy, iz, resx, resy, resz; const float dx, dy, dz )
{
	float denom = 0.0;
	float s = sin(PI*ix/resx)/dx;
	denom += s*s;
	s = sin(PI*iy/resy)/dy;
	denom += s*s;
	s = sin(PI*iz/resz)/dz;
	denom += s*s;
	float fac = 0.0;
	if (denom > 1e-16) {
	  fac = -.25/denom;
	}
	fre *= fac;
	fim *= fac;
}
