#include "math.h";
#include "Util.h";

// Compute the velocity for a normalized two-component wave function
void VelocityOneForm( const string psi1re, psi1im, psi2re, psi2im;
    const int ix, iy, iz, resx, resy, resz ; 
	const float hbar, dx, dy, dz;
    export float wx, wy, wz )
{
	vector ip = set(ix,iy,iz);
	float p1r = volumeindex(0,psi1re,ip);
	float p1i = volumeindex(0,psi1im,ip);
	float p2r = volumeindex(0,psi2re,ip);
	float p2i = volumeindex(0,psi2im,ip);
	
	int ixp = (ix+1) % resx;
	int iyp = (iy+1) % resy;
	int izp = (iz+1) % resz;
	
	// x increment
	vector iq = set(ixp,iy,iz);
	float q1r = volumeindex(0,psi1re,iq);
	float q1i = volumeindex(0,psi1im,iq);
	float q2r = volumeindex(0,psi2re,iq);
	float q2i = volumeindex(0,psi2im,iq);
	float a = p1r*q1r + p1i*q1i + p2r*q2r + p2i*q2i;
	float b = p1r*q1i - p1i*q1r + p2r*q2i - p2i*q2r;
	wx = hbar * atan2(b,a);
	
	// y increment
	iq = set(ix,iyp,iz);
	q1r = volumeindex(0,psi1re,iq);
	q1i = volumeindex(0,psi1im,iq);
	q2r = volumeindex(0,psi2re,iq);
	q2i = volumeindex(0,psi2im,iq);
	a = p1r*q1r + p1i*q1i + p2r*q2r + p2i*q2i;
	b = p1r*q1i - p1i*q1r + p2r*q2i - p2i*q2r;
	wy = hbar * atan2(b,a);
	
	// z increment
	iq = set(ix,iy,izp);
	q1r = volumeindex(0,psi1re,iq);
	q1i = volumeindex(0,psi1im,iq);
	q2r = volumeindex(0,psi2re,iq);
	q2i = volumeindex(0,psi2im,iq);
	a = p1r*q1r + p1i*q1i + p2r*q2r + p2i*q2i;
	b = p1r*q1i - p1i*q1r + p2r*q2i - p2i*q2r;
	wz = hbar * atan2(b,a);
}


// Compute the velocity of a one-component wave function
void VelocityOneForm( const string re, im;
	const int ix, iy, iz, resx, resy, resz ; 
	const float hbar, dx, dy, dz;
	export float wx, wy, wz )
{
	vector ip = set(ix,iy,iz);
	float pr = volumeindex(0,re,ip);
	float pi = volumeindex(0,im,ip);
	if (pr<1e-6 && pi<1e-6) pr=1e-6;
	
	int ixp = (ix+1) % resx;
	int iyp = (iy+1) % resy;
	int izp = (iz+1) % resz;
	
	// x increment
	vector iq = set(ixp,iy,iz);
	float qr = volumeindex(0,re,iq);
	float qi = volumeindex(0,im,iq);
	if (qr<1e-6 && qi<1e-6) qr=1e-6;
	float a = pr*qr + pi*qi;
	float b = pr*qi - pi*qr;
	wx = hbar * atan2(b,a);
	
	// y increment
	iq = set(ix,iyp,iz);
	qr = volumeindex(0,re,iq);
	qi = volumeindex(0,im,iq);
	if (qr<1e-6 && qi<1e-6) qr=1e-6;
	a = pr*qr + pi*qi;
	b = pr*qi - pi*qr;
	wy = hbar * atan2(b,a);
	
	// z increment
	iq = set(ix,iy,izp);
	qr = volumeindex(0,re,iq);
	qi = volumeindex(0,im,iq);
	if (qr<1e-6 && qi<1e-6) qr=1e-6;
	a = pr*qr + pi*qi;
	b = pr*qi - pi*qr;
	wz = hbar * atan2(b,a);
}

// Gauge Transformation for Psi by phase q
void GaugeTransform( export float a,b,c,d; const float q )
{
	float cq = cos(q);
	float sq = sin(q);
	ComplexProd( a, b, cq, sq, a, b );
	ComplexProd( c, d, cq, sq, c, d );
}

// Schroedinger Evolve in Fourier Domain
void SchroedingerFlow(export float a,b,c,d; const float hbar; const float dt; 
    const int ix, iy, iz, resx, resy, resz; 
    const vector dPdx, dPdy, dPdz )
{
	float fac = -2.0 * PI*PI * hbar * dt;
	float ox = (ix-resx/2)/(resx * dPdx.x);
	float oy = (iy-resy/2)/(resy * dPdy.y);
	float oz = (iz-resz/2)/(resz * dPdz.z);
	float lambda = fac * (ox*ox + oy*oy + oz*oz);
	float cl = cos(lambda);
	float sl = sin(lambda);
	ComplexProd( a, b, cl, sl, a, b );
	ComplexProd( c, d, cl, sl, c, d );
}

// Compute S_x
void S3(const float a,b,c,d; export float s3)
{
	s3 = - a*a - b*b + c*c + d*d;
}

// Add Circular Vortex Ring by Biot Savart
void AddCircle(export float re,im; 
	const vector pos, center, normal; const float r)
{
	vector p = pos - center;
	vector n = normalize(normal);
	float a = dot(p,p) - r*r;
	float b = 2 * dot(p,n);
	ComplexProd( re,im,a,b,re,im );
}
// Add Circular Vortex Ring Locally
void AddCircleLocal(export float re,im;
                    const vector pos, center, normal;
                    const float r;
                    const float d)
{
    vector p = pos - center;
    vector n = normalize(normal);
    float alpha = 0.0;
    float z = dot(p,n);
    if ( dot(p,p) - z*z < r*r){
        if ( z>0 && z<=d/2 )
        {
            alpha = -PI*(2*z/d - 1);
        }
        else if (z<=0 && z>= -d/2)
        {
            alpha = -PI*(2*z/d + 1);
        }
    }
    float a = cos(alpha);
    float b = sin(alpha);
    ComplexProd( re, im, a, b, re, im);
}