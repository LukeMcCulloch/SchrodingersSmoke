// Complex Product
void ComplexProd(const float a,b,c,d; export float outRe, outIm )
{
	float tmpRe = a*c - b*d;
	float tmpIm = b*c + a*d;
	outRe = tmpRe;
	outIm = tmpIm;
}

// Normalize
void Normalize(export float a,b,c,d )
{
	float fac = 1.0/sqrt( a*a + b*b + c*c + d*d );
	a *= fac;
	b *= fac;
	c *= fac;
	d *= fac;
}
