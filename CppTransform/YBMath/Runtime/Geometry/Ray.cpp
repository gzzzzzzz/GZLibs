#include "Ray.h"

float Ray::SqrDistToPoint(const Vector3f& P) const {
//Vector v = L.P1 - L.P0; 
	Vector3f v = direction;
//	Vector3f w = P - L.P0; 
	Vector3f w = P - origin; 

	float c1 = Dot(w,v); 
	float c2 = Dot(v,v); 
	float b = c1 / c2; 

//Point Pb = L.P0 + b * v; 
	Vector3f Pb = GetPoint (b);
	return SqrMagnitude(P - Pb); 
}
