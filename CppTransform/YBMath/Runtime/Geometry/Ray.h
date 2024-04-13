#ifndef RAY_H
#define RAY_H

#include "Runtime/Math/FloatConversion.h"
#include "Runtime/Math/Vector3.h"

class Ray
{
public:
	Vector3f	origin;
	Vector3f	direction; // Direction is always normalized
	
public:
	Ray () {}
	Ray (const Vector3f& orig, const Vector3f& dir) { origin = orig; SetDirection (dir); }
		
	const Vector3f& GetDirection ()const		{ return direction; }
	// Direction has to be normalized
	void SetDirection (const Vector3f& dir)	{ AssertIf (!IsNormalized (dir)); direction = dir; }
	void SetApproxDirection (const Vector3f& dir)	{ direction = NormalizeFast (dir); }
	void SetOrigin (const Vector3f& origin)	{ this->origin = origin; }

	const Vector3f& GetOrigin ()const		{ return origin; }
	Vector3f GetPoint (float t) const			{ return origin + t * direction; }

	float SqrDistToPoint (const Vector3f &v) const;
};


#endif
