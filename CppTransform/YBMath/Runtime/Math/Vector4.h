#pragma once

#include "Vector3.h"

class Vector4f
{
public:
	Vector4f() {}
	Vector4f( float inX, float inY, float inZ, float inW ) : x(inX), y(inY), z(inZ), w(inW) {}
	explicit Vector4f( const Vector3f& v, float inW ) : x(v.x), y(v.y), z(v.z), w(inW) {}
	explicit Vector4f( const float* v ) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}

	void Set( float inX, float inY, float inZ, float inW ) { x = inX; y = inY; z = inZ; w = inW; }
	void Set( const float* array ) { x = array[0]; y = array[1]; z = array[2]; w = array[3]; }

	float* GetPtr() 			{ return &x; }
	const float* GetPtr() const { return &x; }

	float& operator[] (int i)						{ DebugAssertIf (i < 0 || i > 3); return (&x)[i]; }
	const float& operator[] (int i)const			{ DebugAssertIf (i < 0 || i > 3); return (&x)[i]; }

	bool operator == (const Vector4f& v) const		{ return x == v.x && y == v.y && z == v.z && w == v.w; }
	bool operator != (const Vector4f& v) const		{ return x != v.x || y != v.y || z != v.z || w != v.w; }
	bool operator == (const float v[4]) const		{ return x == v[0] && y == v[1] && z == v[2] && w == v[3]; }
	bool operator != (const float v[4]) const		{ return x != v[0] || y != v[1] || z != v[2] || w != v[3]; }

	Vector4f operator - () const					{ return Vector4f (-x, -y, -z, -w); }
	/*
	DEFINE_GET_TYPESTRING_IS_ANIMATION_CHANNEL (Vector4f)
	template<class TransferFunction> void Transfer (TransferFunction& transfer);
	*/
	float x;
	float y;
	float z;
	float w;
};


inline Vector4f operator * (const Vector4f& lhs, const Vector4f& rhs)	{ return Vector4f (lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z, lhs.w * rhs.w); }
inline Vector4f operator * (const Vector4f& inV, const float s)			{ return Vector4f (inV.x * s, inV.y * s, inV.z * s, inV.w * s); }
inline Vector4f operator + (const Vector4f& lhs, const Vector4f& rhs)	{ return Vector4f (lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w); }
inline Vector4f operator - (const Vector4f& lhs, const Vector4f& rhs)	{ return Vector4f (lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w); }
inline float Dot (const Vector4f& lhs, const Vector4f& rhs)				{ return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w; }

inline Vector4f Lerp (const Vector4f& from, const Vector4f& to, float t) { return to * t + from * (1.0F - t); }
/*
template<class TransferFunction>
inline void Vector4f::Transfer (TransferFunction& t)
{
	t.AddMetaFlag (kTransferUsingFlowMappingStyle);
	t.Transfer (x, "x");
	t.Transfer (y, "y");
	t.Transfer (z, "z");
	t.Transfer (w, "w");
}
*/

inline
std::ostream& operator << (std::ostream& out, const Vector4f& p)
{
	out << "[" << p.x << ", " << p.y << ", " << p.z << ", " << p.w << "]";
	return out;
}