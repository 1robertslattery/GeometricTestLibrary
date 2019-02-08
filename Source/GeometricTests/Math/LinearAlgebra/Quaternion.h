/*
 * Copyright (c) 2018-2019 Robert Slattery
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#pragma once

#include <cmath>
#include <utility>
#include "Vector.h"

/*
 * Quaternion.h
 *
 * A quaternion is a 4D Vector that contains
 * a scalar component (W) and a 3D Vector component (V): 
 * 
 * { W, V } or { W, (X, Y, Z) }
 *
 * A quaternion contains an axis and an angle. Rotation 
 * values are stored in degrees.
*/

template<typename T>
struct Quaternion
{

public:

	//~ Quaternion Variables. Commonly use 4-byte float or 8-byte double.
	T W;	// 4 bytes
	T X;	// 4 bytes
	T Y;	// 4 bytes
	T Z;	// 4 bytes
	//~ Total Memory Load: 16 bytes

	// Constructor: No Arguments, Zero Rotation
	inline explicit Quaternion<T>() : W(0), X(0), Y(0), Z(0) {}
	
	// Constructor: Passing Four Coordinates
	inline constexpr Quaternion<T>(
		  const T& _W
		, const T& _X
		, const T& _Y
		, const T& _Z) :
		  W(_W)
		, X(_X)
		, Y(_Y)
		, Z(_Z)
	{}
	
	// Constructor: Convert from axis-angle format in radians
	inline explicit Quaternion<T>(const Vector3<T>& Axis, const float Radians)
	{
		const float HalfAngle = 0.5f * Radians;
		W = std::cos(HalfAngle);
		X = Axis.GetX() * std::sin(HalfAngle);
		Y = Axis.GetY() * std::sin(HalfAngle);
		Z = Axis.GetZ() * std::sin(HalfAngle);
	}
	
	~Quaternion<T>() = default;
	
	inline constexpr T GetW() const { return W; }
	inline constexpr T GetX() const { return X; }
	inline constexpr T GetY() const { return Y; }
	inline constexpr T GetZ() const { return Z; }
	
	inline T Magnitude()
	{
		return sqrt(pow(W, 2) + pow(X, 2) + pow(Y, 2) + pow(Z, 2));
	}

	inline T DotProduct(const Quaternion<T>& LHS, const Quaternion<T>& RHS)
	{
		return (LHS.W*RHS.W + LHS.X*RHS.X + LHS.Y*RHS.Y + LHS.Z*RHS.Z);
	}

	inline Quaternion<T> Inverse(const Quaternion<T>& Quat)
	{
		return Quaternion<T>(Quat.W, -Quat.X, -Quat.Y, -Quat.Z);
	}
	
	static const Quaternion<T> ZeroRotation;
	
	//~ Begin Operator Overloads
	inline bool operator==(const Quaternion<T>& Other) const
	{
		return (W == Other.W
			 && X == Other.X 
			 && Y == Other.Y 
			 && Z == Other.Z);
	}

	inline bool operator!=(const Quaternion<T>& Other) const
	{
		return (W != Other.W
			 || X != Other.X
			 || Y != Other.Y 
			 || Z != Other.Z);
	}

	inline Quaternion<T> operator+(const Quaternion<T>& Other) const
	{
		return Quaternion<T>(W + Other.W, X + Other.X, Y + Other.Y, Z + Other.Z);
	}

	inline Quaternion<T> operator+(float Bias) const
	{
		return Quaternion<T>(W + Bias, X + Bias, Y + Bias, Z + Bias);
	}

	inline Quaternion<T> operator-() const
	{
		return { -W, -X, -Y, -Z };
	}

	inline Quaternion<T> operator-(const Quaternion<T>& Other) const
	{
		return Quaternion<T>(W - Other.W, X - Other.X, Y - Other.Y, Z - Other.Z);
	}

	inline Quaternion<T> operator-(float Bias) const
	{
		return Quaternion<T>(W - Bias, X - Bias, Y - Bias, Z - Bias);
	}

	inline Quaternion<T> operator*(float Scalar) const
	{
		return Quaternion<T>(W * Scalar, X * Scalar, Y * Scalar, Z * Scalar);
	}

	inline Quaternion<T> operator*(const Quaternion<T>& Other) const
	{
		return Quaternion<T>(W * Other.W, X * Other.X, Y * Other.Y, Z * Other.Z);
	}

	inline Quaternion<T> operator*=(float Scalar)
	{
		W *= Scalar;
		X *= Scalar;
		Y *= Scalar;
		Z *= Scalar;
		return *this;
	}

	inline Quaternion<T> operator*=(const Quaternion<T>& Other)
	{
		W *= Other.W;
		X *= Other.X;
		Y *= Other.Y;
		Z *= Other.Z;
		return *this;
	}
	//~ End Operator Overloads
		
	// Copy Constructor
	Quaternion<T>(const Quaternion<T>& Other) : 
		  W(Other.W)
		, X(Other.X)
		, Y(Other.Y)
		, Z(Other.Z)
	{}
	
	// Copy-Assignment Operator
	Quaternion<T>& operator=(const Quaternion<T>& Other) 
    { 
		if (this != &Other)
		{
			X = Other.X;
			Y = Other.Y;
			Z = Other.Z;
		}
		
        return *this; 
    }
	
	// Move Constructor
	Quaternion<T>(Quaternion<T>&& Other) : 
		  W(0)
		, X(0)
		, Y(0)
		, Z(0)
	{	*this = std::move(Other);	}
	
	// Move-Assignment Operator
	Quaternion<T>& operator=(Quaternion<T>&& Other) 
    { 
		if (this != &Other)
		{
			W = Other.W;
			X = Other.X;
			Y = Other.Y;
			Z = Other.Z;
			
			Other.W = 0;
			Other.X = 0;
			Other.Y = 0;
			Other.Z = 0;
		}
		
		return *this;
    }
};

template<typename T> Quaternion<T> const Vector3<T>::ZeroRotation = Quaternion<T>(0, 0, 0, 0);

template<typename T>
inline Quaternion<T> operator*(const Quaternion<T>& LHS, float Scalar)
{
	return LHS.operator*(Scalar);
}