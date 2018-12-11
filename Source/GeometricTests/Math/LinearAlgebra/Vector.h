/*
 * Copyright (c) 2018 Robert Slattery
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

#include "CoreMinimal.h"
#include <cstddef>
#include <cmath>

/*
 * Vector.h
 *
 * A vector is a directed line segment with
 * a direction and a magnitude (length). A vector
 * may have any nonnegative length.
*/

// 2D Vectors (Using for Intersection of Two Lines in 2D)
template<typename T>
struct Vector2
{

public:

	//~ 2D Vector Variables. Commonly use 4-byte float or 8-byte double.
	T X;	// 4 bytes
	T Y;	// 4 bytes
	//~ Total Memory Load: 8 bytes
	
	// Constructor: No Arguments
	inline explicit Vector2<T>() : X(0), Y(0) {}
	
	// Constructor: Passing Two Coordinates
	inline Vector2<T>(
		  const T& _X
		, const T& _Y) :
		  X(_X)
		, Y(_Y)
	{}
	
	inline constexpr T GetX() const { return X; }
	inline constexpr T GetY() const { return Y; }
	
	inline void Set(const T& _X, const T& _Y)
	{
		X = _X;
		Y = _Y;
	}
	
	// The magnitude (length) of a vector is the square root of the sum
	// of the squares of the components of the vector
	inline T Magnitude()
	{
		return sqrt(pow(X, 2) + pow(Y, 2));
	}
	
	// Vector Directions
	static const Vector2<T> ZeroVector;
	static const Vector2<T> ForwardVector;
	static const Vector2<T> BackwardVector;
	static const Vector2<T> RightVector;
	static const Vector2<T> LeftVector;
	
	//~ Begin Operator Overloads
	bool operator==(const Vector2<T>& Other) const
	{
		return (X == Other.X 
			 && Y == Other.Y);
	}

	bool operator!=(const Vector2<T>& Other) const
	{
		return (X != Other.X
			 || Y != Other.Y);
	}

	bool operator<(const Vector2<T>& Other) const
	{
		if (X != Other.X) return X < Other.X;
		else if (Y != Other.Y) return Y < Other.Y;
		else return false;
	}

	bool operator>(const Vector2<T>& Other) const
	{
		if (X != Other.X) return X > Other.X;
		else if (Y != Other.Y) return Y > Other.Y;
		else return false;
	}

	inline T& operator[](int Index)
	{
		check(Index >= 0 && Index < 2);	// In C++, use "assert".
		if (Index == 0) return X;
		else return Y;
	}

	inline T operator[](int Index) const
	{
		check(Index >= 0 && Index < 2); // In C++, use "assert".
		if (Index == 0) return X;
		else return Y;
	}

	inline Vector2<T> operator-() const
	{
		return{ -X, -Y };
	}

	inline Vector2<T> operator*(T Scalar) const
	{
		return Vector2<T>(X * Scalar, Y * Scalar);
	}

	inline Vector2<T> operator*(const Vector2<T>& Other) const
	{
		return Vector2<T>(X * Other.X, Y * Other.Y);
	}

	inline Vector2<T> operator*=(T Scalar) const
	{
		X *= Scalar;
		Y *= Scalar;
		return *this;
	}

	inline Vector2<T> operator*=(const Vector2<T>& Other) const
	{
		X *= Other.X;
		Y *= Other.Y;
		return *this;
	}
	//~ End Operator Overloads
		
	// Copy-Assignment Operator returns *this
	Vector2<T>& operator=(const Vector2<T>& Other) 
	{ 
		X = Other.X;
		Y = Other.Y;
		return *this; 
	}
};

// 3D Vectors
template<typename T>
struct Vector3
{

public:

	//~ 3D Vector Variables. Commonly use 4-byte float or 8-byte double.
	T X;	// 4 bytes
	T Y;	// 4 bytes
	T Z;	// 4 bytes
	//~ Total Memory Load: 12 bytes

	// Constructor: No Arguments
	inline explicit Vector3<T>() : X(0), Y(0), Z(0) {}
	
	// Constructor: Passing Three Coordinates
	inline Vector3<T>(
		  const T& _X
		, const T& _Y
		, const T& _Z) :
		  X(_X)
		, Y(_Y)
		, Z(_Z)
	{}
	
	inline constexpr T GetX() const { return X; }
	inline constexpr T GetY() const { return Y; }
	inline constexpr T GetZ() const { return Z; }

	inline void Set(const T& _X, const T& _Y, const T& _Z)
	{
		X = _X;
		Y = _Y;
		Z = _Z;
	}
	
	// The magnitude (length) of a vector is the square root of the sum
	// of the squares of the components of the vector
	inline T Magnitude()
	{
		return sqrt(pow(X, 2) + pow(Y, 2) + pow(Z, 2));
	}
	
	// Vector Directions
	static const Vector3<T> ZeroVector;
	static const Vector3<T> UpVector;
	static const Vector3<T> DownVector;
	static const Vector3<T> ForwardVector;
	static const Vector3<T> BackwardVector;
	static const Vector3<T> RightVector;
	static const Vector3<T> LeftVector;
	
	//~ Begin Operator Overloads
	inline bool operator==(const Vector3<T>& Other) const
	{
		return (X == Other.X 
			 && Y == Other.Y 
			 && Z == Other.Z);
	}

	inline bool operator!=(const Vector3<T>& Other) const
	{
		return (X != Other.X
			 || Y != Other.Y 
			 || Z != Other.Z);
	}

	inline bool operator<(const Vector3<T>& Other) const
	{
		if (X != Other.X) return X < Other.X;
		else if (Y != Other.Y) return Y < Other.Y;
		else if (Z != Other.Z) return Z < Other.Z;
		else return false;
	}

	inline bool operator>(const Vector3<T>& Other) const
	{
		if (X != Other.X) return X > Other.X;
		else if (Y != Other.Y) return Y > Other.Y;
		else if (Z != Other.Z) return Z > Other.Z;
		else return false;
	}

	inline T& operator[](int Index)
	{
		check(Index >= 0 && Index < 3);	// In C++, use "assert".
		if (Index == 0) return X;
		else if (Index == 1) return Y;
		else return Z;
	}

	inline T operator[](int Index) const
	{
		check(Index >= 0 && Index < 3); // In C++, use "assert".
		if (Index == 0) return X;
		else if (Index == 1) return Y;
		else return Z;
	}

	inline Vector3<T> operator+(const Vector3<T>& Other) const
	{
		return Vector3<T>(X + Other.X, Y + Other.Y, Z + Other.Z);
	}

	inline Vector3<T> operator+(float Bias) const
	{
		return Vector3<T>(X + Bias, Y + Bias, Z + Bias);
	}

	inline Vector3<T> operator-() const
	{
		return { -X, -Y, -Z };
	}

	inline Vector3<T> operator-(const Vector3<T>& Other) const
	{
		return Vector3<T>(X - Other.X, Y - Other.Y, Z - Other.Z);
	}

	inline Vector3<T> operator-(float Bias) const
	{
		return Vector3<T>(X - Bias, Y - Bias, Z - Bias);
	}

	inline Vector3<T> operator*(float Scalar) const
	{
		return Vector3<T>(X * Scalar, Y * Scalar, Z * Scalar);
	}

	inline Vector3<T> operator*(const Vector3<T>& Other) const
	{
		return Vector3<T>(X * Other.X, Y * Other.Y, Z * Other.Z);
	}

	inline Vector3<T> operator*=(float Scalar)
	{
		X *= Scalar;
		Y *= Scalar;
		Z *= Scalar;
		return *this;
	}

	inline Vector3<T> operator*=(const Vector3<T>& Other)
	{
		X *= Other.X;
		Y *= Other.Y;
		Z *= Other.Z;
		return *this;
	}

	inline Vector3<T> operator/(float Scalar) const
	{
		const float RScale = 1.f / Scalar;
		return Vector3<T>(X * RScale, Y * RScale, Z * RScale);
	}

	inline Vector3<T> operator/(const Vector3<T>& Other) const
	{
		return Vector3<T>(X / Other.X, Y / Other.Y, Z / Other.Z);
	}
	//~ End Operator Overloads
		
	// Copy-Assignment Operator returns *this
	Vector3<T>& operator=(const Vector3<T>& Other) 
    { 
        X = Other.X;
        Y = Other.Y;
        Z = Other.Z;
        return *this; 
    }
};

template<typename T> Vector2<T> const Vector2<T>::ZeroVector = Vector2<T>(0.0, 0.0);
template<typename T> Vector2<T> const Vector2<T>::ForwardVector = Vector2<T>(0.0, 1.0);
template<typename T> Vector2<T> const Vector2<T>::BackwardVector = Vector2<T>(0.0, -1.0);
template<typename T> Vector2<T> const Vector2<T>::RightVector = Vector2<T>(1.0, 0.0);
template<typename T> Vector2<T> const Vector2<T>::LeftVector = Vector2<T>(-1.0, 0.0);

template<typename T> Vector3<T> const Vector3<T>::ZeroVector = Vector3<T>(0.0, 0.0, 0.0);
template<typename T> Vector3<T> const Vector3<T>::UpVector = Vector3<T>(0.0, 0.0, 1.0);
template<typename T> Vector3<T> const Vector3<T>::DownVector = Vector3<T>(0.0, 0.0, -1.0);
template<typename T> Vector3<T> const Vector3<T>::ForwardVector = Vector3<T>(0.0, 1.0, 0.0);
template<typename T> Vector3<T> const Vector3<T>::BackwardVector = Vector3<T>(0.0, -1.0, 0.0);
template<typename T> Vector3<T> const Vector3<T>::RightVector = Vector3<T>(1.0, 0.0, 0.0);
template<typename T> Vector3<T> const Vector3<T>::LeftVector = Vector3<T>(-1.0, 0.0, 0.0);

template<typename T>
inline Vector2<T> operator*(const Vector2<T>& LHS, T Scalar)
{
	return LHS.operator*(Scalar);
}

template<typename T>
inline Vector3<T> operator*(const Vector3<T>& LHS, float Scalar)
{
	return LHS.operator*(Scalar);
}