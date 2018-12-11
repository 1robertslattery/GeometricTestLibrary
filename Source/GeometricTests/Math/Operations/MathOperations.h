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

#include "../LinearAlgebra/Vector.h"
#include <cstddef>
#include <cmath>

/*
 * MathOperations.h
 *
 * Non-member functions for vector math
*/
namespace MyMathLibrary
{
	// Solves the dot product of two 3D vectors
	template<typename T>
	inline T DotProduct(const Vector3<T>& LHS, const Vector3<T>& RHS, const int Size = 3)
	{
		T Result = 0.0;
		
		for (int i = 0; i < Size; ++i)
			Result += LHS[i] * RHS[i];
		
		return Result;
	}

	// Solves the cross product of two 3D vectors
	template<typename T>
	inline Vector3<T> CrossProduct(const Vector3<T>& LHS, const Vector3<T>& RHS)
	{
		return Vector3<T>(LHS.Y*RHS.Z - LHS.Z*RHS.Y,
						  LHS.Z*RHS.X - LHS.X*RHS.Z,
						  LHS.X*RHS.Y - LHS.Y*RHS.X);
	}

	// Solves the normal of a 3D vector, a.k.a. unit vector
	template<typename T>
	inline Vector3<T> Normalize(const Vector3<T>& Vector)
	{
		Vector3<T> Result = Vector3<T>::ZeroVector;
		
		// Find the Vector Magnitude
		T Magnitude = sqrt(pow(Vector.X, 2) + pow(Vector.Y, 2) + pow(Vector.Z, 2));
		
		Result = Vector / Magnitude;
		
		return Result;
	}

	// 3D Vector Distance Formula
	template <typename T>
	inline T Distance(const Vector3<T>& LHS, const Vector3<T>& RHS)
	{
		T DiffX = LHS.X - RHS.X;
		T DiffY = LHS.Y - RHS.Y;
		T DiffZ = LHS.Z - RHS.Z;
		return sqrt((pow(DiffX, 2) + pow(DiffY, 2) + pow(DiffZ, 2)));
	}
	
	// Solves the best-fit plane normal for a set of points
	template <typename T>
	inline Vector3<T> BestFitNormal(const Vector3<T> Vector[], int N)
	{
		Vector3<T> Result = Vector3<T>::ZeroVector;
		
		Vector3<T> *Previous = &Vector[N - 1];
		
		for (int i = 0; i < N; ++i)
		{
			const Vector3<T> *Current = &Vector[i];
			
			Result.X += (Previous->Z + Current->Z) * (Previous->Y - Current->Y);
			Result.Y += (Previous->X + Current->X) * (Previous->Z - Current->Z);
			Result.Z += (Previous->Y + Current->Y) * (Previous->X - Current->X);
			
			Previous = Current;
		}
		
		Normalize(Result);
		return Result;
	}
	
	// Center of gravity of a triangle
	template <typename T>
	inline Vector3<T> CenterOfGravity(const Vector3<T>& V1, const Vector3<T>& V2, const Vector3<T>& V3)
	{
		return ((V1 + V2 + V3) / 3);
	}
		
	// Incenter of a triangle
	template <typename T>
	inline Vector3<T> InCenter(
		  const Vector3<T>& V1
		, const Vector3<T>& V2
		, const Vector3<T>& V3
		, const Vector3<T>& L1
		, const Vector3<T>& L2
		, const Vector3<T>& L3)
	{
		Vector3<T> Perimeter = L1 + L2 + L3;
		Vector3<T> Numerator = (L1*V1) + (L2*V2) + (L3*V3);
		return (Numerator / Perimeter);
	}
	
	/*
	 * TODO: Circumcenter of a triangle
	 *
	 * Solve for t:
	 *
	 * t = (c2 + c3)*v1 + (c3 + c1)*v2 + (c1 + c2)*v3 / 2*c
	 *
	 * Where:
	 *
	 * v1, v2, v3 = vertices of triangle
	 * e1, e2, e3 = edge of triangle
	 * d1 = -e2*e3
	 * d2 = -e3*e1
	 * d3 = -e1*e2
	 * c1 = d2*d3
	 * c2 = d3*d1
	 * c3 = d1*d2
	 * c = c1 + c2 + c3
	*/
} 

