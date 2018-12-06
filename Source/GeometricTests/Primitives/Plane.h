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
#include "../Math/LinearAlgebra/Vector.h"
#include "../Math/Operations/MathOperations.h"

/*
 * Plane.h
 *
 * An infinite surface defined in 3D space by:
 *
 * ax + by + cz + d = 0
 */
template<typename T>
struct Plane
{

public:
	
	Vector3<T> Normal;
	float D;		// 4 bytes
	
	// Constructor: Using 3 vector points
	explicit Plane<T>(const T& _A, const T& _B, const T& _C, const T& _D)
	{
		Normal.X = (T)_A;
        Normal.Y = (T)_B;
        Normal.Z = (T)_C;
		D = (T)_D;
	}
	
	// Constructor: Using a vector point with a vector normal
	explicit Plane<T>(const Vector3<T>& _Normal, const Vector3<T>& _Point)
	{
		Normal = _Normal;
		D = -MyMathLibrary::DotProduct(_Normal, _Point);
	}
	
	~Plane() = default;
	
	// Don't allow copying or moving
	Plane(Plane const&) = delete;
	Plane& operator=(const Plane&) = delete;
	Plane(Plane&&) = delete;
	Plane& operator=(Plane&&) = delete;
};