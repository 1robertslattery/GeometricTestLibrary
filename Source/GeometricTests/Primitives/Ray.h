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

#include "../Math/LinearAlgebra/Vector.h"

/*
 * Ray.h
 *
 * A ray is a directed line segment defined by:
 *
 * p(t) = p0 + t*d
 *
 * Where:
 *
 * p0 = ray origin vector
 * d = delta vector, which contains length and direction
 *
 * t is restricted to normalized range [0, 1]
*/
template<typename T>
struct Ray
{
	
public:

	// Constructor: No Arguments
	inline explicit Ray<T>() :
		  Origin(Vector3<T>::ZeroVector)
		, Delta(Vector3<T>::ZeroVector)
	{}

	// Constructor: Passing Origin and Delta
	inline explicit Ray<T>(
		  const Vector3<T>& _Origin
		, const Vector3<T>& _Delta) : 
		  Origin(_Origin)
		, Delta(_Delta)
	{}
		
	// Destructor
	~Ray<T>() = default;

	// Copy Constructor
	Ray<T>(Ray<T> const&) = delete;

	// Copy-Assignment Operator
	Ray<T>& operator=(const Ray<T>&) = delete;

	// Move Constructor
	Ray<T>(Ray<T>&&) = delete;

	// Move-Assignment Operator
	Ray<T>& operator=(Ray<T>&&) = delete;
	
	inline constexpr Vector3<T> GetOrigin() const { return Origin; }
	inline constexpr Vector3<T> GetDelta() const { return Delta; }

	inline void Set(const Vector3<T>& _Origin, const Vector3<T>& _Delta) 
	{
		Origin = _Origin;
		Delta = _Delta;
	}
	
private:
	
	//~ Ray Variables
	Vector3<T> Origin;	// p0
	Vector3<T> Delta;	// d
	//~ End Ray Variables
};
