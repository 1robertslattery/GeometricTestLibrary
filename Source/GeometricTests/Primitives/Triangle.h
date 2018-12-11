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
 * Triangle.h
 *
 * Defined by vertices A, B, and C:
 *
 * Triangle(T) = [A, B, C]
 *
 * Area of Triangle = b*h / 2
*/
template<typename T>
struct Triangle
{
	
public:

	// Constructor: Passing A, B, C
	inline explicit Triangle<T>(
		  const Vector3<T>& _VertexA
		, const Vector3<T>& _VertexB
		, const Vector3<T>& _VertexC) : 
		  VertexPosition[0](_VertexA)
		, VertexPosition[1](_VertexB)
		, VertexPosition[2](_VertexC)
	{}
	
	// Destructor
	~Triangle<T>() = default;
	
	// Copy Constructor
	Triangle<T>(Triangle<T> const&) = delete;

	// Copy-Assignment Operator
	Triangle<T>& operator=(const Triangle<T>&) = delete;

	// Move Constructor
	Triangle<T>(Triangle<T>&&) = delete;

	// Move-Assignment Operator
	Triangle<T>& operator=(Triangle<T>&&) = delete;
	
	inline constexpr Vector3<T> GetVertexA() const { return VertexPosition[0]; }
	inline constexpr Vector3<T> GetVertexB() const { return VertexPosition[1]; }
	inline constexpr Vector3<T> GetVertexC() const { return VertexPosition[2]; }
			
private:

	//~ Triangle Variables
	Vector3<T> VertexPosition[3];
	//~ End Triangle Variables
};
