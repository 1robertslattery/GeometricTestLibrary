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
 * Circle.h
 *
 * Defined by a center and radius:
 *
 * (x - h)^2 + (y - k)^2 = r^2
 *
 * Where:
 *
 * (h, k) = center
 * r = radius
*/
template<typename T>
struct Circle
{
	
public:

	// Constructor: No Arguments
	inline explicit Circle<T>() :
		  Radius(0.f)
		, Center(Vector3<T>::ZeroVector)
	{}

	// Constructor: Passing Radius and Center
	inline explicit Circle<T>(
		  float _Radius
		, Vector3<T> _Center = Vector3<T>::ZeroVector) : 
		  Radius(_Radius)
		, Center(_Center)
	{}
	
	// Destructor
	~Circle<T>() = default;

	// Copy Constructor
	Circle<T>(Circle<T> const&) = delete;

	// Copy-Assignment Operator
	Circle<T>& operator=(const Circle<T>&) = delete;

	// Move Constructor
	Circle<T>(Circle<T>&&) = delete;

	// Move-Assignment Operator
	Circle<T>& operator=(Circle<T>&&) = delete;
	
	inline constexpr float GetRadius() const { return Radius; }
	inline constexpr Vector3<T> GetCenter() const { return Center; }

	inline void Set(const float _Radius, const Vector3<T>& _Center) 
	{
		Radius = _Radius;
		Center = _Center;
	}
			
private:

	//~ Circle Variables
	float Radius;			// r (4 bytes)
	Vector3<T> Center;		// (h, k)
	//~ End Circle Variables
	
};
