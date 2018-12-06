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

	explicit Circle<T>(
		  float _Radius
		, Vector3<T> _Center = Vector3<T>::ZeroVector) : 
		  Radius(_Radius)
		, Center(_Center)
	{}
	
	~Circle() = default;

	// Don't allow copying or moving
	Circle(Circle const&) = delete;
	Circle& operator=(const Circle&) = delete;
	Circle(Circle&&) = delete;
	Circle& operator=(Circle&&) = delete;
	
	inline constexpr float GetRadius() const noexcept { return Radius; }
	inline void SetRadius(const float _Radius) { return Radius = _Radius; }
	
	inline constexpr Vector3<T> GetCenter() const noexcept { return Center; }
	inline void SetCenter(const Vector3<T>& _Center) { return Center = _Center; }
			
private:

	//~ Circle Data
	float Radius;			// r
	Vector3<T> Center;		// (h, k)
	//~ Total Memory Load = 4 bytes
	
};
