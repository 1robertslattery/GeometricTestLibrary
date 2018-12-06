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
 * p(t) = p0 + t*d
 *
 * Where:
 *
 * p0 = ray origin vector
 * d = ray direction vector
*/
template<typename T>
struct Ray
{
	
public:

	explicit Ray<T>(
		  const Vector3<T>& _Origin
		, const Vector3<T>& _Direction) : 
		  Origin(_Origin)
		, Direction(_Direction)
	{}
		
	~Ray() = default;

	// No copy or move
	Ray(Ray const&) = delete;
	Ray& operator=(const Ray&) = delete;
	Ray(Ray&&) = delete;
	Ray& operator=(Ray&&) = delete;
	
	inline constexpr Vector3<T> GetOrigin() const noexcept { return Origin; }
	inline void SetOrigin(const Vector3<T>& _Origin) { return Origin = _Origin; }
	
	inline constexpr Vector3<T> GetDirection() const noexcept { return Direction; }
	inline void SetDirection(const Vector3<T>& _Direction) { return Direction = _Direction; }
	
private:
	
	Vector3<T> Origin;		// p0
	Vector3<T> Direction;	// d
	
};
