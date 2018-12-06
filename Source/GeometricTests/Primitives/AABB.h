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
 * AABB.h
 *
 * Axis-Aligned Bounding Box
*/
template<typename T>
struct AABB
{
	
public:

	// Constructor
	explicit AABB<T>(
		  const Vector3<T>& _Min
		, const Vector3<T>& _Max) : 
		  Min(_Min)
		, Max(_Max)
	{}

	~AABB() = default;
	
	// Don't allow copying or moving
	AABB(AABB const&) = delete;
	AABB& operator=(const AABB&) = delete;
	AABB(AABB&&) = delete;
	AABB& operator=(AABB&&) = delete;
	
	inline void Clear() const noexcept
	{
		GetMin().X = GetMin().Y = GetMin().Z = FLT_MAX;
		GetMax().X = GetMax().Y = GetMax().Z = -FLT_MAX;
	}
	
	inline void Add(const Vector3<T>& Point)
	{
		if (Point.X < GetMin().X) GetMin().X = (T)Point.X;
		if (Point.X > GetMax().X) GetMax().X = (T)Point.X;
		
		if (Point.Y < GetMin().X) GetMin().Y = (T)Point.Y;
		if (Point.Y > GetMax().X) GetMax().Y = (T)Point.Y;
		
		if (Point.Z < GetMin().X) GetMin().Z = (T)Point.Z;
		if (Point.Z > GetMax().X) GetMax().Z = (T)Point.Z;
	}
	
	inline constexpr Vector3<T> GetMin() const noexcept { return Min; }
	inline void SetMin(const Vector3<T>& _Min) { return Min = _Min; }
	
	inline constexpr Vector3<T> GetMax() const noexcept { return Max; }
	inline void SetMax(const Vector3<T>& _Max) { return Max = _Max; }
		
private:

	Vector3<T> Min;
	Vector3<T> Max;

};
