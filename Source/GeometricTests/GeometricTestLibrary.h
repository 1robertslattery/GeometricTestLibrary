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

#include "CoreMinimal.h"
#include "UObject/ObjectMacros.h"
#include "UObject/Object.h"
#include "UObjectGlobals.h"
#include "Math/LinearAlgebra/Vector.h"
#include "Math/Operations/MathOperations.h"
#include "Primitives/AABB.h"
#include "Primitives/Circle.h"
#include "Primitives/Plane.h"
#include "Primitives/Ray.h"
#include "Engine.h"
#include <algorithm>
#include <cmath>
#include "GeometricTestLibrary.generated.h"

/*
 * GeometricTestLibrary.h
 *
 * Solves a series of 2D and 3D math problems,
 * including, closest point, intersection
 * and line of sight.
 *
 * Examples:
 *
 * 1. Closest Point in AABB
 * 2. Closest Point on a Ray
 * 3. Closest Point on a Plane
 * 4. Closest Point on a Sphere
 * 5. Intersection AABB-Plane
 * 6. Intersection of Two Lines in 2D
 * 7. Intersection Ray-AABB (Common for Line of Sight Tests)
 * 8. Intersection Ray-Plane
 * 9. Intersection Ray-Triangle
 * 10. Intersection of Two Rays in 3D
 * 11. Dynamic Intersection Sphere-Plane
 * 12. Static Intersection Sphere-Plane
 * 13. Dynamic Intersection of Two Spheres
 * 14. Static Intersection of Two Spheres
 * 15. Intersection of Three Planes
 * 16. Intersection of Two AABBs
 * 17. Reflection Vector
 * 18. Barycentric Coordinates of Triangle in 3D
 *
 * Updated January 10, 2019
*/
UCLASS()
class GEOMETRICTESTS_API UGeometricTestLibrary : public UObject
{
	GENERATED_BODY()

public:

	// Closest Point in AABB
	template<typename T>
	static Vector3<T> ClosestPointInAABB(
		  const Vector3<T>& Point
		, const AABB<T>& AABB_Ref);
		
	// Closest Point on a Ray
	template<typename T>
	static Vector3<T> ClosestPointOnRay(
		  const Vector3<T>& Point
		, const Vector3<T>& RayOrigin
		, const Vector3<T>& RayDelta);
	
	// Closest Point On Plane
	template<typename T>
	static Vector3<T> ClosestPointOnPlane(
		  const Vector3<T>& Point			// q
		, const Vector3<T>& PlaneNormal 	// n
		, float PlaneD);					// p*PlaneNormal = PlaneD	
	
	// Closest Point On Sphere
	template<typename T>
	static Vector3<T> ClosestPointOnSphere(
		  const Vector3<T>& Point
		, const Vector3<T>& SphereCenter
		, float SphereRadius);
		
	// Intersection AABB-Plane
	template<typename T>	
	static int DoesAABBPlaneIntersect(
		  const Vector3<T>& PlaneNormal
		, float D
		, const AABB<T>& AABB_Ref);
	
	// Intersection of Two Lines in 2D
	template<typename T>
	static Vector2<T> DoLinesIntersect(
		  const Vector2<T>& Origin1
		, const Vector2<T>& Delta1
		, const Vector2<T>& Origin2
		, const Vector2<T>& Delta2);
	
	// Intersection Ray-AABB
	template<typename T>
	static float DoesRayAABBIntersect(
		  const Vector3<T>& RayOrigin
		, const Vector3<T>& RayDelta
		, const AABB<T>& AABB_Ref
		, Vector3<T>* ReturnNormal);
	
	// Intersection Ray-Plane
	template<typename T>
	static bool DoesRayPlaneIntersect(
		  const Vector3<T>& RayOrigin
		, const Vector3<T>& RayDelta
		, const Vector3<T>& SurfaceNormal
		, float PlaneD);
		
	// Intersection Ray-Sphere
	template<typename T>	
	static bool DoesRaySphereIntersect(
		  const Vector3<T>& SphereCenter	
		, float SphereRadius					
		, const Vector3<T>& RayOrigin
		, const Vector3<T>& RayDelta);
		
	// Intersection Ray-Triangle
	template<typename T>	
	static float DoesRayTriangleIntersect(
		  const Vector3<T>& RayOrigin		// origin of ray
		, const Vector3<T>& RayDelta		// direction and length of ray
		, const Vector3<T>& Vertex1			// triangle vertices
		, const Vector3<T>& Vertex2			// .
		, const Vector3<T>& Vertex3			// .
		, float MinT);						// closest intersection found so far 
	
	// Intersection Two 3D Rays	
	template<typename T>
	static bool DoRaysIntersect(
		  const Vector3<T>& RayOrigin1		// p1
		, const Vector3<T>& RayDelta1		// d1
		, const Vector3<T>& RayOrigin2		// p2
		, const Vector3<T>& RayDelta2);		// d2
	
	// Dynamic Intersection Sphere-Plane
	template<typename T>
	static bool DoesSpherePlaneIntersect_Dynamic(
		  const Vector3<T>& PlaneNormal	// must be normalized first
		, float PlaneD							
		, const Vector3<T>& SphereDeltaVector
		, const Vector3<T>& SphereCenter
		, float SphereRadius);
		
	// Static Intersection Sphere-Plane
	template<typename T>
	static int DoesSpherePlaneIntersect_Static(
		  const Vector3<T>& PlaneNormal		// must be normalized first
		, float PlaneD						// p*PlaneNormal = PlaneD
		, const Vector3<T>& SphereCenter	// center of sphere
		, float SphereRadius);				// radius of sphere 
	
	// Dynamic Intersection of Two Spheres
	template<typename T>
	static bool DoSpheresIntersect_Dynamic(
		  const Vector3<T>& StationaryDeltaVector
		, const Vector3<T>& MovingDeltaVector
		, const Vector3<T>& StationarySphereCenter
		, const Vector3<T>& MovingSphereCenter
		, float StationarySphereRadius
		, float MovingSphereRadius);
	
	// Static Intersection of Two Spheres
	template<typename T>
	static bool DoSpheresIntersect_Static(
		  const Vector3<T>& SphereCenter1
		, float SphereRadius1
		, const Vector3<T>& SphereCenter2
		, float SphereRadius2);
		
	// Intersection of Three Planes
	template<typename T>
	static Vector3<T> PointOfIntersection(
		  const Vector3<T>& PlaneNormal1	// must be normalized first
		, float PlaneD1						// p*PlaneNormal = PlaneD
		, const Vector3<T>& PlaneNormal2	// must be normalized first
		, float PlaneD2						// p*PlaneNormal2 = PlaneD2	
		, const Vector3<T>& PlaneNormal3	// must be normalized first
		, float PlaneD3);					// p*PlaneNormal3 = PlaneD3
			
	// Intersection of Two AABBs
	template<typename T>
	static bool DoAABBsIntersect(
		  const AABB<T>& A
		, const AABB<T>& B);
		
	// Reflection Vector 
	template<typename T>
	static Vector3<T> SolveReflectionVector(
		  const Vector3<T>& SurfaceNormal
		, const Vector3<T>& LightDelta);
		
	// Barycentric Coordinates of Triangle in 3D
	template<typename T>
	static Vector3<T> SolveBarycentricCoordinates3D(
		  const Vector3<T> Vertices[3]		// vertices of the triangle
		, const Vector3<T>& Point			// p
		, float Barycentric[3]);			// barycentric coordinates
		
protected:

	// Constructor
	UGeometricTestLibrary() = default;
	
	// Destructor
	~UGeometricTestLibrary() = default;

	// UE4: UObject handles copying and moving
};

/*
 * Closest Point in an AABB (Axis-Aligned Bounding Box)
 *
 * To find the closest point on an AABB, we just clamp
 * the X, Y and Z positions of the point to the AABB
 *
 * q = given point
*/
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::ClosestPointInAABB(
	  const Vector3<T>& Point	// q
	, const AABB<T>& AABB_Ref) 
{
	Vector3<T> Result = Vector3<T>::ZeroVector;
				
	// Print q
	GEngine->AddOnScreenDebugMessage(1, 30.f, FColor::Red, TEXT("\nPoint q is..."));
	GEngine->AddOnScreenDebugMessage(
		  2, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Point.X, Point.Y, Point.Z));

	// Solve by "pushing" q onto B along each axis
	if (Point.X < AABB_Ref.GetMin().X)		Result.X = (T)AABB_Ref.GetMin().X;
	else if (Point.X > AABB_Ref.GetMax().X)	Result.X = (T)AABB_Ref.GetMax().X;
	
	if (Point.Y < AABB_Ref.GetMin().Y) 	   	Result.Y = (T)AABB_Ref.GetMin().Y;
	else if (Point.Y > AABB_Ref.GetMax().Y) Result.Y = (T)AABB_Ref.GetMax().Y;
	
	if (Point.Z < AABB_Ref.GetMin().Z) 		Result.Z = (T)AABB_Ref.GetMin().Z;
	else if (Point.Z > AABB_Ref.GetMax().Z) Result.Z = (T)AABB_Ref.GetMax().Z;
	
	// Print Result
	GEngine->AddOnScreenDebugMessage(3, 30.f, FColor::Red, TEXT("Closest point in AABB is..."));
	GEngine->AddOnScreenDebugMessage(
		  4, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
	
	// If Result is already inside the box, this returns the original point
	return Result;
}	

/* 
 * Closest Point on a Ray
 * 
 * Given a point q, we wish to find the point q`
 * which is the result of projecting q onto a ray R
 * 
 * t = Ď * v = Ď * (q - p0)
 * q` = p(t) = p0 + t*Ď = p0 + (Ď * (q - p0))Ď
 *
 * Where:
 *
 * q = given point
 * p0 = ray origin
 * Ď = normalized ray delta vector, a.k.a. unit vector
 * v = q - p0
*/
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::ClosestPointOnRay(
	  const Vector3<T>& Point			// q
	, const Vector3<T>& RayOrigin		// p0
	, const Vector3<T>& RayDelta)		// D (unnormalized)
{
	Vector3<T> Result = Vector3<T>::ZeroVector;

	// Print q
	GEngine->AddOnScreenDebugMessage(1, 30.f, FColor::Red, TEXT("\nPoint q is..."));
	GEngine->AddOnScreenDebugMessage(
		  2, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Point.X, Point.Y, Point.Z));
	
	// Solve Ď
	Vector3<T> D = MyMathLibrary::Normalize(RayDelta);
	
	// Solve v = q - p0
	Vector3<T> v = Point - RayOrigin;
	
	// Solve t
	float t = MyMathLibrary::DotProduct(D, v);

	// Solve R, the parametric ray
	Vector3<T> R = RayOrigin + (D*t);
	
	// If t < 0 or t > length of ray, then q` is not within the ray,
	// in which case, the closest point to q on R will 
	// be the ray origin (if t < 0) or endpoint (if t > length of ray)
	if (t < 0)
	{
		GEngine->AddOnScreenDebugMessage(
			  3, 30.f, FColor::Red
			, TEXT("\nResult not within the ray! Closest point is ray origin!, GeometricTestLibrary.h:323\n"));

		return Result;
	}
	else if (t > R.Magnitude()) // t > length of ray
	{
		GEngine->AddOnScreenDebugMessage(
			  4, 30.f, FColor::Red
			, TEXT("\nResult not within the ray! Closest point is the endpoint!, GeometricTestLibrary.h:331\n"));

		return Result;
	}
		
	// Solve q`
	Result = R;

	// Print q`
	GEngine->AddOnScreenDebugMessage(5, 30.f, FColor::Red, TEXT("\nClosest point to q is..."));
	GEngine->AddOnScreenDebugMessage(
		  6, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
	
	// Return closest point on parametic ray
	return Result;
}		

/* 
 * Closest Point on a Plane
 *
 * Given a point q, we wish to find the point q`
 * which is the result of projecting q onto a plane P
 *
 * q` = q + (d - q * ň)ň
 *
 * Where:
 *
 * q = given point
 * ň = normalized plane normal, a.k.a. unit vector
 * d = p*n = the plane equation, the distance from the plane to the point
 *
 * NOTE: We will flip d - q and then take inverse dot product
*/
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::ClosestPointOnPlane(
	  const Vector3<T>& Point			// q
	, const Vector3<T>& PlaneNormal		// n (unnormalized)
	, float PlaneD)						// p*PlaneNormal = PlaneD			
{
	Vector3<T> Result = Vector3<T>::ZeroVector;

	// Print q
	GEngine->AddOnScreenDebugMessage(1, 30.f, FColor::Red, TEXT("Point q is..."));
	GEngine->AddOnScreenDebugMessage(
		  2, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Point.X, Point.Y, Point.Z));
	
	// Solve ň
	Vector3<T> NormalizedPlaneNormal = MyMathLibrary::Normalize(PlaneNormal);
	
	// Solve q - d, flipping sign of d
	Vector3<T> Diff = Point - PlaneD;

	// Solve q - d * ň, with inverse dot product
	// 
	// Get the shortest distance between a point and a plane. The output is signed so it holds information
	// as to which side of the plane normal the point is.
	float Distance = -MyMathLibrary::DotProduct(Diff, NormalizedPlaneNormal);
		
	// Solve (d - q * ň)ň
	NormalizedPlaneNormal *= Distance;
		
	// Solve q`
	Result = Point + NormalizedPlaneNormal;
	
	// Print q`
	GEngine->AddOnScreenDebugMessage(3, 30.f, FColor::Red, TEXT("Closest Point to q is..."));
	GEngine->AddOnScreenDebugMessage(
		  4, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
	
	// Return closest point on plane
	return Result;
}

/*
 * Closest Point on a Sphere
 *
 * Given a point q, we wish to find q`, which
 * is the closest point on the circle to q
 *
 * q` = q + d * (||d|| - r / ||d||)
 *
 * Where:
 *
 * q = Given Point
 * r = radius of sphere
 * c = center of sphere
 * d = c - q
 * ||d|| = Magnitude of d
*/
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::ClosestPointOnSphere(
	  const Vector3<T>& Point			// q
	, const Vector3<T>& SphereCenter	// c
	, float SphereRadius)				// r
{
	Vector3<T> Result = Vector3<T>::ZeroVector;

	// Print q
	GEngine->AddOnScreenDebugMessage(1, 30.f, FColor::Red, TEXT("\nPoint q is..."));
	GEngine->AddOnScreenDebugMessage(
		  2, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Point.X, Point.Y, Point.Z));
		
	// Solve d = c - q
	Vector3<T> d = SphereCenter - Point; 
	
	// Solve ||d||
	float dMag = d.Magnitude();
	
	// If ||d|| < r, then q is inside the sphere
	if (dMag < SphereRadius)
	{
		GEngine->AddOnScreenDebugMessage(
			  3, 30.f, FColor::Red
			, TEXT("\nFailure! Point q is inside circle!, GeometricTestLibrary.h:452"));

		return Result;
	}

	// Solve ||d|| - r
	float numerator = dMag - SphereRadius;

	// Solve d * (||d|| - r / ||d||)
	Vector3<T> rhs = d * (numerator / dMag);
	
	// Solve q`
	Result = Point + rhs;
		
	// Print q`
	GEngine->AddOnScreenDebugMessage(4, 30.f, FColor::Red, TEXT("\nClosest Point to q is..."));
	GEngine->AddOnScreenDebugMessage(
		  5, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
		
	// Return closest point on a circle or sphere
	return Result;
}

/*
 * Intersection of an AABB and a plane
 *
 * Common for collision detection and ray tracing
 *
 * <0	Box is completely on the back side of the plane
 * >0	Box is completely on the front side of the plane
 * 0	Box intersects the plane
*/
template<typename T>
inline static int UGeometricTestLibrary::DoesAABBPlaneIntersect(
	  const Vector3<T>& PlaneNormal
	, float D						// the D from Plane
	, const AABB<T>& AABB_Ref) 
{
	// Inspect the normal and solve the min and max D values
	float MinD = 0.f;
	float MaxD = 0.f;
		
	if (PlaneNormal.X > 0.0f) 
	{
		MinD = PlaneNormal.X * AABB_Ref.GetMin().X;
		MaxD = PlaneNormal.X * AABB_Ref.GetMax().X;
	}
	else
	{
		MinD = PlaneNormal.X * AABB_Ref.GetMax().X;
		MaxD = PlaneNormal.X * AABB_Ref.GetMin().X;
	}
	
	if (PlaneNormal.Y > 0.0f) 
	{
		MinD += PlaneNormal.Y * AABB_Ref.GetMin().Y;
		MaxD += PlaneNormal.Y * AABB_Ref.GetMax().Y;
	}
	else
	{
		MinD += PlaneNormal.Y * AABB_Ref.GetMax().Y;
		MaxD += PlaneNormal.Y * AABB_Ref.GetMin().Y;
	}
	
	if (PlaneNormal.Z > 0.0f) 
	{
		MinD += PlaneNormal.Z * AABB_Ref.GetMin().Z;
		MaxD += PlaneNormal.Z * AABB_Ref.GetMax().Z;
	}
	else
	{
		MinD += PlaneNormal.Z * AABB_Ref.GetMax().Z;
		MaxD += PlaneNormal.Z * AABB_Ref.GetMin().Z;
	}
	
	// Check if completely on the front side of the plane
	if (MinD >= D) return 1;
	
	// Check if completely on the back side of the plane
	if (MaxD <= D) return -1;
	
	// Return intersection of AABB and plane
	return 0;
}		

/*
 * Intersection of two lines in 2D
 *
 * Given lines defined in 2D by:
 *
 * d1 = a1*x + b1*y
 * d2 = a2*x + b2*y
 *
 * We can solve for the point of intersection.
 *
 * x = b2*d1 - b1*d2 / a1*b2 - a2*b1
 * y = a1*d2 - a2*d1 / a1*b2 - a2*b1
*/
template<typename T>
inline static Vector2<T> UGeometricTestLibrary::DoLinesIntersect(
	  const Vector2<T>& Origin1
	, const Vector2<T>& Delta1
	, const Vector2<T>& Origin2
	, const Vector2<T>& Delta2) 
{
	Vector2<T> Result = Vector2<T>::ZeroVector;
	
	// Line AB
	// d1 = a1*x + b1*y
	const float a1 = Delta1.Y - Origin1.Y;
	const float b1 = Origin1.X - Delta1.X;
	const float d1 = a1*(Origin1.X) + b1*(Origin1.Y);
	
	// Line CD
	// d2 = a2*x + b2*y
	const float a2 = Delta2.Y - Origin2.Y;
	const float b2 = Origin2.X - Delta2.X;
	const float d2 = a2*(Origin2.X) + b2*(Origin2.Y);
		
	// Solve a1*b2 - a2*b1
	const float Determinant = (a1*b2) - (a2*b1);
	
	// Solve b2*d1 - b1*d2
	const float NumeratorX = ((b2*d1) - (b1*d2));
	
	// Solve a1*d2 - a2*d1
	const float NumeratorY = ((a1*d2) - (a2*d1));
	
	// If the lines are coincident, there are an infinite number of solutions
	if (NumeratorX == 0.f && NumeratorY == 0.f && Determinant == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nLines are coincident! Infinte solutions!, GeometricTestLibrary.h:587"));

		return Result;
	}
	
	// If the denominator of both equations is zero, there are no solutions and no intersection
	if (Determinant == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  2, 30.f, FColor::Red
			, TEXT("\nLines are parallel! No solutions and no intersection!, GeometricTestLibrary.h:597"));

		return Result;
	}
	else
	{
		// Solve x = b2*d1 - b1*d2 / a1*b2 - a2*b1
		const float x = NumeratorX / Determinant;
		
		// Solve y = a1*d2 - a2*d1 / a1*b2 - a2*b1
		const float y = NumeratorY / Determinant;
		 
		// Check if the x and y coordinates are within both lines
		if (x < std::min(Origin1.X, Delta1.X) || x > std::max(Origin1.X, Delta1.X) ||
			x < std::min(Origin2.X, Delta2.X) || x > std::max(Origin2.X, Delta2.X))
		{
			GEngine->AddOnScreenDebugMessage(
				  2, 30.f, FColor::Red
				, TEXT("\nFAILURE! X coordinates outside of lines!, GeometricTestLibrary.h:615\n"));

			return Result;
		}

		if (y < std::min(Origin1.Y, Delta1.Y) || y > std::max(Origin1.Y, Delta1.Y) ||
			y < std::min(Origin2.Y, Delta2.Y) || y > std::max(Origin2.Y, Delta2.Y))
		{
			GEngine->AddOnScreenDebugMessage(
				  2, 30.f, FColor::Red
				, TEXT("\nFAILURE! Y coordinates outside of lines!, GeometricTestLibrary.h:625\n"));

			return Result;
		}

		Result.Set((T)x, (T)y);
	}
	
	// Print Result
	GEngine->AddOnScreenDebugMessage(3, 30.f, FColor::Red, TEXT("\nPoint of intersection of two lines is..."));
	GEngine->AddOnScreenDebugMessage(
		  4, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f\n")
		, Result.X, Result.Y));

	// Returns the point of intersection of two lines
	return Result;
}

/*
 * Intersection of a Ray and an AABB
 *
 * Common for line of sight or tests that require
 * trivial rejection on complex objects.
 *
 * Determine which side of the box will be intersected
 * and then perform ray-plane intersection test on 
 * that side. If the point of intersection with the
 * plane is within AABB, there is intersection. 
 * Otherwise, there is not.
 *
 * 1, -1	No Intersection
 * 0		Intersection
*/
template<typename T>
inline static float UGeometricTestLibrary::DoesRayAABBIntersect(
	  const Vector3<T>& RayOrigin		// origin of ray
	, const Vector3<T>& RayDelta		// length and direction of ray
	, const AABB<T>& AABB_Ref			// AABB data
	, Vector3<T>* ReturnNormal)			// normal to return
{
	const float NoIntersection = FLT_MAX;	// Bogus value
	
	// Check for point inside box, trivial reject,
	// and determine parametric distance to each front face
	bool IsInside = true;
	
	float xt = 0.f, xn = 0.f;
	
	if (RayOrigin.X < AABB_Ref.GetMin().X)
	{
		xt = AABB_Ref.GetMin().X - RayOrigin.X;
		
		if (xt > RayDelta.X) return NoIntersection;
		
		xt /= RayDelta.X;
		IsInside = false;
		xn = -1.0f;
	}
	else if (RayOrigin.X > AABB_Ref.GetMax().X)
	{
		xt = AABB_Ref.GetMax().X - RayOrigin.X;
		
		if (xt < RayDelta.X) return NoIntersection;
		
		xt /= RayDelta.X;
		IsInside = false;
		xn = 1.0f;
	}
	else
		xt = -1.0f;
		
	float yt = 0.f, yn = 0.f;
	
	if (RayOrigin.Y < AABB_Ref.GetMin().Y)
	{
		yt = AABB_Ref.GetMin().Y - RayOrigin.Y;
		
		if (yt > RayDelta.Y) return NoIntersection;
		
		yt /= RayDelta.Y;
		IsInside = false;
		yn = -1.0f;
	}
	else if (RayOrigin.Y > AABB_Ref.GetMax().Y)
	{
		yt = AABB_Ref.GetMax().Y - RayOrigin.Y;
		
		if (yt < RayDelta.Y) return NoIntersection;
		
		yt /= RayDelta.Y;
		IsInside = false;
		yn = 1.0f;
	}
	else
		yt = -1.0f;
		
	float zt = 0.f, zn = 0.f;
	
	if (RayOrigin.Z < AABB_Ref.GetMin().Z)
	{
		zt = AABB_Ref.GetMin().Z - RayOrigin.Z;
		
		if (zt > RayDelta.Z) return NoIntersection;
		
		zt /= RayDelta.Z;
		IsInside = false;
		zn = -1.0f;
	}
	else if (RayOrigin.Z > AABB_Ref.GetMax().Z)
	{
		zt = AABB_Ref.GetMax().Z - RayOrigin.Z;
		
		if (zt < RayDelta.Z) return NoIntersection;
		
		zt /= RayDelta.Z;
		IsInside = false;
		zn = 1.0f;
	}
	else
		zt = -1.0f;
		
	// If ray origin is inside box, then return intersection
	if (IsInside)
	{
		if (ReturnNormal != nullptr)
			*ReturnNormal = -RayDelta;	// unnormalized
				
		return 0.0f;
	}
	
	// Select the farthest plane, the plane of intersection
	int PlaneID = 0;
	float t = xt;
	
	if (yt > t)
	{
		PlaneID = 1;
		t = yt;
	}
	
	if (zt > t)
	{
		PlaneID = 2;
		t = zt;
	}
	
	switch (PlaneID)
	{
		case 0:		// intersect with yz plane
		{
			float y = RayOrigin.Y + RayDelta.Y*t;
			if (y < AABB_Ref.GetMin().Y || y > AABB_Ref.GetMax().Y) return NoIntersection;
			
			float z = RayOrigin.Z + RayDelta.Z*t;
			if (z < AABB_Ref.GetMin().Z || z > AABB_Ref.GetMax().Z) return NoIntersection;
			
			if (ReturnNormal != nullptr)
				ReturnNormal->Set(xn, 0.0f, 0.0f);
		} break;
			
		case 1:		// intersect with xz plane
		{
			float x = RayOrigin.X + RayDelta.X*t;
			if (x < AABB_Ref.GetMin().X || x > AABB_Ref.GetMax().X) return NoIntersection;
			
			float z = RayOrigin.Z + RayDelta.Z*t;
			if (z < AABB_Ref.GetMin().Z || z > AABB_Ref.GetMax().Z) return NoIntersection;
			
			if (ReturnNormal != nullptr)
				ReturnNormal->Set(0.0f, yn, 0.0f);
		} break;
		
		case 2:		// intersect with xy plane
		{
			float x = RayOrigin.X + RayDelta.X*t;
			if (x < AABB_Ref.GetMin().X || x > AABB_Ref.GetMax().X) return NoIntersection;
			
			float y = RayOrigin.Y + RayDelta.Y*t;
			if (y < AABB_Ref.GetMin().Y || y > AABB_Ref.GetMax().Y) return NoIntersection;
			
			if (ReturnNormal != nullptr)
				ReturnNormal->Set(0.0f, 0.0f, zn);
		} break;
	}
	
	// Return parametric point of intersection
	return t;	
}

/*
 * Intersection of a ray and a plane
 *
 * Solve for t at the point of intersection:
 *
 * t = d - p0*n / ň*n
 *
 * Where:
 *
 * p0 = ray origin vector
 * n = surface normal of plane, i.e., cross product of two non-parallel edges of a polygon
 * ň = normalized ray delta vector, a.k.a. unit vector
 * d = p*n, the plane equation, the distance from the plane to a point
*/
template<typename T>
inline static bool UGeometricTestLibrary::DoesRayPlaneIntersect(
	  const Vector3<T>& RayOrigin
	, const Vector3<T>& RayDelta
	, const Vector3<T>& SurfaceNormal	// must be normalized
	, float PlaneD)
{
	bool Result = false;

	// Solve normalized surface normal
	Vector3<T> NormalizedSurfaceNormal = MyMathLibrary::Normalize(SurfaceNormal);

	// Solve d - p0
	Vector3<T> Diff = RayOrigin - PlaneD;

	// Solve ň
	Vector3<T> DeltaNormal = MyMathLibrary::Normalize(RayDelta);
	
	// Solve numerator (d - p0*n)
	float numerator = -MyMathLibrary::DotProduct(Diff, NormalizedSurfaceNormal);
					
	// Solve ň*n
	float denominator = MyMathLibrary::DotProduct(DeltaNormal, NormalizedSurfaceNormal);
	
	// If denominator is zero, then ray is parallel to the plane
	// and there is no intersection
	if (denominator == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nRay parallel to plane! No intersection!, GeometricTestLibrary.h:859\n"));

		return Result;
	}
	
	// Solve d - p0*n / ň*n
	float t = numerator / denominator;

	// Solve R, the parametric ray
	Vector3<T> R = RayOrigin + (DeltaNormal*t);
	
	// TODO: Allow intersection only with front of plane (i.e., denominator < 0)
	// Thus, intersection only if the ray points in a direction opposite 
	// to the normal of the plane.

	// Return true if t is within range
	// If t < 0 or t > length of ray, intersection does not occur
	if (t < 0 || t > R.Magnitude())
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nFAILURE! No intersection!, GeometricTestLibrary.h:880\n"));

		Result = false;
	}
	else
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nSUCCESS! Ray and Plane intersect!, GeometricTestLibrary.h:888\n"));

		Result = true;
	}

	return Result;
}

/*
 * Intersection of a ray and a sphere
 *
 * Solve for t at the point of intersection:
 *
 * t = a - sqrt(r^2 - e^2 + a^2)
 *
 * Where:
 *
 * a = e * Ď, the length of this vector
 * Ď = normalized ray delta vector, a.k.a. unit vector
 * p0 = ray origin
 * r = radius of sphere
 * c = center of sphere
 * e = c - p0
*/
template<typename T>	
inline static bool UGeometricTestLibrary::DoesRaySphereIntersect(
	  const Vector3<T>& SphereCenter	
	, float SphereRadius					
	, const Vector3<T>& RayOrigin
	, const Vector3<T>& RayDelta)
{
	bool Result = false;
	
	// Solve Ď
	Vector3<T> D = MyMathLibrary::Normalize(RayDelta);
	
	// Solve e, e^2
	Vector3<T> e = SphereCenter - RayOrigin;
	float e2 = MyMathLibrary::DotProduct(e, e);
	
	// Solve r^2
	float r2 = pow(SphereRadius, 2);
	
	// if e^2 < r^2, the ray origin is inside the sphere
	if (e2 < r2)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("Intersection at point of ray origin!, GeometricTestLibrary.h:936"));

		Result = true;
		return Result;
	}
	
	// Solve a, a^2
	float a = MyMathLibrary::DotProduct(e, D);
	float a2 = pow(a, 2);
	
	// Solve sqrt(r^2 - e^2 + a^2)
	float SqrtVal = (r2 - e2) + a2;
	float f = sqrt(SqrtVal);
	
	// If the argument square root is negative,
	// then the ray does not intersect.
	if (SqrtVal < 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  2, 30.f, FColor::Red
			, TEXT("\nFAILURE! No intersection!, Square root value is negative!, GeometricTestLibrary.h:956\n"));

		return Result;
	}
	
	// Solve a - sqrt(r^2 - e^2 + a^2)
	float t = a - f;

	// Solve R, the parametric ray
	Vector3<T> R = RayOrigin + (D*t);
	
	// If t < 0 or t > length of ray, intersection does not occur
	if (t < 0.f || t > R.Magnitude())
	{
		GEngine->AddOnScreenDebugMessage(
			  3, 30.f, FColor::Red
			, TEXT("\nFAILURE! No intersection!, GeometricTestLibrary.h:972\n"));

		Result = false;
	}
	else
	{
		GEngine->AddOnScreenDebugMessage(
			  4, 30.f, FColor::Red
			, TEXT("\nSUCCESS! Ray and sphere intersected at t = 0!, GeometricTestLibrary.h:980\n"));

		Result = true;
	}
		
	return Result; 	
}

/* 
 * Intersection of a Ray and a Triangle
 *
 * Solve the point where the ray intersects the 
 * plane containing the triangle. Then test to 
 * see whether point is inside triangle by 
 * solving barycentric coordinates of the point.
 *
 * We are testing collisions where the ray 
 * approaches the triangle from the front side.
*/
template<typename T>	
inline static float UGeometricTestLibrary::DoesRayTriangleIntersect(
	  const Vector3<T>& RayOrigin		// origin of ray
	, const Vector3<T>& RayDelta		// direction and length of ray
	, const Vector3<T>& Vertex1			// triangle vertices
	, const Vector3<T>& Vertex2			// ...
	, const Vector3<T>& Vertex3			// ...
	, float MinT)						// closest intersection found so far, start at 1.0f
{
	const float NoIntersection = FLT_MAX;
	
	// Solve clockwise edge vectors
	Vector3<T> e1 = Vertex2 - Vertex1;
	Vector3<T> e2 = Vertex3 - Vertex2;
	
	// Solve Surface Normal (unnormalized), i.e., the cross product of two non-parallel edges
	Vector3<T> n = MyMathLibrary::CrossProduct(e1, e2);
	
	// Solve gradient which tells us how steep of an 
	// angle we are approaching the front side of the triangle
	float dot = MyMathLibrary::DotProduct(n, RayDelta);
	
	// Check for a ray that is parallel to the triangle, 
	// or not pointing towards the front face of the triangle.
	// This will also reject degenerate triangles and rays as well.
	if (!(dot < 0.0f))
		return NoIntersection;
	
	// Solve d value for the plane equation
	float d = MyMathLibrary::DotProduct(n, Vertex1);
	
	// Solve parametric point of intersection with the plane
	// containing the triangle, checking at the earliest
	// possible stages for trivial rejection
	float t = d - MyMathLibrary::DotProduct(n, RayOrigin);
	
	// Is ray origin on the backside of the polygon?
	if (!(t <= 0.0f))
		return NoIntersection;
	
	// Closer intersection already found?
	if (!(t >= dot * MinT))
		return NoIntersection;
	
	// Ray intersects the plane
	t /= dot;
	check(t >= 0.0f && t <= MinT);	// "check" is from UE4. In C++, use "assert".
	
	// Solve p, the 3D point of intersection
	Vector3<T> p = RayOrigin + (RayDelta*t);
	
	// Dominant axis to select which plane to project onto
	float u0 = 0.f, u1 = 0.f, u2 = 0.f;
	float v0 = 0.f, v1 = 0.f, v2 = 0.f;
	
	if (std::fabs(n.X) > std::fabs(n.Y))
	{
		if (std::fabs(n.X) > std::fabs(n.Z))
		{
			u0 = p.Y - Vertex1.Y;
			u1 = Vertex2.Y - Vertex1.Y;
			u2 = Vertex3.Y - Vertex1.Y;
			
			v0 = p.Z - Vertex1.Z;
			v1 = Vertex2.Z - Vertex1.Z;
			v2 = Vertex3.Z - Vertex1.Z;
		}
		else
		{
			u0 = p.X - Vertex1.X;
			u1 = Vertex2.X - Vertex1.X;
			u2 = Vertex3.X - Vertex1.X;
			
			v0 = p.Y - Vertex1.Y;
			v1 = Vertex2.Y - Vertex1.Y;
			v2 = Vertex3.Y - Vertex1.Y;
		}
	}
	else
	{
		if (std::fabs(n.Y) > std::fabs(n.Z))
		{
			u0 = p.X - Vertex1.X;
			u1 = Vertex2.X - Vertex1.X;
			u2 = Vertex3.X - Vertex1.X;
			
			v0 = p.Z - Vertex1.Z;
			v1 = Vertex2.Z - Vertex1.Z;
			v2 = Vertex3.Z - Vertex1.Z;
		}
		else
		{
			u0 = p.X - Vertex1.X;
			u1 = Vertex2.X - Vertex1.X;
			u2 = Vertex3.X - Vertex1.X;
			
			v0 = p.Y - Vertex1.Y;
			v1 = Vertex2.Y - Vertex1.Y;
			v2 = Vertex3.Y - Vertex1.Y;
		}
	}
	
	// Solve Denominator
	float Denominator = u1*v2 - v1*u2;
	if (!(Denominator != 0.0f))	return NoIntersection;
	Denominator = 1.0f / Denominator;

	// Solve barycentric coordinates, checking for out of range at each step
	float alpha = (u0*v2 - v0*u2) * Denominator;
	if (!(alpha >= 0.0f)) return NoIntersection;
	
	float beta = (u1*v0 - v1*u0) * Denominator;
	if (!(beta >= 0.0f)) return NoIntersection;
	
	float gamma = 1.0f - alpha - beta;
	if (!(gamma >= 0.0f)) return NoIntersection;
	
	// Return parametric point of intersection
	return t;	
}

/*
 * Intersection of two rays in 3D
 *
 * Given two rays in 3D defined parametrically by:
 *
 * r1(t1) = p1 + t1*d1
 * r2(t2) = p2 + t1*d2
 *
 * We can solve for the point of intersection.
 *
 * t1 = ((p2 - p1) x d2) * (d1 x d2) / ||d1 x d2||^2
 * t2 = ((p2 - p1) x d1) * (d1 x d2) / ||d1 x d2||^2
 *
 * Where:
 * 
 * x = cross product
 * p1, p2 = origin vectors
 * d1, d2 = delta vectors
 * ||d1 x d2||^2 = the magnitude of the cross product of the delta vectors, squared
 *
 * NOTE: the range of t1 and t2 is not restricted.
 *
 * There will be 1 solution, no solutions, or infinite
 * solutions (coincident lines). 3D has another case called "skew",
 * meaning the lines do not share a common plane.
*/
template<typename T>
inline static bool UGeometricTestLibrary::DoRaysIntersect(
	  const Vector3<T>& RayOrigin1				// p1
	, const Vector3<T>& RayDelta1				// d1
	, const Vector3<T>& RayOrigin2				// p2
	, const Vector3<T>& RayDelta2)				// d2
{
	// Solve p2 - p1
	Vector3<T> OriginDistance = RayOrigin2 - RayOrigin1;
	
	// Solve ((p2 - p1) x d2)
	Vector3<T> CrossProdT1 = MyMathLibrary::CrossProduct(OriginDistance, RayDelta2);
		
	// Solve ((p2 - p1) x d1)
	Vector3<T> CrossProdT2 = MyMathLibrary::CrossProduct(OriginDistance, RayDelta1);
	
	// Solve ||d1 x d2||^2
	Vector3<T> CrossProdDirs = MyMathLibrary::CrossProduct(RayDelta1, RayDelta2);
	float CrossProdMagnitude = CrossProdDirs.Magnitude();
	float Denominator = pow(CrossProdMagnitude, 2);
	
	// Solve ((p2 - p1) x d2) * (d1 x d2)
	float Numerator1 = MyMathLibrary::DotProduct(CrossProdT1, CrossProdDirs);
	
	// Solve ((p2 - p1) x d1) * (d1 x d2)
	float Numerator2 = MyMathLibrary::DotProduct(CrossProdT2, CrossProdDirs);
	
	// If the lines are coincident, there are an infinite number of solutions
	// All numerators and denominators are zero in this case
	if (Numerator1 == 0.f && Numerator2 == 0.f && Denominator == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nRays are coincident! Infinite number of solutions!, GeometricTestLibrary.h:1179"));
		
		return false;
	}
	
	// If the lines are parallel, then the cross product of d1 and d2 is the zero vector
	// Therefore, the denominator of both equations is zero
	if (Denominator == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  2, 30.f, FColor::Red
			, TEXT("\nRays are parallel! No solutions and no intersection!, GeometricTestLibrary.h:1190"));

		return false;
	}
		
	// Solve ((p2 - p1) x d2) * (d1 x d2) / ||d1 x d2||^2
	float t1 = Numerator1 / Denominator;
	
	// Solve ((p2 - p1) x d1) * (d1 x d2) / ||d1 x d2||^2
	float t2 = Numerator2 / Denominator;
	
	// Solve R1, R2 for point of intersection
	Vector3<T> R1 = RayOrigin1 + (RayDelta1*t1);
	Vector3<T> R2 = RayOrigin2 + (RayDelta2*t2);

	// Return true if the two rays intersect at t1 == 0 && t2 == 0
	if (t1 < 0 || t1 > R1.Magnitude() || t2 < 0 || t2 > R2.Magnitude())
	{
		GEngine->AddOnScreenDebugMessage(
			2, 30.f, FColor::Red
			, TEXT("\nRays are skew! They are neither parallel nor do they intersect!, GeometricTestLibrary.h:1210"));

		return false;
	}
	else
	{
		GEngine->AddOnScreenDebugMessage(
			2, 30.f, FColor::Red
			, TEXT("\nSUCCESS! Rays intersect!, GeometricTestLibrary.h:1218"));

		return true;
	}
}

/*
 * Dynamic (moving) sphere and static (stationary) plane intersection:
 * 
 * Solve for t, the point of first intersection
 *
 * t = d - c*ň + r \ Ď*ň
 *
 * Where:
 *
 * d = p*ň, the plane equation, the distance from the plane to a point
 * c = center of sphere
 * r = radius of sphere
 * Ď = normalized delta vector, a.k.a. unit vector
 * ň = normalized surface normal of plane
 * c - r*ň = point of contact
 * c + t*Ď = motion of the center of the sphere
 *
 * NOTE: We are flipping d - c, using the inverse dot product to get signed distance
*/
template<typename T>
inline static bool UGeometricTestLibrary::DoesSpherePlaneIntersect_Dynamic(
	  const Vector3<T>& PlaneNormal			// must be normalized first
	, float PlaneD
	, const Vector3<T>& SphereDeltaVector
	, const Vector3<T>& SphereCenter
	, float SphereRadius)
{
	// Solve normalized plane normal
	Vector3<T> NormalizedPlaneNormal = MyMathLibrary::Normalize(PlaneNormal);

	// Solve d - c
	Vector3<T> Diff = SphereCenter - PlaneD;

	// Solve d - c*ň, flipped
	float Dot = -MyMathLibrary::DotProduct(Diff, NormalizedPlaneNormal);

	// Solve Ď
	Vector3<T> DeltaNormal = MyMathLibrary::Normalize(SphereDeltaVector);

	// Solve Ď*ň
	float Denominator = MyMathLibrary::DotProduct(DeltaNormal, NormalizedPlaneNormal);

	// If Denominator is zero, there is no intersection
	if (Denominator == 0.f)
		return false;

	// Solve t = d - c*ň + r \ Ď*ň
	float t = (Dot + SphereRadius) / Denominator;

	// Solve l, the motion of the sphere
	Vector3<T> l = SphereCenter + (DeltaNormal*t);

	// If t < 0 or t > length of total relative displacement of sphere, intersection does not occur
	if (t < 0 || t > l.Magnitude())
	{
		GEngine->AddOnScreenDebugMessage(
			  -1, 30.f, FColor::Red
			, TEXT("\nNo intersection!, GeometricTestLibrary.h:1281\n"));

		return false;
	}
	else
	{
		GEngine->AddOnScreenDebugMessage(
			  -1, 30.f, FColor::Red
			, TEXT("\nSUCCESS! Sphere and plane intersect during time in question!, GeometricTestLibrary.h:1289\n"));

		return true;
	}
}

/*
 * Static intersection of a sphere and a plane
 *
 * Given a sphere and plane, determine which side of 
 * the plane the sphere is on.
 *
 * <0	Sphere is completely on the back
 * >0	Sphere is completely on the front
 * 0 	Sphere straddles plane
*/
template<typename T>
inline static int UGeometricTestLibrary::DoesSpherePlaneIntersect_Static(
	  const Vector3<T>& PlaneNormal		// must be normalized first
	, float PlaneD						// p*PlaneNormal = PlaneD
	, const Vector3<T>& SphereCenter	// center of sphere
	, float SphereRadius)				// radius of sphere 
{
	/*
	 * Solved by calculating distance from center
	 * of sphere to the plane:
	 *
	 * d = n*c - (p*n)
	*/

	// Solve normalized plane normal
	Vector3<T> NormalizedPlaneNormal = MyMathLibrary::Normalize(PlaneNormal);

	// Solve d
	float Dot = MyMathLibrary::DotProduct(NormalizedPlaneNormal, SphereCenter);
	float d = Dot - PlaneD;
	
	// On the front side?
	if (d >= SphereRadius)
		return 1;
	
	// On the back side?
	if (d <= SphereRadius)
		return -1;
	
	// Return intersection of sphere and the plane
	return 0;	
}

/*
 * Dynamic intersection of two circles or spheres
 *
 * Solve for t, the point of first intersection:
 *
 * t = e*Ď - sqrt((e*Ď)^2 + r^2 - e*e)
 *
 * Where:
 *
 * Ď = normalized delta vector 
 * D = moving delta vector - stationary delta vector
 * e = stationary sphere center - moving sphere center 
 * ||e|| = magnitude of e
 * r = sum of the sphere radii
 *
 * If ||e|| < r, then the spheres are intersecting at t = 0
 * If t < 0 or t > l, intersection does not occur
 * If the square root argument is negative, there is no intersection
*/
template<typename T>
inline static bool UGeometricTestLibrary::DoSpheresIntersect_Dynamic(
	  const Vector3<T>& StationaryDeltaVector
	, const Vector3<T>& MovingDeltaVector
	, const Vector3<T>& StationarySphereCenter
	, const Vector3<T>& MovingSphereCenter		// Defined at t = 0
	, float StationarySphereRadius
	, float MovingSphereRadius) 
{
	bool Result = false;
	
	// Solve D, Ď
	Vector3<T> D = MovingDeltaVector - StationaryDeltaVector;
	Vector3<T> NormalizedD = MyMathLibrary::Normalize(D);
	
	// Solve e, ||e||
	Vector3<T> e = StationarySphereCenter - MovingSphereCenter;
	float eMag = e.Magnitude();
	
	// Solve r, r^2
	float r = StationarySphereRadius + MovingSphereRadius;
	float RadiusSquared = pow(r, 2);
	
	// Solve e*Ď, (e*Ď)^2
	float LHS = MyMathLibrary::DotProduct(e, NormalizedD);
	float DotSquared = pow(LHS, 2);
	
	// Solve e*e
	float eDot = MyMathLibrary::DotProduct(e, e);
	
	// Solve sqrt((e * Ď)^2 + r^2 - e*e)
	float SqrtVal = (DotSquared + RadiusSquared) - eDot;
	float RHS = sqrt(SqrtVal);
	
	// If the square root argument is negative, there is no intersection
	if (SqrtVal < 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nFAILURE! No intersection! Square root value is negative!, GeometricTestLibrary.h:1396\n"));

		Result = false;
		return Result;
	}
	
	// Solve t
	float t = LHS - RHS;

	// Solve l, the position of the center of moving sphere at time t
	Vector3<T> l = MovingSphereCenter + (NormalizedD*t);
	
	// If t < 0 or t > length of total relative displacement of spheres, intersection does not occur
	if (t < 0.f || t > l.Magnitude())
	{
		GEngine->AddOnScreenDebugMessage(
			  2, 30.f, FColor::Red
			, TEXT("\nFAILURE! No intersection!, GeometricTestLibrary.h:1413\n"));

		Result = false;
	}
	else
	{
		// If ||e|| < r, then the spheres are intersecting at t = 0 during the time in question
		if (eMag < r)
		{
			GEngine->AddOnScreenDebugMessage(
				  3, 30.f, FColor::Red
				, TEXT("\nSUCCESS! Spheres intersected at t = 0!, GeometricTestLibrary.h:1424\n"));

			Result = true;
		}
	}
	
	return Result;	
}

/*
 * Static intersection of two circles or spheres
 *
 * Check if d^2 < (r1 + r2)^2
 *
 * Where:
 *
 * d = distance between the sphere centers
 * r1, r2 = radii of spheres
*/
template<typename T>
inline static bool UGeometricTestLibrary::DoSpheresIntersect_Static(
	  const Vector3<T>& SphereCenter1
	, float SphereRadius1
	, const Vector3<T>& SphereCenter2
	, float SphereRadius2) 
{
	// Solve d
	float Distance = MyMathLibrary::Distance(SphereCenter1, SphereCenter2);
	
	// Solve squares
	float DistanceSquared = pow(Distance, 2);
	float SumRadii = SphereRadius1 + SphereRadius2;
	float RadiiSquared = pow(SumRadii, 2);
	
	// Solve d^2 < (r1 + r2)^2
	if (DistanceSquared < RadiiSquared)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nSUCCESS! Spheres intersect!, GeometricTestLibrary.h:1463\n"));

		return true;
	}
	else
	{
		GEngine->AddOnScreenDebugMessage(
			  2, 30.f, FColor::Red
			, TEXT("\nFAILURE! Spheres do not intersect!, GeometricTestLibrary.h:1471\n"));

		return false;
	}	
}

/*
 * Intersection of three planes
 *
 * Solve for t, the point of intersection:
 *
 * t = d1*(n2 x n3) + d2*(n3 x n1) + d3*(n1 x n2) / (n1 x n2) * n3
 *
 * Where:
 *
 * x = cross product
 * n = plane normal, a.k.a. surface normal, a.k.a. cross product of two non-parallel edges of polygon
 * d = p*n = The Plane Equation (ax + by + cz + d = 0), the distance from the plane to the point (i.e., d1 = p*n1) 
 *
 * If any pair of planes is parallel, the point of intersection
 * either does not exist or is not unique.  The denominator of 
 * the equation will be zero in this case.
*/
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::PointOfIntersection(
	  const Vector3<T>& PlaneNormal1	// must be normalized first
	, float PlaneD1						// p*PlaneNormal = PlaneD1
	, const Vector3<T>& PlaneNormal2	// must be normalized first
	, float PlaneD2						// p*PlaneNormal2 = PlaneD2	
	, const Vector3<T>& PlaneNormal3	// must be normalized first
	, float PlaneD3)					// p*PlaneNormal3 = PlaneD3	
{
	Vector3<T> Result = Vector3<T>::ZeroVector;

	// Normalized Plane Normals
	Vector3<T> NormalizedPlaneNormal1 = MyMathLibrary::Normalize(PlaneNormal1);
	Vector3<T> NormalizedPlaneNormal2 = MyMathLibrary::Normalize(PlaneNormal2);
	Vector3<T> NormalizedPlaneNormal3 = MyMathLibrary::Normalize(PlaneNormal3);
	
	// Solve Cross Products
	Vector3<T> CrossProd1 = MyMathLibrary::CrossProduct(NormalizedPlaneNormal1, NormalizedPlaneNormal2);
	Vector3<T> CrossProd2 = MyMathLibrary::CrossProduct(NormalizedPlaneNormal2, NormalizedPlaneNormal3);
	Vector3<T> CrossProd3 = MyMathLibrary::CrossProduct(NormalizedPlaneNormal3, NormalizedPlaneNormal1);
	
	// Solve d1*(n2 x n3)
	CrossProd2 *= PlaneD1;

	// Solve d2*(n3 x n1)
	CrossProd3 *= PlaneD2;
	
	// Solve d3*(n1 x n2)
	CrossProd1 *= PlaneD3;

	// Solve Numerator
	Vector3<T> Numerator = CrossProd2 + CrossProd3 + CrossProd1;
	
	// Solve (n1 x n2) * n3
	float Denominator = MyMathLibrary::DotProduct(CrossProd1, NormalizedPlaneNormal3);
	
	// If any pair of planes is parallel, the point of intersection
	// either does not exist or is not unique. The denominator of 
	// the equation will be zero in this case.
	if (Denominator == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nFAILURE! A pair of planes are parallel! No intersection!, GeometricTestLibrary.h:1537"));

		return Result;
	}

	// Solve t
	Result = Numerator / Denominator;

	// Print Result
	GEngine->AddOnScreenDebugMessage(2, 30.f, FColor::Red, TEXT("\nPoint of intersection of three planes is..."));
	GEngine->AddOnScreenDebugMessage(
		  3, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
	
	// Return the point of intersection
	return Result;
}

// Intersection of Two AABBs
template<typename T>
inline static bool UGeometricTestLibrary::DoAABBsIntersect(
	  const AABB<T>& A
	, const AABB<T>& B) 
{
	/*
     * Separating Axis Test
	 *
	 * If there is no overlap on a particular axis,
	 * then the two AABBs do not intersect
	*/

	if (A.GetMin().X >= B.GetMax().X) return false;
	if (A.GetMax().X <= B.GetMin().X) return false;
	
	if (A.GetMin().Y >= B.GetMax().Y) return false;
	if (A.GetMax().Y <= B.GetMin().Y) return false;
	
	if (A.GetMin().Z >= B.GetMax().Z) return false;
	if (A.GetMax().Z <= B.GetMin().Z) return false;
	
	// Overlap on all three axes, so their intersection must be non-empty
	return true;	
}

/*
 * Reflection Vector ("Perfect Mirror Bounce")
 *
 * Formula:
 *
 * r = 2(n * l)n - l
 *
 * But in order to "bounce", we must negate the r:
 *
 * r = -2(n * l)n + l
 *
 * Where:
 *
 * n = normalized surface normal
 * l = normalized delta of light source
*/
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::SolveReflectionVector(
	  const Vector3<T>& SurfaceNormal
	, const Vector3<T>& LightDelta) 
{
	Vector3<T> Result = Vector3<T>::ZeroVector;

	// Solve n
	const Vector3<T> NormalizedSurfaceNormal = MyMathLibrary::Normalize(SurfaceNormal);
	
	// Solve l
	const Vector3<T> LightNormal = MyMathLibrary::Normalize(LightDelta);
	
	// Solve -2(n * l)
	const float Dot = -2 * MyMathLibrary::DotProduct(NormalizedSurfaceNormal, LightNormal);
	
	// Solve n + l
	const Vector3<T> Sum = NormalizedSurfaceNormal + LightNormal;
	
	// Solve r
	Result = Sum * Dot;
	
	// Print r
	GEngine->AddOnScreenDebugMessage(1, 30.f, FColor::Red, TEXT("\nReflection vector of l about n is..."));
	GEngine->AddOnScreenDebugMessage(
		  2, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
	
	// Return the reflection vector
	return Result;
}

/*
 * TODO: Barycentric Coordinates of a Triangle in 2D
 * 
 * Solve for (b1, b2, b3):
 *
 * b1 = (py - y3)*(x2 - x3) + (y2 - y3)*(x3 - px) / (y1 - y3)*(x2 - x3) + (y2 - y3)*(x3 - x1)
 * b2 = (py - y1)*(x3 - x1) + (y3 - y1)*(x1 - px) / (y1 - y3)*(x2 - x3) + (y2 - y3)*(x3 - x1)
 * b3 = (py - y2)*(x1 - x2) + (y1 - y2)*(x2 - px) / (y1 - y3)*(x2 - x3) + (y2 - y3)*(x3 - x1)
 *
 * Where:
 *
 * p = given point
 * x1, x2, x3 = 
 * y1, y2, y3 = 
*/

// Barycentric Coordinates of a Triangle in 3D
template<typename T>
inline static Vector3<T> UGeometricTestLibrary::SolveBarycentricCoordinates3D(
	  const Vector3<T> Vertices[3]		// vertices of the triangle
	, const Vector3<T>& Point			// p
	, float Barycentric[3])				// barycentric coordinates
{
	Vector3<T> Result = Vector3<T>::ZeroVector;
	
	// Solve two clockwise edge vectors
	Vector3<T> d1 = Vertices[1] - Vertices[0];
	Vector3<T> d2 = Vertices[2] - Vertices[1];
	
	// Solve n, the surface normal, unnormalized
	Vector3<T> n = MyMathLibrary::CrossProduct(d1, d2);
	
	// Locate dominant axis of normal and select plane of projection
	float u1, u2, u3, u4;
	float v1, v2, v3, v4;
	
	if ((std::fabs(n.X) >= std::fabs(n.Y)) && (std::fabs(n.X) >= std::fabs(n.Z)))
	{
		// Discard X and project onto yz plane
		u1 = Vertices[0].Y - Vertices[2].Y;
		u2 = Vertices[1].Y - Vertices[2].Y;
		u3 = Point.Y - Vertices[0].Y;
		u4 = Point.Y - Vertices[2].Y;
		
		v1 = Vertices[0].Z - Vertices[2].Z;
		v2 = Vertices[1].Z - Vertices[2].Z;
		v3 = Point.Z - Vertices[0].Z;
		v4 = Point.Z - Vertices[2].Z;
	}
	else if (std::fabs(n.Y) >= std::fabs(n.Z))
	{
		// Discard Y and project onto xz plane
		u1 = Vertices[0].Z - Vertices[2].Z;
		u2 = Vertices[1].Z - Vertices[2].Z;
		u3 = Point.Z - Vertices[0].Z;
		u4 = Point.Z - Vertices[2].Z;
		
		v1 = Vertices[0].X - Vertices[2].X;
		v2 = Vertices[1].X - Vertices[2].X;
		v3 = Point.X - Vertices[0].X;
		v4 = Point.X - Vertices[2].X;
	}
	else
	{
		// Discard Z and project onto xy plane
		u1 = Vertices[0].X - Vertices[2].X;
		u2 = Vertices[1].X - Vertices[2].X;
		u3 = Point.X - Vertices[0].X;
		u4 = Point.X - Vertices[2].X;
		
		v1 = Vertices[0].Y - Vertices[2].Y;
		v2 = Vertices[1].Y - Vertices[2].Y;
		v3 = Point.Y - Vertices[0].Y;
		v4 = Point.Y - Vertices[2].Y;
	}
	
	// Solve denominator
	float Denominator = v1*u2 - v2*u1;
	
	// If denominator is zero, there are no coordinates because the triangle has zero area
	if (Denominator == 0.f)
	{
		GEngine->AddOnScreenDebugMessage(
			  1, 30.f, FColor::Red
			, TEXT("\nFAILURE! No coordinates!, GeometricTestLibrary.h:1714\n"));

		return Result;
	}
	
	// Solve barycentric coordinates
	float OneOverDenom = 1 / Denominator;
	Barycentric[0] = (v4*u2 - v2*u4) * OneOverDenom;
	Barycentric[1] = (v1*u3 - v3*u1) * OneOverDenom;
	Barycentric[2] = 1.f - Barycentric[0] - Barycentric[1];
			
	Result.Set((float)Barycentric[0], (float)Barycentric[1], (float)Barycentric[2]);
	
	// Print Result
	GEngine->AddOnScreenDebugMessage(2, 30.f, FColor::Red, TEXT("\nBarycentric coordinates of triangle are..."));
	GEngine->AddOnScreenDebugMessage(
		  3, 30.f, FColor::Red
		, FString::Printf(TEXT("X: %f, Y: %f, Z: %f\n")
		, Result.X, Result.Y, Result.Z));
	
	// Return barycentric coordinates of triangle
	return Result;
}

/*
 * These previous problems are courtesy of 
 * "3D Math Primer for Graphics and Game Development" by Dunn & Parberry:
 *
 * Intersection of AABB-Plane
 * Intersection Ray-AABB
 * Intersection of Ray-Triangle
 * Intersection of Two AABBs
 * Barycentric Coordinates of a Triangle in 3D
*/