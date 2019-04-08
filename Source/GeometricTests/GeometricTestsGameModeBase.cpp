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

#include "GeometricTestsGameModeBase.h"
#include "GeometricTestLibrary.h"
#include "Math/LinearAlgebra/Vector.h"
#include "Math/Operations/MathOperations.h"
#include "Primitives/AABB.h"
#include "Primitives/Circle.h"
#include "Primitives/Plane.h"
#include "Primitives/Ray.h"
#include "Engine.h"

void AGeometricTestsGameModeBase::BeginPlay()
{
	Super::BeginPlay();
	
	// Closest Point in AABB
	//TEST_ClosestPointInAABB();
	
	// Closest Point on Ray
	//TEST_ClosestPointOnRay();
	
	// Closest Point on Plane
	//TEST_ClosestPointOnPlane();
	
	// Closest Point on Sphere
	//TEST_ClosestPointOnSphere();
	
	// Intersection AABB-Plane
	//TEST_DoesAABBPlaneIntersect();
	
	// Intersection of Two Lines in 2D
	//TEST_DoLinesIntersect();
	
	// Intersection Ray-AABB
	//TEST_DoesRayAABBIntersect();
	
	// Intersection Ray-Plane
	//TEST_DoesRayPlaneIntersect();
	
	// Intersection Ray-Sphere
	//TEST_DoesRaySphereIntersect();  
	
	// Intersection Ray-Triangle
	//TEST_DoesRayTriangleIntersect(); 
	
	// Intersection of Two Rays in 3D
	//TEST_DoRaysIntersect();
	
	// Dynamic Intersection Sphere-Plane
	TEST_DoesSpherePlaneIntersect_Dynamic();
	
	// Static Intersection Sphere-Plane
	//TEST_DoesSpherePlaneIntersect_Static();
	
	// Dynamic Intersection of Two Circles or Spheres
	//TEST_DoSpheresIntersect_Dynamic();
	
	// Static Intersection of Two Circles or Spheres
	//TEST_DoSpheresIntersect_Static();
	
	// Intersection of Three Planes
	//TEST_DoThreePlanesIntersect();
	
	// Intersection of Two AABBs
	//TEST_DoAABBsIntersect();
	
	// Reflection Vector
	//TEST_SolveReflectionVector();
	
	// Barycentric Coordinates in 3D
	//TEST_SolveBarycentricCoordinates3D();
}

void AGeometricTestsGameModeBase::TEST_ClosestPointInAABB()
{
	const Vector3<float> Point = Vector3<float>(6.f, 3.f, -4.f);
	
	const Vector3<float> MinBounds = Vector3<float>(-1.f, -1.f, -1.f);
	const Vector3<float> MaxBounds = Vector3<float>(5.f, 5.f, 5.f);
	AABB<float> MyAABB(MinBounds, MaxBounds);
	
	UGeometricTestLibrary::ClosestPointInAABB(Point, MyAABB);
}

void AGeometricTestsGameModeBase::TEST_ClosestPointOnRay()
{
	const Vector3<float> Point = Vector3<float>(12.f, 6.f, 4.f);
	//const Vector3<float> Point = Vector3<float>(5.f, 5.f, 6.f);
	
	const Vector3<float> RayOrigin = Vector3<float>(2.f, 2.f, 4.f);
	const Vector3<float> RayDelta = Vector3<float>(10.f, 15.f, 4.f);
	Ray<float> MyRay(RayOrigin, RayDelta);
		
	UGeometricTestLibrary::ClosestPointOnRay(Point, MyRay.GetOrigin(), MyRay.GetDelta());
}

void AGeometricTestsGameModeBase::TEST_ClosestPointOnPlane()
{
	const Vector3<float> Point = Vector3<float>(12.f, 6.f, 4.f);

	const float A = 2.f;
	const float B = 3.f;
	const float C = 1.f;
	const float D = 3.f;
	Plane<float> MyPlane(A, B, C, D);
		
	UGeometricTestLibrary::ClosestPointOnPlane(Point, MyPlane.Normal, MyPlane.D);
}

void AGeometricTestsGameModeBase::TEST_ClosestPointOnSphere()  
{
	const Vector3<float> Point = Vector3<float>(12.f, 6.f, 4.f);
	
	const Vector3<float> SphereCenter = Vector3<float>(0.f, 0.f, 0.f);
	const float SphereRadius = 10.f;
	Circle<float> MyCircle(SphereRadius, SphereCenter);
	
	UGeometricTestLibrary::ClosestPointOnSphere(Point, MyCircle.GetCenter(), MyCircle.GetRadius());
}

void AGeometricTestsGameModeBase::TEST_DoesAABBPlaneIntersect()  
{
	const float A = 2.f;
	const float B = 3.f;
	const float C = 1.f;
	const float D = 3.f;
	Plane<float> MyPlane(A, B, C, D);
		
	const Vector3<float> MinBounds = Vector3<float>(-1.f, -1.f, -1.f);
	const Vector3<float> MaxBounds = Vector3<float>(5.f, 5.f, 5.f);
	AABB<float> MyAABB(MinBounds, MaxBounds);
	
	UGeometricTestLibrary::DoesAABBPlaneIntersect(MyPlane.Normal, MyPlane.D, MyAABB);
}

void AGeometricTestsGameModeBase::TEST_DoLinesIntersect()  
{
	const Vector2<float> Origin1 = Vector2<float>(1.f, 1.f);
	const Vector2<float> Delta1 = Vector2<float>(7.f, 15.f);
	
	//const Vector2<float> Origin2 = Vector2<float>(-1.f, 5.f);
	//const Vector2<float> Delta2 = Vector2<float>(-1.f, 2.f);
	const Vector2<float> Origin2 = Vector2<float>(-1.f, 5.f);
	const Vector2<float> Delta2 = Vector2<float>(8.f, 12.f);
	
	UGeometricTestLibrary::DoLinesIntersect(Origin1, Delta1, Origin2, Delta2);
}

void AGeometricTestsGameModeBase::TEST_DoesRayAABBIntersect()  
{
	Vector3<float>* ReturnNormal = nullptr;

	const Vector3<float> RayOrigin = Vector3<float>(-1.f, 2.f, 5.f);
	const Vector3<float> RayDelta = Vector3<float>(-1.f, 5.f, 7.f);
	Ray<float> MyRay(RayOrigin, RayDelta);
	
	const Vector3<float> MinBounds = Vector3<float>(-1.f, -1.f, -1.f);
	const Vector3<float> MaxBounds = Vector3<float>(5.f, 5.f, 5.f);
	AABB<float> MyAABB(MinBounds, MaxBounds);
	
	UGeometricTestLibrary::DoesRayAABBIntersect(MyRay.GetOrigin(), MyRay.GetDelta(), MyAABB, ReturnNormal);

	/*
	 * NOTE: AABB can be replaced by any box, square, capsule, etc. 
	*/
}

void AGeometricTestsGameModeBase::TEST_DoesRayPlaneIntersect()  
{
	const Vector3<float> RayOrigin = Vector3<float>(1.f, 1.f, 1.f);
	const Vector3<float> RayDelta = Vector3<float>(7.f, 15.f, 18.f);
	Ray<float> MyRay(RayOrigin, RayDelta);

	const float A = -1.f;
	const float B = 2.f;
	const float C = 1.f;
	const float D = 2.f;
	Plane<float> MyPlane(A, B, C, D);
	
	UGeometricTestLibrary::DoesRayPlaneIntersect(
		  MyRay.GetOrigin()
		, MyRay.GetDelta()
		, MyPlane.Normal
		, MyPlane.D);
}

void AGeometricTestsGameModeBase::TEST_DoesRaySphereIntersect()  
{
	const Vector3<float> RayOrigin = Vector3<float>(0.f, 0.f, 0.f);
	const Vector3<float> RayDelta = Vector3<float>(5.f, 8.f, 5.f);
	Ray<float> MyRay(RayOrigin, RayDelta);

	const Vector3<float> SphereCenter = Vector3<float>(3.f, 3.f, 3.f);
	const float SphereRadius = 5.f;
	Circle<float> MyCircle(SphereRadius, SphereCenter);
	
	UGeometricTestLibrary::DoesRaySphereIntersect(
		  MyCircle.GetCenter()
		, MyCircle.GetRadius()
		, MyRay.GetOrigin()
		, MyRay.GetDelta());
}

void AGeometricTestsGameModeBase::TEST_DoesRayTriangleIntersect()  
{
	const Vector3<float> RayOrigin = Vector3<float>(1.f, 1.f, 1.f);
	const Vector3<float> RayDelta = Vector3<float>(7.f, 15.f, 18.f);
	Ray<float> MyRay(RayOrigin, RayDelta);
	
	const Vector3<float> Vertex1 = Vector3<float>(2.f, 3.f, 1.f);
	const Vector3<float> Vertex2 = Vector3<float>(3.f, 5.f, 2.f);
	const Vector3<float> Vertex3 = Vector3<float>(2.f, 2.f, 2.f);
	
	float MinT = 1.f;
	
	UGeometricTestLibrary::DoesRayTriangleIntersect(
		  MyRay.GetOrigin()
		, MyRay.GetDelta()
		, Vertex1
		, Vertex2
		, Vertex3
		, MinT);
}

void AGeometricTestsGameModeBase::TEST_DoRaysIntersect()  
{
	const Vector3<float> RayOrigin1 = Vector3<float>(10.f, 15.f, 1.f);
	const Vector3<float> RayDelta1 = Vector3<float>(49.f, 25.f, 1.f);
	
	const Vector3<float> RayOrigin2 = Vector3<float>(20.f, 10.f, 1.f);
	const Vector3<float> RayDelta2 = Vector3<float>(32.f, 32.f, 1.f);
	
	Ray<float> Ray1(RayOrigin1, RayDelta1);
	Ray<float> Ray2(RayOrigin2, RayDelta2);
	
	UGeometricTestLibrary::DoRaysIntersect(
		  Ray1.GetOrigin()
		, Ray1.GetDelta()
		, Ray2.GetOrigin()
		, Ray2.GetDelta());
}

void AGeometricTestsGameModeBase::TEST_DoesSpherePlaneIntersect_Dynamic()  
{
	const float A = 2.f;
	const float B = 3.f;
	const float C = 1.f;
	const float D = 3.f;
	Plane<float> MyPlane(A, B, C, D);
	
	const Vector3<float> SphereCenter = Vector3<float>(3.f, 3.f, 3.f);
	const float SphereRadius = 5.f;
	Circle<float> MyCircle(SphereRadius, SphereCenter);

	Vector3<float> SphereDeltaVector = Vector3<float>(2.f, 3.5f, 0.f);
	//Vector3<float> SphereDeltaVector = Vector3<float>(11.f, 3.5f, -9.f);
	
	UGeometricTestLibrary::DoesSpherePlaneIntersect_Dynamic(
		  MyPlane.Normal
		, MyPlane.D
		, SphereDeltaVector
		, MyCircle.GetCenter()
		, MyCircle.GetRadius());
}

void AGeometricTestsGameModeBase::TEST_DoesSpherePlaneIntersect_Static()  
{
	const float A = 2.f;
	const float B = 3.f;
	const float C = 1.f;
	const float D = 3.f;
	Plane<float> MyPlane(A, B, C, D);
		
	const Vector3<float> SphereCenter = Vector3<float>(-5.f, -5.f, -5.f);
	const float SphereRadius = 5.f;
	Circle<float> MyCircle(SphereRadius, SphereCenter);
	
	UGeometricTestLibrary::DoesSpherePlaneIntersect_Static(
		  MyPlane.Normal
		, MyPlane.D
		, MyCircle.GetCenter()
		, MyCircle.GetRadius());
}

void AGeometricTestsGameModeBase::TEST_DoSpheresIntersect_Dynamic()  
{
	Vector3<float> StationaryDeltaVector = Vector3<float>(10.f, 2.f, -9.f);
	Vector3<float> MovingDeltaVector = Vector3<float>(11.f, 3.5f, -9.f);
	
	const Vector3<float> StationarySphereCenter = Vector3<float>(-5.f, -5.f, -5.f);
	Vector3<float> MovingSphereCenter = Vector3<float>::ZeroVector;
	
	const float StationarySphereRadius = 5.f;
	const float MovingSphereRadius = 6.f;

	Circle<float> MovingCircle(MovingSphereRadius, MovingSphereCenter);
	Circle<float> StationaryCircle(StationarySphereRadius, StationarySphereCenter);
	
	UGeometricTestLibrary::DoSpheresIntersect_Dynamic(
		  StationaryDeltaVector
		, MovingDeltaVector
		, StationaryCircle.GetCenter()
		, MovingCircle.GetCenter()
		, StationaryCircle.GetRadius()
		, MovingCircle.GetRadius());
}

void AGeometricTestsGameModeBase::TEST_DoSpheresIntersect_Static()  
{
	const Vector3<float> SphereCenter1 = Vector3<float>::ZeroVector;
	const Vector3<float> SphereCenter2 = Vector3<float>(-5.f, -5.f, -5.f);
	
	const float SphereRadius1 = 12.f;
	const float SphereRadius2 = 2.5f;
	
	Circle<float> Circle1(SphereRadius1, SphereCenter1);
	Circle<float> Circle2(SphereRadius2, SphereCenter2);
	
	UGeometricTestLibrary::DoSpheresIntersect_Static(
		  Circle1.GetCenter()
		, Circle1.GetRadius()
		, Circle2.GetCenter()
		, Circle2.GetRadius());
}

void AGeometricTestsGameModeBase::TEST_DoThreePlanesIntersect()  
{
	float A1 = 2.f;
	float B1 = 3.f;
	float C1 = 1.f;
	float D1 = 3.f;
	Plane<float> Plane1(A1, B1, C1, D1);
	
	float A2 = -1.f;
	float B2 = 2.f;
	float C2 = 1.f;
	float D2 = 2.f;
	Plane<float> Plane2(A2, B2, C2, D2);
	
	float A3 = 6.f;
	float B3 = -5.f;
	float C3 = 2.f;
	float D3 = -8.f;
	Plane<float> Plane3(A3, B3, C3, D3);
		
	UGeometricTestLibrary::PointOfIntersection(
		  Plane1.Normal
		, Plane1.D
		, Plane2.Normal
		, Plane2.D
		, Plane3.Normal
		, Plane3.D);
}

void AGeometricTestsGameModeBase::TEST_DoAABBsIntersect()  
{
	const Vector3<float> MinBoundsA = Vector3<float>(-1.f, -1.f, -1.f);
	const Vector3<float> MaxBoundsA = Vector3<float>(5.f, 5.f, 5.f);
	AABB<float> A(MinBoundsA, MaxBoundsA);
	
	const Vector3<float> MinBoundsB = Vector3<float>(1.f, 2.f, 1.f);
	const Vector3<float> MaxBoundsB = Vector3<float>(3, 2.f, 3.f);
	AABB<float> B(MinBoundsB, MaxBoundsB);
	
	UGeometricTestLibrary::DoAABBsIntersect(A, B);
}

void AGeometricTestsGameModeBase::TEST_SolveReflectionVector()
{
	Vector3<float> LightDelta = Vector3<float>(-3.f, -5.f, 7.f);
	
	const float A = 6.f;
	const float B = -5.f;
	const float C = 2.f;
	const float D = -8.f;
	Plane<float> MyPlane(A, B, C, D);
		
	UGeometricTestLibrary::SolveReflectionVector(MyPlane.Normal, LightDelta);
}

void AGeometricTestsGameModeBase::TEST_SolveBarycentricCoordinates3D()
{
	const Vector3<float> Vertices[3] = { Vector3<float>(-3.f, -5.f, 7.f)
									   , Vector3<float>(1.f, 2.f, 1.f)
									   , Vector3<float>(3, 2.f, 3.f) };
	
	const Vector3<float> Point = Vector3<float>(0.f, 6.f, 4.f);
	
	float Barycentric[3] = { 1, 1, 1 };
		
	UGeometricTestLibrary::SolveBarycentricCoordinates3D(Vertices, Point, Barycentric);
}

