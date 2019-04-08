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
#include "GameFramework/GameModeBase.h"
#include "GeometricTestsGameModeBase.generated.h"
 
/**
 * 
 */
UCLASS()
class GEOMETRICTESTS_API AGeometricTestsGameModeBase : public AGameModeBase
{
	GENERATED_BODY()

public:
	
	AGeometricTestsGameModeBase() = default;
	~AGeometricTestsGameModeBase() = default;
	
	// UE4: UObject copying and moving disabled
	
protected:

	virtual void BeginPlay() override;
	
private:
	
	/*
	 * Closest Point Tests
	*/
	void TEST_ClosestPointInAABB();
	void TEST_ClosestPointOnRay();
	void TEST_ClosestPointOnPlane();
	void TEST_ClosestPointOnSphere();
	
	/*
	 * Intersection Tests
	*/
	void TEST_DoesAABBPlaneIntersect();
	void TEST_DoLinesIntersect();
	void TEST_DoesRayAABBIntersect();
	void TEST_DoesRayPlaneIntersect();
	void TEST_DoesRaySphereIntersect();
	void TEST_DoesRayTriangleIntersect();
	void TEST_DoRaysIntersect();
	void TEST_DoesSpherePlaneIntersect_Dynamic();
	void TEST_DoesSpherePlaneIntersect_Static();
	void TEST_DoSpheresIntersect_Dynamic();
	void TEST_DoSpheresIntersect_Static();
	void TEST_DoThreePlanesIntersect();
	void TEST_DoAABBsIntersect();
		
	/*
	 * Reflection Vector
	*/
	void TEST_SolveReflectionVector();
	
	/*
	 * Barycentric Coordinates
	*/
	void TEST_SolveBarycentricCoordinates3D();
};
