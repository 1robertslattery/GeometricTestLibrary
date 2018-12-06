# Geometric Test Library
Solves 2D and 3D math problems, including, `closest point`, `intersection`, `line of sight` and `reflection vector`.

## OVERVIEW

This project aims to teach 3D math. Everything can be found in the header file `GeometricTestLibrary.h`.  
Since these tests are commonly used in video games, this repository includes the Unreal Engine 4 project `GeometricTests.uproject`. 
These tests are designed to be used with any `C++` application or game engine, however. The project is released under the Apache 2.0 License.

**Examples:**
- 1. Closest Point in AABB
- 2. Closest Point on a Ray
- 3. Closest Point on a Plane
- 4. Closest Point on a Sphere
- 5. Intersection AABB-Plane
- 6. Intersection of Two Lines in 2D
- 7. Intersection Ray-AABB
- 8. Intersection Ray-Plane
- 9. Intersection Ray-Triangle
- 10. Intersection of Two Rays in 3D
- 11. Dynamic Intersection Sphere-Plane
- 12. Static Intersection Sphere-Plane
- 13. Dynamic Intersection of Two Spheres
- 14. Static Intersection of Two Spheres
- 15. Intersection of Three Planes
- 16. Intersection of Two AABBs
- 17. Reflection Vector
- 18. Barycentric Coordinates of Triangle in 3D

## FILES AND FOLDERS

| Files and Folders	| Description						|
| --------------------- |:-----------------------------------------------------:|
| `Config`		| UE4 project configuration files	|
| `Content`		| UE4 project files	|
| `Source`			| Project source						|
| `CHANGELOG`		| Log to track changes in respository			|
| `LICENSE`		| Apache 2.0 License			|
| `GeometricTests.uproject`		| UE4 project			|
| `README`		| This file						|

### Config

| Files		| Description						|
| ----------------------------- |:-----------------------------------------------------:|
| `DefaultEditor.ini`		| UE4 Config File				|
| `DefaultEngine.ini`			| UE4 Config File				|
| `DefaultGame.ini`			| UE4 Config File					|

### Content/Maps

| Files		| Description						|
| ----------------------------- |:-----------------------------------------------------:|
| `EmptyLevel.umap`		| UE4 Level				|

### Source

| Files and Folders		| Description						|
| ----------------------------- |:-----------------------------------------------------:|
| `GeometricTests`		| Geometric Tests Source				|
| `GeometricTests.Target.cs`			| UE4 Build System				|
| `GeometricTestsEditor.Target.cs`			| UE4 Build System					|

### Source/GeometricTests

| Files and Folders				| Description						|
| ----------------------------- |:-----------------------------------------------------:|
| `Math`		| Math Directory			|
| `Primitives`		| Geometric Primitives Directory			|
| `GeometricTestLibrary.h`			| Geometric Test Library						|
| `GeometricTest.Build.cs`			| UE4 Build System							|
| `GeometricTests.h`			| Header: Module implementation for UE4 Project			|
| `GeometricTests.cpp`			| Source: Module implementation for UE4 Project			|
| `GeometricTestsGameModeBase.h`			| Header: UE4 Game Mode Base		|
| `GeometricTestsGameModeBase.cpp`			| Source: UE4 Game Mode Base				|

### Source/GeometricTests/Math/LinearAlgebra

| Files		| Description						|
| ------------- |:-----------------------------------------------------:|
| `Vector.h`	| Templated structs that handle vector math						|

### Source/GeometricTests/Math/Operations

| Files		| Description						|
| ------------- |:-----------------------------------------------------:|
| `MathOperations.h`	| Common math functions					|

### Source/GeometricTests/Primitives

| Files		| Description						|
| ------------- |:-----------------------------------------------------:|
| `AABB.h`	| Axis-Aligned Bounding Box Data				|
| `Circle.h`	| Circle Data				|
| `Plane.h`	| Plane Data				|
| `Ray.h`	| Ray Data				|

## LICENSE

> Apache 2.0 License
>
> Copyright (c) 2018 Robert Slattery
>
> Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
>
> http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the License for the specific language governing permissions and limitations under the License.
