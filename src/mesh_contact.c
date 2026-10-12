// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

#include "contact.h"
#include "manifold.h"
#include "physics_world.h"
#include "qsort.h"
#include "shape.h"
#include "voxel_grid.h"

#include "box3d/types.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

// This guards against excessive memory usage and complex collision
#define B3_MAX_MESH_CONTACT_TRIANGLES 256
#define B3_MAX_POINTS_PER_TRIANGLE 32

#if B3_ENABLE_VALIDATION
static bool b3IsSorted( const int* array, int count )
{
	for ( int i = 0; i < count - 1; ++i )
	{
		if ( array[i] >= array[i + 1] )
		{
			return false;
		}
	}

	return true;
}
#endif

typedef struct b3TriangleQueryContext
{
	int* indices;
	int capacity;
	int count;
} b3TriangleQueryContext;

static bool b3CollectTriangleIndicesCallback( b3Vec3 a, b3Vec3 b, b3Vec3 c, int triangleIndex, void* context )
{
	B3_UNUSED( a, b, c );
	b3TriangleQueryContext* triangleContext = (b3TriangleQueryContext*)context;
	if ( triangleContext->count == triangleContext->capacity )
	{
		return false;
	}
	triangleContext->indices[triangleContext->count] = triangleIndex;
	triangleContext->count += 1;
	return triangleContext->count < triangleContext->capacity;
}

static int b3QueryMeshTriangles( int* indices, int capacity, const b3Mesh* mesh, b3AABB bounds )
{
	b3TriangleQueryContext context = {
		.indices = indices,
		.capacity = capacity,
		.count = 0,
	};

	b3QueryMesh( mesh, bounds, b3CollectTriangleIndicesCallback, &context );
	return context.count;
}

static int b3QueryHeightFieldTriangles( int* indices, int capacity, const b3HeightFieldData* heightField, b3AABB bounds )
{
	b3TriangleQueryContext context = {
		.indices = indices,
		.capacity = capacity,
		.count = 0,
	};

	b3QueryHeightField( heightField, bounds, b3CollectTriangleIndicesCallback, &context );
	return context.count;
}

static void b3RefreshCache( b3Contact* contact, const b3Shape* shapeA, b3WorldTransform xfA, const b3AABB* bounds )
{
	B3_ASSERT( shapeA->type == b3_meshShape || shapeA->type == b3_heightShape );

	b3MeshContact* meshContact = &contact->meshContact;

	// If the dynamic body didn't move out of the cached query bounds we are done!
	if ( b3AABB_Contains( meshContact->queryBounds, *bounds ) )
	{
		if ( shapeA->type == b3_meshShape )
		{
			for ( int i = 0; i < contact->meshContact.triangleCache.count; ++i )
			{
				B3_ASSERT( 0 <= contact->meshContact.triangleCache.data[i].triangleIndex &&
						   contact->meshContact.triangleCache.data[i].triangleIndex < shapeA->mesh.data->triangleCount );
			}
		}

		return;
	}

	// Enlarge to the query bounds to absorb small movement
	float radius = B3_MAX_AABB_MARGIN + B3_SPECULATIVE_DISTANCE;
	b3Vec3 extension = { radius, radius, radius };
	meshContact->queryBounds.lowerBound = b3Sub( bounds->lowerBound, extension );
	meshContact->queryBounds.upperBound = b3Add( bounds->upperBound, extension );

	// Query triangles
	int triangleCapacity = B3_MAX_MESH_CONTACT_TRIANGLES;

	int triangleIndices[B3_MAX_MESH_CONTACT_TRIANGLES];

	// Bounds are in world space. Convert to the local mesh frame. The broadphase bounds are float,
	// so the demoted mesh transform is the matching float world frame (exact in float mode).
	b3Transform meshTransform = b3ToRelativeTransform( xfA, b3Pos_zero );
	b3AABB localBounds = b3AABB_Transform( b3InvertTransform( meshTransform ), meshContact->queryBounds );
	int triangleCount;
	if ( shapeA->type == b3_meshShape )
	{
		triangleCount = b3QueryMeshTriangles( triangleIndices, triangleCapacity, &shapeA->mesh, localBounds );
	}
	else
	{
		B3_ASSERT( shapeA->type == b3_heightShape );
		triangleCount = b3QueryHeightFieldTriangles( triangleIndices, triangleCapacity, shapeA->heightField, localBounds );
	}

	if ( triangleCount == triangleCapacity )
	{
		static bool s_once = false;
		if ( s_once == false )
		{
			b3Log( "WARNING: complex mesh detected, triangle buffer capacity of %d reached", triangleCapacity );
			s_once = true;
		}
	}

	// Triangle indices must be sorted to match caches.
	B3_VALIDATE( b3IsSorted( triangleIndices, triangleCount ) );

	// Create new contact cache and match with old one
	b3ContactCache contactCache[B3_MAX_MESH_CONTACT_TRIANGLES];

	int index2 = 0;
	for ( int index1 = 0; index1 < triangleCount; ++index1 )
	{
		contactCache[index1] = (b3ContactCache){ 0 };

		while ( index2 < contact->meshContact.triangleCache.count &&
				contact->meshContact.triangleCache.data[index2].triangleIndex < triangleIndices[index1] )
		{
			index2 += 1;
		}

		if ( index2 < contact->meshContact.triangleCache.count &&
			 contact->meshContact.triangleCache.data[index2].triangleIndex == triangleIndices[index1] )
		{
			contactCache[index1] = contact->meshContact.triangleCache.data[index2].cache;
		}
	}

	// Save new cache
	b3Array_Resize( contact->meshContact.triangleCache, triangleCount );
	for ( int i = 0; i < triangleCount; ++i )
	{
		contact->meshContact.triangleCache.data[i] = (b3TriangleCache){ triangleIndices[i], contactCache[i] };

		if ( shapeA->type == b3_meshShape )
		{
			B3_ASSERT( 0 <= contact->meshContact.triangleCache.data[i].triangleIndex &&
					   contact->meshContact.triangleCache.data[i].triangleIndex < shapeA->mesh.data->triangleCount );
		}
	}
}

typedef struct b3TentativeTriangle
{
	float squaredDistance;
	int index;
} b3TentativeTriangle;

#define B3_MAX_EDGE_COUNT 64

typedef struct b3FoundEdges
{
	uint64_t keys[B3_MAX_EDGE_COUNT];
	int count;
} b3FoundEdges;

static inline bool b3AddEdge( b3FoundEdges* edges, int vertex1, int vertex2 )
{
	uint64_t i1 = (uint64_t)b3MinInt( vertex1, vertex2 );
	uint64_t i2 = (uint64_t)b3MaxInt( vertex1, vertex2 );
	uint64_t key = i1 << 32 | i2;

	int count = edges->count;
	for ( int i = 0; i < count; ++i )
	{
		if ( edges->keys[i] == key )
		{
			return false;
		}
	}

	if ( count == B3_MAX_EDGE_COUNT )
	{
		// This will lead to a potential ghost collision
		return true;
	}

	edges->keys[count] = key;
	edges->count += 1;

	return true;
}

static inline bool b3FindEdge( b3FoundEdges* edges, int vertex1, int vertex2 )
{
	uint64_t i1 = (uint64_t)b3MinInt( vertex1, vertex2 );
	uint64_t i2 = (uint64_t)b3MaxInt( vertex1, vertex2 );
	uint64_t key = i1 << 32 | i2;

	int count = edges->count;
	for ( int i = 0; i < count; ++i )
	{
		if ( edges->keys[i] == key )
		{
			return true;
		}
	}

	return false;
}

#if 0
// Two triangles share an edge iff they share at least two vertex indices.
static inline bool b3TrianglesShareEdge( int a1, int a2, int a3, int b1, int b2, int b3 )
{
	int matches = 0;
	matches += ( a1 == b1 || a1 == b2 || a1 == b3 );
	matches += ( a2 == b1 || a2 == b2 || a2 == b3 );
	matches += ( a3 == b1 || a3 == b2 || a3 == b3 );
	return matches >= 2;
}
#endif

#define B3_MAX_VERTEX_COUNT 64

typedef struct b3FoundVertices
{
	int keys[B3_MAX_VERTEX_COUNT];
	int count;
} b3FoundVertices;

static inline bool b3AddVertex( b3FoundVertices* vertices, int vertex )
{
	int key = vertex;

	int count = vertices->count;
	for ( int i = 0; i < count; ++i )
	{
		if ( vertices->keys[i] == key )
		{
			return false;
		}
	}

	if ( count == B3_MAX_VERTEX_COUNT )
	{
		// This will lead to a potential ghost collision
		return true;
	}

	vertices->keys[count] = key;
	vertices->count += 1;

	return true;
}

// Returns true if (score, separation) should replace (bestScore, bestSeparation).
static inline bool b3IsBetterCullCandidate( float score, float separation, float bestScore, float bestSeparation, float scoreTol,
											float separationTol )
{
	if ( score > bestScore + scoreTol )
	{
		return true;
	}
	if ( score < bestScore - scoreTol )
	{
		return false;
	}

	// Break the tie using separation
	return separation < bestSeparation - separationTol;
}

typedef struct b3Point2D
{
	b3Vec2 p;
	float separation;
	int originalIndex;
} b3Point2D;

static int b3CullPoints( b3Point2D* points, int count )
{
	if ( count <= 1 )
	{
		return count;
	}

	float tol = 0.25f * B3_LINEAR_SLOP;
	float tolSqr = tol * tol;
	float separationTol = B3_LINEAR_SLOP;

	b3Point2D finalPoints[4];
	int count1 = count;

	// Step 1: the two points with the largest distance, ties broken by deepest combined separation
	float bestScore = 0.0f;
	float bestSeparation = FLT_MAX;
	int bestIndex1 = B3_NULL_INDEX;
	int bestIndex2 = B3_NULL_INDEX;

	for ( int i = 0; i < count1; ++i )
	{
		b3Vec2 p1 = points[i].p;
		for ( int j = i + 1; j < count1; ++j )
		{
			float score = b3DistanceSquared2( p1, points[j].p );
			// Separation sum heuristic
			float separation = points[i].separation + points[j].separation;

			if ( b3IsBetterCullCandidate( score, separation, bestScore, bestSeparation, tolSqr, separationTol ) )
			{
				bestIndex1 = i;
				bestIndex2 = j;
				bestScore = score;
				bestSeparation = separation;
			}
		}
	}

	if ( bestScore < tolSqr )
	{
		// Choose deepest point
		int deepestIndex = 0;
		for ( int i = 1; i < count1; ++i )
		{
			if ( points[i].separation < points[deepestIndex].separation )
			{
				deepestIndex = i;
			}
		}

		if ( deepestIndex != 0 )
		{
			points[0] = points[deepestIndex];
		}
		return 1;
	}

	finalPoints[0] = points[bestIndex1];
	finalPoints[1] = points[bestIndex2];

	// Cull
	points[bestIndex2] = points[count1 - 1];
	points[bestIndex1] = points[count1 - 2];
	count1 -= 2;

	if ( count1 == 0 )
	{
		points[0] = finalPoints[0];
		points[1] = finalPoints[1];
		return 2;
	}

	// First anchor point
	b3Vec2 a = finalPoints[0].p;

	// Second anchor point
	b3Vec2 b = finalPoints[1].p;
	b3Vec2 ba = b3Sub2( b, a );
	// float length = b3Length2( ba );
	// float areaTol = tol * length;

	// Step 2: find the point with the maximum triangular area, ties broken by deepest separation
	bestScore = 0.0f;
	bestSeparation = FLT_MAX;
	int bestIndex = B3_NULL_INDEX;
	float bestSignedArea = 0.0f;
	for ( int i = 0; i < count1; ++i )
	{
		b3Vec2 p = points[i].p;
		float signedArea = b3Cross2( ba, b3Sub2( p, a ) );
		float score = b3AbsFloat( signedArea );

		if ( b3IsBetterCullCandidate( score, points[i].separation, bestScore, bestSeparation, tolSqr, separationTol ) )
		{
			bestSignedArea = signedArea;
			bestScore = score;
			bestSeparation = points[i].separation;
			bestIndex = i;
		}
	}

	if ( bestIndex == B3_NULL_INDEX )
	{
		// All points collinear
		points[0] = finalPoints[0];
		points[1] = finalPoints[1];
		return 2;
	}

	// Store best point
	finalPoints[2] = points[bestIndex];

	if ( count1 == 1 )
	{
		points[0] = finalPoints[0];
		points[1] = finalPoints[1];
		points[2] = finalPoints[2];
		return 3;
	}

	// Cull
	points[bestIndex] = points[count1 - 1];
	count1 -= 1;

	// Step 4: get the point that adds the most area outside the current triangle

	// Third anchor
	b3Vec2 c = finalPoints[2].p;

	// Ensure CCW ordering
	if ( bestSignedArea < 0.0f )
	{
		B3_SWAP( b, c );
		ba = b3Sub2( b, a );
	}

	b3Vec2 cb = b3Sub2( c, b );
	b3Vec2 ac = b3Sub2( a, c );

	bestScore = 0.0f;
	bestSeparation = FLT_MAX;
	bestIndex = B3_NULL_INDEX;
	for ( int i = 0; i < count1; ++i )
	{
		b3Vec2 p = points[i].p;
		float u1 = b3Cross2( b3Sub2( p, a ), ba );
		float u2 = b3Cross2( b3Sub2( p, b ), cb );
		float u3 = b3Cross2( b3Sub2( p, c ), ac );
		float score = b3MaxFloat( u1, b3MaxFloat( u2, u3 ) );

		// Use the area tolerance for collinear points and hysteresis
		if ( b3IsBetterCullCandidate( score, points[i].separation, bestScore, bestSeparation, tolSqr, separationTol ) )
		{
			bestScore = score;
			bestSeparation = points[i].separation;
			bestIndex = i;
		}
	}

	if ( bestIndex == B3_NULL_INDEX )
	{
		// No additional area
		points[0] = finalPoints[0];
		points[1] = finalPoints[1];
		points[2] = finalPoints[2];
		return 3;
	}

	// Store best point
	finalPoints[3] = points[bestIndex];

	// Full quad
	points[0] = finalPoints[0];
	points[1] = finalPoints[1];
	points[2] = finalPoints[2];
	points[3] = finalPoints[3];
	return 4;
}

// A cluster of many points, such as the box pairs of two flat voxel grid surfaces, costs the exact reduction a number of
// distance evaluations that grows with the square of the count. Above this count the points are first cut to the extremes
// along eight directions of the cluster plane and the deepest point, which keeps the outline and the depth.
#define B3_CLUSTER_PREFILTER_COUNT 48

static int b3PrefilterClusterPoints( b3LocalManifoldPoint* points, int count, b3Vec3 u, b3Vec3 v )
{
	// Four axes and their opposites: x, y, x + y, x - y
	int highest[4], lowest[4];
	float highValue[4], lowValue[4];
	for ( int k = 0; k < 4; ++k )
	{
		highest[k] = lowest[k] = 0;
		highValue[k] = -FLT_MAX;
		lowValue[k] = FLT_MAX;
	}

	int deepest = 0;
	b3Vec3 origin = points[0].point;
	for ( int i = 0; i < count; ++i )
	{
		b3Vec3 d = b3Sub( points[i].point, origin );
		float x = b3Dot( d, u ), y = b3Dot( d, v );
		float values[4] = { x, y, x + y, x - y };
		for ( int k = 0; k < 4; ++k )
		{
			if ( values[k] > highValue[k] )
			{
				highValue[k] = values[k];
				highest[k] = i;
			}

			if ( values[k] < lowValue[k] )
			{
				lowValue[k] = values[k];
				lowest[k] = i;
			}
		}

		if ( points[i].separation < points[deepest].separation )
		{
			deepest = i;
		}
	}

	int candidates[9];
	int candidateCount = 0;
	int all[9] = { highest[0], lowest[0], highest[1], lowest[1], highest[2], lowest[2], highest[3], lowest[3], deepest };
	for ( int n = 0; n < 9; ++n )
	{
		bool found = false;
		for ( int m = 0; m < candidateCount; ++m )
		{
			found = found || candidates[m] == all[n];
		}

		if ( found == false )
		{
			candidates[candidateCount++] = all[n];
		}
	}

	b3LocalManifoldPoint kept[9];
	for ( int n = 0; n < candidateCount; ++n )
	{
		kept[n] = points[candidates[n]];
	}

	memcpy( points, kept, candidateCount * sizeof( b3LocalManifoldPoint ) );
	return candidateCount;
}

static int b3ReduceCluster( b3LocalManifoldPoint* points, int count1, b3Vec3 normal, b3Arena arena )
{
	int targetCount = 1;
	if ( count1 <= targetCount )
	{
		return count1;
	}

	b3Vec3 u = b3Perp( normal );
	b3Vec3 v = b3Cross( normal, u );
	if ( count1 > B3_CLUSTER_PREFILTER_COUNT )
	{
		count1 = b3PrefilterClusterPoints( points, count1, u, v );
	}

	b3Point2D* pts = b3Bump( &arena, count1 * sizeof( b3Point2D ) );
	b3Vec3 origin = points[0].point;

	for ( int i = 0; i < count1; ++i )
	{
		b3Vec3 d = b3Sub( points[i].point, origin );
		pts[i].p = (b3Vec2){ b3Dot( d, u ), b3Dot( d, v ) };
		pts[i].separation = points[i].separation;
		pts[i].originalIndex = i;
	}

	int count2 = b3CullPoints( pts, count1 );
	B3_ASSERT( count2 <= B3_MAX_MANIFOLD_POINTS );

	b3LocalManifoldPoint finalPoints[B3_MAX_MANIFOLD_POINTS];
	for ( int i = 0; i < count2; ++i )
	{
		int index = pts[i].originalIndex;
		B3_ASSERT( 0 <= index && index < count1 );
		finalPoints[i] = points[index];
	}

	memcpy( points, finalPoints, count2 * sizeof( b3LocalManifoldPoint ) );
	return count2;
}

typedef struct b3Cluster
{
	b3Vec3 manifoldNormal;
	b3Vec3 triangleNormal;
	b3LocalManifoldPoint* points;
	int pointCapacity;
	int pointCount;
} b3Cluster;

// Clusters accepted manifolds by normal, reduces each cluster to one manifold and matches them with the previous manifolds
// so impulses carry over. The local manifolds are in frame B. Shared by meshes, height fields and voxel grids. Returns false when
// there is nothing to solve, with the contact manifolds freed.
static bool b3BuildClusterManifolds( b3World* world, b3Contact* contact, b3LocalManifold** acceptedManifolds,
									 int acceptedManifoldCount, b3WorldTransform xfA, b3WorldTransform xfB, float restOffset,
									 b3Arena arena )
{
	if ( acceptedManifoldCount == 0 )
	{
		if ( contact->manifoldCount > 0 )
		{
			b3FreeManifolds( world, contact->manifolds, contact->manifoldCount );
			contact->manifolds = NULL;
			contact->manifoldCount = 0;
		}
		return false;
	}

	b3Cluster* clusters = b3Bump( &arena, acceptedManifoldCount * sizeof( b3Cluster ) );
	int* clusterMemberships = b3Bump( &arena, acceptedManifoldCount * sizeof( int ) );

	// Cluster tolerance is tighter than the warm starting manifold matching tolerance. These
	// serve different purposes.
	const float clusterThreshold = 0.996f;
	int clusterCount = 0;
	int clusterPointCount = 0;
	for ( int i = 0; i < acceptedManifoldCount; ++i )
	{
		clusterMemberships[i] = B3_NULL_INDEX;

		const b3LocalManifold* manifold = acceptedManifolds[i];
		clusterPointCount += manifold->pointCount;

		// Cluster based on the triangle normal and contact normal.
		// The first cluster found is accepted because the tolerance is tight.
		// todo consider requiring the triangles to be connect by an edge.
		// todo consider looking for the best cluster instead of the first one within tolerance
		// This bool is here to allow quick testing with and without clustering.
		bool allowClustering = true;
		b3Vec3 manifoldNormal = manifold->normal;
		b3Vec3 triangleNormal = manifold->triangleNormal;
		int clusterIndex = B3_NULL_INDEX;
		for ( int j = 0; j < clusterCount && allowClustering; ++j )
		{
			float cosManifoldAngle = b3Dot( clusters[j].manifoldNormal, manifoldNormal );
			float cosTriangleAngle = b3Dot( clusters[j].triangleNormal, triangleNormal );
			if ( cosManifoldAngle <= clusterThreshold || cosTriangleAngle <= clusterThreshold )
			{
				continue;
			}

#if 0
			// todo there could be later triangles that create the connection
			// then failure to cluster breaks greedy impulse warm starting
			bool edgeConnected = false;

			for ( int k = 0; k < i; ++k )
			{
				if ( clusterMemberships[k] != j )
				{
					continue;
				}

				const b3LocalManifold* other = acceptedManifolds[k];
				if ( b3TrianglesShareEdge( manifold->i1, manifold->i2, manifold->i3, other->i1, other->i2, other->i3 ) )
				{
					edgeConnected = true;
					break;
				}
			}

			if ( edgeConnected )
			{
				clusterIndex = j;
				break;
			}
#else

			// Found a cluster
			clusterIndex = j;
			break;
#endif
		}

		if ( clusterIndex != B3_NULL_INDEX )
		{
			clusterMemberships[i] = clusterIndex;
			clusters[clusterIndex].pointCapacity += manifold->pointCount;
		}
		else
		{
			clusters[clusterCount].manifoldNormal = manifoldNormal;
			clusters[clusterCount].triangleNormal = triangleNormal;
			clusters[clusterCount].pointCapacity = manifold->pointCount;
			clusterMemberships[i] = clusterCount;
			clusterCount += 1;
		}
	}

	if ( clusterPointCount == 0 )
	{
		return false;
	}

	// Setup clusters
	b3LocalManifoldPoint* clusterPoints = b3Bump( &arena, clusterPointCount * sizeof( b3LocalManifoldPoint ) );
	int pointOffset = 0;

	for ( int i = 0; i < clusterCount; ++i )
	{
		b3Cluster* cluster = clusters + i;
		cluster->points = clusterPoints + pointOffset;
		cluster->pointCount = 0;
		pointOffset += cluster->pointCapacity;
	}

	// Populate clusters
	for ( int i = 0; i < acceptedManifoldCount; ++i )
	{
		int clusterIndex = clusterMemberships[i];
		if ( clusterIndex == B3_NULL_INDEX )
		{
			continue;
		}

		B3_ASSERT( 0 <= clusterIndex && clusterIndex < clusterCount );

		b3LocalManifold* am = acceptedManifolds[i];
		b3Cluster* cm = clusters + clusterIndex;
		for ( int j = 0; j < am->pointCount; ++j )
		{
			B3_ASSERT( cm->pointCount < cm->pointCapacity );
			b3LocalManifoldPoint* ap = am->points + j;
			b3LocalManifoldPoint* cp = cm->points + cm->pointCount;

			cp->triangleIndex = am->triangleIndex;
			cp->point = ap->point;
			cp->separation = ap->separation;
			cp->pair = ap->pair;
			cm->pointCount += 1;
		}
	}

	// Simplify clusters
	for ( int i = 0; i < clusterCount; ++i )
	{
		b3Cluster* cm = clusters + i;
		B3_ASSERT( cm->pointCount == cm->pointCapacity );
		int reducedCount = b3ReduceCluster( cm->points, cm->pointCount, cm->triangleNormal, arena );
		cm->pointCount = reducedCount;
	}

	// Make a temporary copy of previous manifolds
	int oldManifoldCount = contact->manifoldCount;
	b3Manifold* oldManifolds = NULL;
	if ( oldManifoldCount > 0 )
	{
		oldManifolds = b3Bump( &arena, oldManifoldCount * sizeof( b3Manifold ) );
		memcpy( oldManifolds, contact->manifolds, oldManifoldCount * sizeof( b3Manifold ) );
	}

	// Resize manifolds if needed
	if ( oldManifoldCount != clusterCount )
	{
		b3FreeManifolds( world, contact->manifolds, contact->manifoldCount );
		contact->manifolds = b3AllocateManifolds( world, clusterCount );
		contact->manifoldCount = (uint16_t)clusterCount;
	}
	else
	{
		// Mem zero manifolds
		memset( contact->manifolds, 0, contact->manifoldCount * sizeof( b3Manifold ) );
	}

	bool* consumed = NULL;
	if ( oldManifoldCount > 0 )
	{
		consumed = b3Bump( &arena, oldManifoldCount * sizeof( bool ) );
		memset( consumed, 0, oldManifoldCount * sizeof( bool ) );
	}

	b3Matrix3 matrixB = b3MakeMatrixFromQuat( xfB.q );
	b3Vec3 offsetA = b3SubPos( xfB.p, xfA.p );

	const float normalMatchTolerance = 0.995f;
	for ( int i = 0; i < clusterCount; ++i )
	{
		b3Cluster* cm = clusters + i;
		int pointCount = cm->pointCount;
		B3_ASSERT( 0 < pointCount && pointCount <= B3_MAX_MANIFOLD_POINTS );

		b3Manifold* manifold = contact->manifolds + i;
		manifold->pointCount = pointCount;
		manifold->normal = b3MulMV( matrixB, cm->manifoldNormal );

		b3Vec3 clusterNormal = b3MulMV( matrixB, cm->manifoldNormal );
		float bestDot = normalMatchTolerance;
		int bestIndex = B3_NULL_INDEX;

		for ( int j = 0; j < oldManifoldCount; ++j )
		{
			if ( consumed[j] == true )
			{
				continue;
			}

			float dot = b3Dot( oldManifolds[j].normal, clusterNormal );
			if ( dot > bestDot )
			{
				bestIndex = j;
				bestDot = dot;
			}
		}

		b3Manifold* matchedManifold = NULL;
		if ( bestIndex != B3_NULL_INDEX )
		{
			matchedManifold = oldManifolds + bestIndex;
			manifold->frictionImpulse = matchedManifold->frictionImpulse;
			manifold->rollingImpulse = matchedManifold->rollingImpulse;
			manifold->twistImpulse = matchedManifold->twistImpulse;
			consumed[bestIndex] = true;
		}

		for ( int j = 0; j < pointCount; ++j )
		{
			const b3LocalManifoldPoint* source = cm->points + j;
			b3ManifoldPoint* target = manifold->points + j;

			// Contact points are computed in frame B
			target->anchorB = b3MulMV( matrixB, source->point );
			target->anchorA = b3Add( target->anchorB, offsetA );
			target->separation = source->separation - restOffset;
			target->featureId = b3MakeFeatureId( source->pair );
			target->triangleIndex = source->triangleIndex;

			// Preserve normal impulse if possible
			if ( matchedManifold != NULL )
			{
				int oldPointCount = matchedManifold->pointCount;
				for ( int k = 0; k < oldPointCount; ++k )
				{
					b3ManifoldPoint* oldPt = matchedManifold->points + k;

					if ( target->featureId == oldPt->featureId && target->triangleIndex == oldPt->triangleIndex )
					{
						target->normalImpulse = oldPt->normalImpulse;
						target->persisted = true;

						// claimed
						oldPt->triangleIndex = B3_NULL_INDEX;
						break;
					}
				}
			}
		}
	}

	return true;
}

// A mesh or height field A and one convex shape B, collided triangle by triangle. A voxel grid contact uses it once for every
// box of the grid, each box being a hull in the frame of its cell.
typedef struct b3MeshConvexInput
{
	b3TaskContext* context;
	const b3Shape* meshShape;
	b3ShapeType typeB;
	const b3HullData* hullB;
	const b3Capsule* capsuleB;
	const b3Sphere* sphereB;

	// From the mesh frame into the frame of B
	b3Transform transformAtoB;
	bool isFast;
	bool enableSpeculative;

	// The faces of hull B that solid covers, for a box of a voxel grid. A contact that would push through one is answered on
	// the exposed faces instead.
	uint8_t coveredB;
} b3MeshConvexInput;

// Finds the manifolds of the triangles in the caches against B. The accepted manifolds point into manifoldBuffer and
// pointBuffer, the caller's memory, because they outlive this call. Returns their count.
static int b3CollideMeshTriangles( const b3MeshConvexInput* input, b3TriangleCache* triangleCaches, int triangleCount,
								   b3LocalManifold** acceptedManifolds, b3LocalManifold* manifoldBuffer,
								   b3LocalManifoldPoint* pointBuffer, int pointBufferCapacity, b3Arena arena )
{
	const b3Shape* shapeA = input->meshShape;
	b3TaskContext* context = input->context;
	bool isFast = input->isFast;
	bool enableSpeculative = input->enableSpeculative;
	uint8_t coveredB = input->coveredB;

	int acceptedManifoldCount = 0;
	b3LocalManifold** tentativeManifolds = b3Bump( &arena, triangleCount * sizeof( b3LocalManifold* ) );
	int tentativeManifoldCount = 0;
	b3TentativeTriangle* tentativeTriangles = b3Bump( &arena, triangleCount * sizeof( b3TentativeTriangle ) );
	int tentativeTriangleCount = 0;

	b3FoundEdges foundEdges;
	b3FoundVertices foundVertices;
	foundEdges.count = 0;
	foundVertices.count = 0;

	// This transform converts from mesh frame into the shapeB frame
	b3Transform transformAtoB = input->transformAtoB;
	b3Matrix3 relativeMatrix = b3MakeMatrixFromQuat( transformAtoB.q );
	float linearSlop = B3_LINEAR_SLOP;

	int totalPointCount = 0;
	int manifoldCount = 0;

	const b3HullData* hullB = input->hullB;

	for ( int index = 0; index < triangleCount && totalPointCount + 3 < pointBufferCapacity; ++index )
	{
		int triangleIndex = triangleCaches[index].triangleIndex;

		b3Triangle triangle;
		if ( shapeA->type == b3_meshShape )
		{
			triangle = b3GetMeshTriangle( &shapeA->mesh, triangleIndex );
		}
		else
		{
			B3_ASSERT( shapeA->type == b3_heightShape );
			triangle = b3GetHeightFieldTriangle( shapeA->heightField, triangleIndex );
		}

		// Transform triangle into the shape frame
		b3Vec3 vertices[3];
		vertices[0] = b3Add( b3MulMV( relativeMatrix, triangle.vertices[0] ), transformAtoB.p );
		vertices[1] = b3Add( b3MulMV( relativeMatrix, triangle.vertices[1] ), transformAtoB.p );
		vertices[2] = b3Add( b3MulMV( relativeMatrix, triangle.vertices[2] ), transformAtoB.p );

		if ( shapeA->type == b3_meshShape && shapeA->mesh.data->doubleSided )
		{
			// A double sided mesh collides from the side the other shape is on: wind the triangle to face it
			b3Vec3 centerB;
			switch ( input->typeB )
			{
				case b3_capsuleShape:
					centerB = b3Lerp( input->capsuleB->center1, input->capsuleB->center2, 0.5f );
					break;

				case b3_hullShape:
					centerB = input->hullB->center;
					break;

				default:
					centerB = input->sphereB->center;
					break;
			}

			b3Vec3 faceNormal = b3Cross( b3Sub( vertices[1], vertices[0] ), b3Sub( vertices[2], vertices[0] ) );
			if ( b3Dot( faceNormal, b3Sub( centerB, vertices[0] ) ) < 0.0f )
			{
				b3Vec3 swappedVertex = vertices[1];
				vertices[1] = vertices[2];
				vertices[2] = swappedVertex;

				int swappedIndex = triangle.i2;
				triangle.i2 = triangle.i3;
				triangle.i3 = swappedIndex;

				// Edges 1 and 3 trade places, and what was concave from the front is convex from the back
				int flags = triangle.flags;
				int swappedFlags = 0;
				swappedFlags |= ( flags & b3_inverseConcaveEdge3 ) ? b3_concaveEdge1 : 0;
				swappedFlags |= ( flags & b3_inverseConcaveEdge2 ) ? b3_concaveEdge2 : 0;
				swappedFlags |= ( flags & b3_inverseConcaveEdge1 ) ? b3_concaveEdge3 : 0;
				swappedFlags |= ( flags & b3_concaveEdge3 ) ? b3_inverseConcaveEdge1 : 0;
				swappedFlags |= ( flags & b3_concaveEdge2 ) ? b3_inverseConcaveEdge2 : 0;
				swappedFlags |= ( flags & b3_concaveEdge1 ) ? b3_inverseConcaveEdge3 : 0;
				triangle.flags = swappedFlags;
			}
		}

		b3ContactCache* cache = &triangleCaches[index].cache;
		int pointCapacity = pointBufferCapacity - totalPointCount;
		b3LocalManifold* manifold = manifoldBuffer + manifoldCount;
		manifold->points = pointBuffer + totalPointCount;
		manifold->pointCount = 0;
		manifold->triangleFlags = triangle.flags;
		manifold->feature = b3_featureNone;

		switch ( input->typeB )
		{
			case b3_capsuleShape:
				b3CollideTriangleAndCapsule( manifold, pointCapacity, vertices, input->capsuleB, &cache->simplexCache );
				break;

			case b3_hullShape:
				// Cached edge contact is dangerous at high speed because the hull can rotate around the edge and tunnel
				// through the triangle.
				if ( isFast && cache->satCache.type == b3_edgePairAxis )
				{
					cache->satCache = (b3SATCache){ 0 };
				}

				b3CollideTriangleAndHull( manifold, pointCapacity, vertices[0], vertices[1], vertices[2], triangle.flags, hullB,
										  &cache->satCache, enableSpeculative );
				context->satCallCount += 1;
				context->satCacheHitCount += cache->satCache.hit;

				if ( coveredB != 0 && manifold->pointCount > 0 &&
					 b3IsNormalIntoCoveredFace( coveredB, b3Neg( manifold->normal ) ) )
				{
					// The box would be pushed through a face solid covers: answer on its exposed faces
					manifold->pointCount = 0;
					manifold->feature = b3_featureNone;
					b3SATCache scratch = { 0 };
					b3CollideTriangleAndHullFaces( manifold, pointCapacity, vertices[0], vertices[1], vertices[2], triangle.flags,
												   hullB, ~coveredB & 0x3f, &scratch, enableSpeculative );
				}
				break;

			case b3_sphereShape:
				b3CollideTriangleAndSphere( manifold, pointCapacity, vertices, input->sphereB );
				break;

			default:
				B3_ASSERT( false );
				return 0;
		}

		int manifoldPointCount = manifold->pointCount;

		if ( manifoldPointCount > 0 )
		{
			B3_ASSERT( manifold->feature != b3_featureNone );

			manifoldCount += 1;
			totalPointCount += manifoldPointCount;
			manifold->triangleIndex = triangleIndex;
			manifold->triangleNormal = b3MakeNormalFromPoints( vertices[0], vertices[1], vertices[2] );
			manifold->i1 = triangle.i1;
			manifold->i2 = triangle.i2;
			manifold->i3 = triangle.i3;

			if ( manifold->feature == b3_featureTriangleFace || B3_FORCE_GHOST_COLLISIONS )
			{
				(void)b3AddEdge( &foundEdges, triangle.i1, triangle.i2 );
				(void)b3AddEdge( &foundEdges, triangle.i2, triangle.i3 );
				(void)b3AddEdge( &foundEdges, triangle.i3, triangle.i1 );
				(void)b3AddVertex( &foundVertices, triangle.i1 );
				(void)b3AddVertex( &foundVertices, triangle.i2 );
				(void)b3AddVertex( &foundVertices, triangle.i3 );

				acceptedManifolds[acceptedManifoldCount++] = manifold;
			}
			else if ( manifold->feature == b3_featureHullFace )
			{
				float cosNormalAngle = b3Dot( manifold->triangleNormal, manifold->normal );
				if ( cosNormalAngle > 0.5f )
				{
					(void)b3AddEdge( &foundEdges, triangle.i1, triangle.i2 );
					(void)b3AddEdge( &foundEdges, triangle.i2, triangle.i3 );
					(void)b3AddEdge( &foundEdges, triangle.i3, triangle.i1 );
					(void)b3AddVertex( &foundVertices, triangle.i1 );
					(void)b3AddVertex( &foundVertices, triangle.i2 );
					(void)b3AddVertex( &foundVertices, triangle.i3 );

					acceptedManifolds[acceptedManifoldCount++] = manifold;
				}
				else
				{
					float minSeparation = manifold->points[0].separation;
					for ( int i = 1; i < manifoldPointCount; ++i )
					{
						minSeparation = b3MinFloat( minSeparation, manifold->points[i].separation );
					}

					if ( minSeparation < -2.0f * linearSlop )
					{
						// Deep overlap
						(void)b3AddEdge( &foundEdges, triangle.i1, triangle.i2 );
						(void)b3AddEdge( &foundEdges, triangle.i2, triangle.i3 );
						(void)b3AddEdge( &foundEdges, triangle.i3, triangle.i1 );
						(void)b3AddVertex( &foundVertices, triangle.i1 );
						(void)b3AddVertex( &foundVertices, triangle.i2 );
						(void)b3AddVertex( &foundVertices, triangle.i3 );
						acceptedManifolds[acceptedManifoldCount++] = manifold;
					}
					else
					{
						b3TentativeTriangle tentativeTriangle = { .squaredDistance = manifold->squaredDistance,
																  .index = tentativeManifoldCount };
						tentativeTriangles[tentativeTriangleCount++] = tentativeTriangle;
						tentativeManifolds[tentativeManifoldCount++] = manifold;
					}
				}
			}
			else
			{
				b3TentativeTriangle tentativeTriangle = { .squaredDistance = manifold->squaredDistance,
														  .index = tentativeManifoldCount };
				tentativeTriangles[tentativeTriangleCount++] = tentativeTriangle;
				tentativeManifolds[tentativeManifoldCount++] = manifold;
			}
		}
	}

	B3_ASSERT( acceptedManifoldCount <= triangleCount );
	B3_ASSERT( tentativeManifoldCount <= triangleCount );
	B3_ASSERT( tentativeTriangleCount <= triangleCount );

	if ( input->typeB == b3_sphereShape )
	{
		// Sort triangles so the closest triangles are processed first
		{
#define LESS( i, j ) tentativeTriangles[(int)i].squaredDistance < tentativeTriangles[(int)j].squaredDistance
#define SWAP( i, j )                                                                                                             \
	do                                                                                                                           \
	{                                                                                                                            \
		b3TentativeTriangle tmp = tentativeTriangles[(int)i];                                                                    \
		tentativeTriangles[(int)i] = tentativeTriangles[(int)j];                                                                 \
		tentativeTriangles[(int)j] = tmp;                                                                                        \
	}                                                                                                                            \
	while ( 0 )
			QSORT( tentativeTriangleCount, LESS, SWAP );
#undef LESS
#undef SWAP
		}

		// Add tentative manifolds in sorted order. Avoid adding manifolds that generate ghost collisions.
		for ( int i = 0; i < tentativeTriangleCount; ++i )
		{
			b3LocalManifold* m = tentativeManifolds[tentativeTriangles[i].index];

			bool addedEdge1 = b3AddEdge( &foundEdges, m->i1, m->i2 );
			bool addedEdge2 = b3AddEdge( &foundEdges, m->i2, m->i3 );
			bool addedEdge3 = b3AddEdge( &foundEdges, m->i3, m->i1 );
			bool addedVertex1 = b3AddVertex( &foundVertices, m->i1 );
			bool addedVertex2 = b3AddVertex( &foundVertices, m->i2 );
			bool addedVertex3 = b3AddVertex( &foundVertices, m->i3 );

			b3TriangleFeature feature = m->feature;
			bool shouldCollide = false;
			switch ( feature )
			{
				case b3_featureNone:
				case b3_featureTriangleFace:
					B3_ASSERT( false );
					break;

				case b3_featureEdge1:
					shouldCollide = addedEdge1;
					break;

				case b3_featureEdge2:
					shouldCollide = addedEdge2;
					break;

				case b3_featureEdge3:
					shouldCollide = addedEdge3;
					break;

				case b3_featureVertex1:
					shouldCollide = addedVertex1;
					break;

				case b3_featureVertex2:
					shouldCollide = addedVertex2;
					break;

				case b3_featureVertex3:
					shouldCollide = addedVertex3;
					break;

				default:
					B3_ASSERT( false );
					break;
			}

			if ( shouldCollide == true )
			{
				acceptedManifolds[acceptedManifoldCount++] = m;
			}
		}
	}
	else
	{
		// Problem: hull can tunnel if time of impact is at concave edge
		// Example: flat box sliding down a ramp to a flat bottom
		// Solution: only ignore flat edges
		for ( int i = 0; i < tentativeManifoldCount; ++i )
		{
			b3LocalManifold* m = tentativeManifolds[i];
			int triangleFlags = m->triangleFlags;

			if ( ( triangleFlags & b3_allFlatEdges ) == b3_allFlatEdges )
			{
				continue;
			}

			if ( ( triangleFlags & b3_flatEdge1 ) == b3_flatEdge1 )
			{
				if ( b3FindEdge( &foundEdges, m->i1, m->i2 ) )
				{
					continue;
				}
			}

			if ( ( triangleFlags & b3_flatEdge2 ) == b3_flatEdge2 )
			{
				if ( b3FindEdge( &foundEdges, m->i2, m->i3 ) )
				{
					continue;
				}
			}

			if ( ( triangleFlags & b3_flatEdge3 ) == b3_flatEdge3 )
			{
				if ( b3FindEdge( &foundEdges, m->i3, m->i1 ) )
				{
					continue;
				}
			}

			acceptedManifolds[acceptedManifoldCount++] = m;
		}
	}

	B3_ASSERT( acceptedManifoldCount <= triangleCount );
	return acceptedManifoldCount;
}

bool b3ComputeMeshManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA, const int* materialMap,
							 b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast, b3Arena arena )
{
	B3_ASSERT( shapeA->type == b3_meshShape || shapeA->type == b3_heightShape );
	B3_UNUSED( isFast );

	b3TaskContext* context = b3Array_Get( world->taskContexts, workerIndex );

	b3RefreshCache( contact, shapeA, xfA, &shapeB->aabb );

	// Collide with triangles and build manifolds
	b3MeshContact* meshContact = &contact->meshContact;
	int triangleCount = meshContact->triangleCache.count;

	// This should push apart shapes after a time of impact event.
	// In the past I've called this `polygon skin`, but PhysX and Unreal
	// call it `rest offset` which seems appropriate in this case.
	// It leads to a small visual gap but seems to improve the quality of mesh
	// collision, especially for hull versus mesh.
	float restOffset = B3_MESH_REST_OFFSET;

	b3MeshConvexInput input = { 0 };
	input.context = context;
	input.meshShape = shapeA;
	input.typeB = shapeB->type;
	input.hullB = shapeB->type == b3_hullShape ? shapeB->hull : NULL;
	input.capsuleB = &shapeB->capsule;
	input.sphereB = &shapeB->sphere;
	// This transform converts from mesh frame into the shapeB frame
	input.transformAtoB = b3InvMulWorldTransforms( xfB, xfA );
	input.isFast = isFast;
	input.enableSpeculative = contact->flags & b3_enableSpeculativePoints;
	input.coveredB = 0;

	b3LocalManifold** acceptedManifolds = b3Bump( &arena, triangleCount * sizeof( b3LocalManifold* ) );

	// Make room for clip points
	int pointBufferCapacity = B3_MAX_POINTS_PER_TRIANGLE * triangleCount;
	b3LocalManifoldPoint* pointBuffer = b3Bump( &arena, pointBufferCapacity * sizeof( b3LocalManifoldPoint ) );
	b3LocalManifold* manifoldBuffer = b3Bump( &arena, triangleCount * sizeof( b3LocalManifold ) );

	int acceptedManifoldCount = b3CollideMeshTriangles( &input, meshContact->triangleCache.data, triangleCount, acceptedManifolds,
														manifoldBuffer, pointBuffer, pointBufferCapacity, arena );

	if ( b3BuildClusterManifolds( world, contact, acceptedManifolds, acceptedManifoldCount, xfA, xfB, restOffset, arena ) == false )
	{
		return false;
	}

	int clusterCount = contact->manifoldCount;

	const b3SurfaceMaterial* materialsA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );
	b3Vec3 tangentVelocityA = b3Vec3_zero;

	// Update friction and restitution if the mesh has per triangle material
	if ( shapeA->materialCount > 0 )
	{
		float friction = 0.0f;
		float restitution = 0.0f;
		float sampleCount = 0.0f;

		const uint8_t* materialIndices;
		if ( shapeA->type == b3_meshShape )
		{
			materialIndices = b3GetMeshMaterialIndices( shapeA->mesh.data );
		}
		else
		{
			materialIndices = b3GetHeightFieldMaterialIndices( shapeA->heightField );
		}

		for ( int i = 0; i < clusterCount; ++i )
		{
			b3Manifold* manifold = contact->manifolds + i;
			int pointCount = manifold->pointCount;
			for ( int j = 0; j < pointCount; ++j )
			{
				int triangleIndex = manifold->points[j].triangleIndex;
				int materialIndex;
				if ( shapeA->type == b3_meshShape )
				{
					materialIndex = materialIndices[triangleIndex];

					if ( materialMap != NULL )
					{
						materialIndex = materialMap[materialIndex];
					}
				}
				else
				{
					materialIndex = materialIndices[triangleIndex >> 1];
				}

				materialIndex = b3ClampInt( materialIndex, 0, shapeA->materialCount - 1 );
				b3SurfaceMaterial material = materialsA[materialIndex];
				friction += world->frictionCallback( material.friction, material.userMaterialId, materialB->friction,
													 materialB->userMaterialId );
				restitution += world->restitutionCallback( material.restitution, material.userMaterialId, materialB->restitution,
														   materialB->userMaterialId );

				tangentVelocityA = b3Add( tangentVelocityA, material.tangentVelocity );

				sampleCount += 1.0f;
			}
		}

		if ( sampleCount > 0.0f )
		{
			float invCount = 1.0f / sampleCount;
			contact->friction = invCount * friction;
			contact->restitution = invCount * restitution;
			tangentVelocityA = b3MulSV( invCount, tangentVelocityA );
		}

		B3_ASSERT( b3IsValidFloat( contact->friction ) && contact->friction >= 0.0f );
		B3_ASSERT( b3IsValidFloat( contact->restitution ) && contact->restitution >= 0.0f );
	}
	else
	{
		// Keep these updated in case the values on the shapes are modified
		contact->friction = world->frictionCallback( materialsA[0].friction, materialsA[0].userMaterialId, materialB->friction,
													 materialB->userMaterialId );
		contact->restitution = world->restitutionCallback( materialsA[0].restitution, materialsA[0].userMaterialId,
														   materialB->restitution, materialB->userMaterialId );
		tangentVelocityA = materialsA[0].tangentVelocity;
	}

	tangentVelocityA = b3RotateVector( xfA.q, tangentVelocityA );

	float radiusB = 0.0f;
	if ( shapeB->type == b3_sphereShape )
	{
		radiusB = shapeB->sphere.radius;
	}
	else if ( shapeB->type == b3_capsuleShape )
	{
		radiusB = shapeB->capsule.radius;
	}
	else if ( shapeB->type == b3_hullShape )
	{
		radiusB = shapeB->hull->innerRadius;
	}

	contact->rollingResistance = materialB->rollingResistance * radiusB;

	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );
	return true;
}

// Voxel grid contacts
//
// A voxel grid contact is a mesh contact whose cache is keyed by box instead of triangle: the key is ( cell << boxBits ) | box,
// ascending, and each entry holds the SAT or simplex cache of that box against shape B. Each box is collided through the
// convex code on its prebuilt hull. A box with covered faces never pushes through them: see b3CollideHullFaces. The
// manifolds then go through the same clustering as a mesh, so boxes that make one flat surface give one manifold of at most
// four points, and impulses carry over by (feature id, key).

// The bound on boxes in one grid contact. Over it the boxes closest to the other shape are kept.
#define B3_MAX_VOXEL_CONTACT_BOXES 2048

typedef struct b3VoxelBoxCandidate
{
	int key;

	// How far the box is from the other shape's bounds, negative when they overlap
	float gap;
} b3VoxelBoxCandidate;

typedef struct b3VoxelKeyContext
{
	const b3VoxelGrid* grid;
	b3AABB region;
	b3VoxelBoxCandidate* candidates;
	int capacity;
	int count;
	int worst;
	bool overflowed;
} b3VoxelKeyContext;

// The largest separation along an axis between two boxes: negative by the least overlap when they overlap on every axis
static inline float b3BoxGap( b3AABB a, b3AABB b )
{
	float gapX = b3MaxFloat( a.lowerBound.x - b.upperBound.x, b.lowerBound.x - a.upperBound.x );
	float gapY = b3MaxFloat( a.lowerBound.y - b.upperBound.y, b.lowerBound.y - a.upperBound.y );
	float gapZ = b3MaxFloat( a.lowerBound.z - b.upperBound.z, b.lowerBound.z - a.upperBound.z );
	return b3MaxFloat( gapX, b3MaxFloat( gapY, gapZ ) );
}

static bool b3CollectVoxelKeysCallback( int cell, int box, b3AABB bounds, void* context )
{
	b3VoxelKeyContext* keyContext = context;

	// Boxes covered on every face are inside solid and can never be touched
	if ( b3IsVoxelBoxEnclosed( keyContext->grid, cell, box ) )
	{
		return true;
	}

	b3VoxelBoxCandidate candidate = { b3VoxelGridKey( keyContext->grid, cell, box ), b3BoxGap( bounds, keyContext->region ) };
	if ( keyContext->count < keyContext->capacity )
	{
		keyContext->candidates[keyContext->count] = candidate;
		keyContext->count += 1;
		if ( keyContext->count == keyContext->capacity )
		{
			keyContext->worst = 0;
			for ( int i = 1; i < keyContext->count; ++i )
			{
				if ( keyContext->candidates[i].gap > keyContext->candidates[keyContext->worst].gap )
				{
					keyContext->worst = i;
				}
			}
		}
		return true;
	}

	// Full: the box closest to the other shape replaces the farthest one kept
	keyContext->overflowed = true;
	if ( candidate.gap < keyContext->candidates[keyContext->worst].gap )
	{
		keyContext->candidates[keyContext->worst] = candidate;
		keyContext->worst = 0;
		for ( int i = 1; i < keyContext->count; ++i )
		{
			if ( keyContext->candidates[i].gap > keyContext->candidates[keyContext->worst].gap )
			{
				keyContext->worst = i;
			}
		}
	}
	return true;
}

static void b3SortVoxelCandidates( b3VoxelBoxCandidate* candidates, int count )
{
#define LESS( i, j ) candidates[(int)i].key < candidates[(int)j].key
#define SWAP( i, j )                                                                                                             \
	do                                                                                                                           \
	{                                                                                                                            \
		b3VoxelBoxCandidate tmp = candidates[(int)i];                                                                            \
		candidates[(int)i] = candidates[(int)j];                                                                                 \
		candidates[(int)j] = tmp;                                                                                                \
	}                                                                                                                            \
	while ( 0 )
	QSORT( count, LESS, SWAP );
#undef LESS
#undef SWAP
}

// Refresh the boxes the other shape may touch. The query bounds are kept in the grid frame, so a grid that moves refreshes its
// boxes as one that stays does.
static void b3RefreshVoxelCache( b3Contact* contact, const b3VoxelGrid* grid, b3WorldTransform xfA, const b3AABB* bounds,
								 b3Arena arena )
{
	b3MeshContact* meshContact = &contact->meshContact;

	b3Transform gridTransform = b3ToRelativeTransform( xfA, b3Pos_zero );
	b3AABB localBounds = b3AABB_Transform( b3InvertTransform( gridTransform ), *bounds );

	// If the other shape didn't move out of the cached query bounds we are done
	if ( b3AABB_Contains( meshContact->queryBounds, localBounds ) )
	{
		return;
	}

	// Enlarge to the query bounds to absorb small movement
	float radius = B3_MAX_AABB_MARGIN + B3_SPECULATIVE_DISTANCE;
	b3Vec3 extension = { radius, radius, radius };
	meshContact->queryBounds.lowerBound = b3Sub( localBounds.lowerBound, extension );
	meshContact->queryBounds.upperBound = b3Add( localBounds.upperBound, extension );

	b3VoxelBoxCandidate* candidates = b3Bump( &arena, B3_MAX_VOXEL_CONTACT_BOXES * sizeof( b3VoxelBoxCandidate ) );
	b3VoxelKeyContext keyContext = { grid, localBounds, candidates, B3_MAX_VOXEL_CONTACT_BOXES, 0, 0, false };
	b3QueryVoxelGrid( grid, meshContact->queryBounds, b3CollectVoxelKeysCallback, &keyContext );

	if ( keyContext.overflowed )
	{
		static bool s_once = false;
		if ( s_once == false )
		{
			b3Log( "WARNING: dense voxel grid detected, box buffer capacity of %d reached", B3_MAX_VOXEL_CONTACT_BOXES );
			s_once = true;
		}

		// The query walks the cells in order, so only a full buffer that took replacements needs sorting
		b3SortVoxelCandidates( candidates, keyContext.count );
	}

	// Keys are ascending, so match with the old cache by merging
	int count = keyContext.count;
	b3ContactCache* contactCache = b3Bump( &arena, ( count + 1 ) * sizeof( b3ContactCache ) );

	int index2 = 0;
	for ( int index1 = 0; index1 < count; ++index1 )
	{
		contactCache[index1] = (b3ContactCache){ 0 };

		while ( index2 < meshContact->triangleCache.count &&
				meshContact->triangleCache.data[index2].triangleIndex < candidates[index1].key )
		{
			index2 += 1;
		}

		if ( index2 < meshContact->triangleCache.count && meshContact->triangleCache.data[index2].triangleIndex == candidates[index1].key )
		{
			contactCache[index1] = meshContact->triangleCache.data[index2].cache;
		}
	}

	// Save new cache
	b3Array_Resize( meshContact->triangleCache, count );
	for ( int i = 0; i < count; ++i )
	{
		meshContact->triangleCache.data[i] = (b3TriangleCache){ candidates[i].key, contactCache[i] };
	}
}

static bool b3ComputeVoxelGridPairManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
											 b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast,
											 b3Arena arena );

bool b3ComputeVoxelGridManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA, b3WorldTransform xfA,
								  const b3Shape* shapeB, b3WorldTransform xfB, bool isFast, b3Arena arena )
{
	B3_ASSERT( shapeA->type == b3_voxelGridShape );

	if ( shapeB->type == b3_voxelGridShape )
	{
		return b3ComputeVoxelGridPairManifolds( world, workerIndex, contact, shapeA, xfA, shapeB, xfB, isFast, arena );
	}

	const b3VoxelGrid* grid = shapeA->voxelGrid;
	b3TaskContext* context = b3Array_Get( world->taskContexts, workerIndex );

	b3RefreshVoxelCache( contact, grid, xfA, &shapeB->aabb, arena );

	b3MeshContact* meshContact = &contact->meshContact;
	int boxCount = meshContact->triangleCache.count;
	b3TriangleCache* boxCaches = meshContact->triangleCache.data;

	b3LocalManifold** acceptedManifolds = b3Bump( &arena, ( boxCount + 1 ) * sizeof( b3LocalManifold* ) );
	int acceptedManifoldCount = 0;
	b3LocalManifold* manifoldBuffer = b3Bump( &arena, ( boxCount + 1 ) * sizeof( b3LocalManifold ) );
	b3LocalManifoldPoint* pointBuffer = b3Bump( &arena, ( boxCount + 1 ) * B3_MAX_MANIFOLD_POINTS * sizeof( b3LocalManifoldPoint ) );

	// The transform from the grid frame into the frame of B
	b3Transform transformAtoB = b3InvMulWorldTransforms( xfB, xfA );

	const b3HullData* hullB = shapeB->type == b3_hullShape ? shapeB->hull : NULL;

	int currentCell = B3_NULL_INDEX;
	const b3VoxelGridModule* module = NULL;
	b3Transform transformCellToB = b3Transform_identity;
	b3Transform transformBtoCell = b3Transform_identity;

	for ( int index = 0; index < boxCount; ++index )
	{
		int key = boxCaches[index].triangleIndex;
		int cell = b3VoxelGridKeyCell( grid, key );
		int box = b3VoxelGridKeyBox( grid, key );

		if ( cell != currentCell )
		{
			currentCell = cell;
			module = b3GetVoxelGridCellModule( grid, cell );
			B3_ASSERT( module != NULL );

			// Hulls live in the cell frame
			b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );
			transformCellToB = transformAtoB;
			transformCellToB.p = b3Add( transformAtoB.p, b3RotateVector( transformAtoB.q, corner ) );
			transformBtoCell = b3InvertTransform( transformCellToB );
		}

		if ( box >= module->boxCount )
		{
			// The cache is stale after an edit that was not announced. Drop the box instead of reading past the module.
			B3_ASSERT( false );
			continue;
		}

		const b3HullData* hullA = &module->hulls[box].base;
		b3ContactCache* cache = &boxCaches[index].cache;

		b3LocalManifoldPoint points[B3_MAX_POINTS_PER_TRIANGLE];
		b3LocalManifold manifold = { 0 };
		manifold.points = points;
		manifold.feature = b3_featureNone;
		int pointCapacity = B3_MAX_POINTS_PER_TRIANGLE;

		switch ( shapeB->type )
		{
			case b3_capsuleShape:
				b3CollideHullAndCapsule( &manifold, pointCapacity, hullA, &shapeB->capsule, transformBtoCell, &cache->simplexCache );
				break;

			case b3_hullShape:
				// Cached edge contact is dangerous at high speed because the hull can rotate around the edge and tunnel
				// through the box.
				if ( isFast && cache->satCache.type == b3_edgePairAxis )
				{
					cache->satCache = (b3SATCache){ 0 };
				}

				b3CollideHulls( &manifold, pointCapacity, hullA, hullB, transformBtoCell, &cache->satCache );
				context->satCallCount += 1;
				context->satCacheHitCount += cache->satCache.hit;
				break;

			case b3_sphereShape:
				b3CollideHullAndSphere( &manifold, pointCapacity, hullA, &shapeB->sphere, transformBtoCell, &cache->simplexCache );
				break;

			default:
				B3_ASSERT( false );
				return false;
		}

		if ( manifold.pointCount == 0 )
		{
			continue;
		}

		// A contact that would push through a face solid covers is answered on the exposed faces instead
		uint8_t covered = b3GetVoxelBoxCoveredFaces( grid, cell, box );
		if ( covered != 0 && b3IsNormalIntoCoveredFace( covered, manifold.normal ) )
		{
			int exposed = ~covered & 0x3f;
			manifold.pointCount = 0;
			switch ( shapeB->type )
			{
				case b3_capsuleShape:
					b3CollideHullFacesAndCapsule( &manifold, pointCapacity, hullA, exposed, &shapeB->capsule, transformBtoCell );
					break;

				case b3_hullShape:
					b3CollideHullFaces( &manifold, pointCapacity, hullA, exposed, hullB, transformBtoCell, &cache->satCache );
					break;

				default:
					b3CollideHullFacesAndSphere( &manifold, pointCapacity, hullA, exposed, &shapeB->sphere, transformBtoCell );
					break;
			}

			if ( manifold.pointCount == 0 )
			{
				continue;
			}
		}

		// Keep the points in frame B, like a mesh manifold
		b3LocalManifold* accepted = manifoldBuffer + acceptedManifoldCount;
		int pointCount = b3MinInt( manifold.pointCount, B3_MAX_MANIFOLD_POINTS );
		B3_ASSERT( manifold.pointCount <= B3_MAX_MANIFOLD_POINTS );
		*accepted = (b3LocalManifold){ 0 };
		accepted->points = pointBuffer + acceptedManifoldCount * B3_MAX_MANIFOLD_POINTS;
		accepted->pointCount = pointCount;
		accepted->normal = b3RotateVector( transformCellToB.q, manifold.normal );
		accepted->triangleNormal = accepted->normal;
		accepted->triangleIndex = key;
		accepted->feature = b3_featureHullFace;
		for ( int i = 0; i < pointCount; ++i )
		{
			accepted->points[i] = manifold.points[i];
			accepted->points[i].point = b3TransformPoint( transformCellToB, manifold.points[i].point );
		}

		acceptedManifolds[acceptedManifoldCount++] = accepted;
	}

	// No rest offset: a grid is boxes, so bodies rest on it exactly as they rest on a plain box
	if ( b3BuildClusterManifolds( world, contact, acceptedManifolds, acceptedManifoldCount, xfA, xfB, 0.0f, arena ) == false )
	{
		return false;
	}

	// One material for the whole grid
	const b3SurfaceMaterial* materialA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );

	// Keep these updated in case the values on the shapes are modified
	contact->friction =
		world->frictionCallback( materialA->friction, materialA->userMaterialId, materialB->friction, materialB->userMaterialId );
	contact->restitution = world->restitutionCallback( materialA->restitution, materialA->userMaterialId, materialB->restitution,
													   materialB->userMaterialId );

	float radiusB = 0.0f;
	if ( shapeB->type == b3_sphereShape )
	{
		radiusB = shapeB->sphere.radius;
	}
	else if ( shapeB->type == b3_capsuleShape )
	{
		radiusB = shapeB->capsule.radius;
	}
	else if ( shapeB->type == b3_hullShape )
	{
		radiusB = shapeB->hull->innerRadius;
	}

	contact->rollingResistance = materialB->rollingResistance * radiusB;

	b3Vec3 tangentVelocityA = b3RotateVector( xfA.q, materialA->tangentVelocity );
	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );
	return true;
}

// Voxel grid against voxel grid
//
// Both grids are sets of boxes. The boxes of grid A that grid B can reach are listed, and each one asks grid B for the boxes
// within speculative distance of it. A pair of boxes is collided as two hulls. A box covered on every face is left out, and a
// contact that would push through a covered face of either box is answered on the exposed faces of that box, the other box
// being clipped to them. That answer is final: its normal is a face normal of one box and tilted in the frame of the other, so
// it is not tested against the covered faces of the other box. The accepted manifolds go through the same
// clustering as every other grid contact, so one flat landing of two shards is a manifold of at most four points however many
// box pairs touch. Pairs are keyed by a hash of both box keys, kept sorted in the contact cache so the separating axes carry
// over.

// The bound on box pairs in one contact. Over it the pairs closest to touching are kept.
#define B3_MAX_VOXEL_PAIRS 4096

typedef struct b3VoxelSideBox
{
	int key;
	uint8_t covered;

	// In the frame of grid B
	b3AABB bounds;
} b3VoxelSideBox;

typedef struct b3VoxelSideContext
{
	const b3VoxelGrid* grid;
	b3Transform toOtherFrame;
	b3VoxelSideBox* boxes;
	int count;
} b3VoxelSideContext;

// With no buffer this counts the boxes, an upper bound of what a second pass lists
static bool b3CollectVoxelSideCallback( int cell, int box, b3AABB bounds, void* context )
{
	b3VoxelSideContext* side = context;
	if ( side->boxes == NULL )
	{
		side->count += 1;
		return true;
	}

	// A box covered on every face is inside solid and can never be touched
	uint8_t covered = b3GetVoxelBoxCoveredFaces( side->grid, cell, box );
	if ( covered == 0x3f )
	{
		return true;
	}

	b3VoxelSideBox* entry = side->boxes + side->count;
	side->count += 1;
	entry->key = b3VoxelGridKey( side->grid, cell, box );
	entry->covered = covered;
	entry->bounds = b3AABB_Transform( side->toOtherFrame, bounds );
	return true;
}

typedef struct b3VoxelPair
{
	int keyA;
	int keyB;
	int hash;
	float gap;
	uint8_t coveredA;
	uint8_t coveredB;
} b3VoxelPair;

typedef struct b3VoxelPairContext
{
	const b3VoxelGrid* gridB;
	const b3VoxelSideBox* boxA;

	// Null while counting
	b3VoxelPair* pairs;
	int capacity;
	int count;
	int worst;
	bool overflowed;
} b3VoxelPairContext;

static void b3FindWorstVoxelPair( b3VoxelPairContext* pairContext )
{
	pairContext->worst = 0;
	for ( int i = 1; i < pairContext->count; ++i )
	{
		if ( pairContext->pairs[i].gap > pairContext->pairs[pairContext->worst].gap )
		{
			pairContext->worst = i;
		}
	}
}

static bool b3CollectVoxelPairCallback( int cell, int box, b3AABB bounds, void* context )
{
	b3VoxelPairContext* pairContext = context;

	float gap = b3BoxGap( pairContext->boxA->bounds, bounds );
	if ( gap > B3_SPECULATIVE_DISTANCE )
	{
		return true;
	}

	// Faces that touch overlap on two axes. Boxes that only meet along an edge or at a corner, such as the neighbours of the
	// box under a flat landing, are left to the boxes that cover them: they would only add ghost contacts to the solve.
	const b3AABB* a = &pairContext->boxA->bounds;
	float tolerance = 0.5f * B3_LINEAR_SLOP;
	int overlapAxes = ( b3MinFloat( a->upperBound.x, bounds.upperBound.x ) - b3MaxFloat( a->lowerBound.x, bounds.lowerBound.x ) > tolerance ) +
					  ( b3MinFloat( a->upperBound.y, bounds.upperBound.y ) - b3MaxFloat( a->lowerBound.y, bounds.lowerBound.y ) > tolerance ) +
					  ( b3MinFloat( a->upperBound.z, bounds.upperBound.z ) - b3MaxFloat( a->lowerBound.z, bounds.lowerBound.z ) > tolerance );
	if ( overlapAxes < 2 )
	{
		return true;
	}

	uint8_t covered = b3GetVoxelBoxCoveredFaces( pairContext->gridB, cell, box );
	if ( covered == 0x3f )
	{
		return true;
	}

	if ( pairContext->pairs == NULL )
	{
		pairContext->count = b3MinInt( pairContext->count + 1, pairContext->capacity );
		return true;
	}

	b3VoxelPair pair = { pairContext->boxA->key, b3VoxelGridKey( pairContext->gridB, cell, box ), 0, gap, pairContext->boxA->covered,
						 covered };
	if ( pairContext->count < pairContext->capacity )
	{
		pairContext->pairs[pairContext->count] = pair;
		pairContext->count += 1;
		if ( pairContext->count == pairContext->capacity )
		{
			b3FindWorstVoxelPair( pairContext );
		}
	}
	else
	{
		// Full: the pair closest to touching replaces the farthest one kept
		pairContext->overflowed = true;
		if ( gap < pairContext->pairs[pairContext->worst].gap )
		{
			pairContext->pairs[pairContext->worst] = pair;
			b3FindWorstVoxelPair( pairContext );
		}
	}

	return true;
}

static inline int b3VoxelPairHash( int keyA, int keyB )
{
	uint32_t h = (uint32_t)keyA * 0x9E3779B1u;
	h ^= (uint32_t)keyB + 0x7F4A7C15u + ( h << 6 ) + ( h >> 2 );
	h *= 0x85EBCA6Bu;
	h ^= h >> 15;
	return (int)( h & 0x7fffffffu );
}

// Collide two boxes of voxel grids. The manifold is in the frame of box A with the normal from A to B.
static void b3CollideVoxelBoxes( b3LocalManifold* manifold, int capacity, const b3HullData* hullA, uint8_t coveredA,
								 const b3HullData* hullB, uint8_t coveredB, b3Transform transformBtoA, b3SATCache* cache )
{
	b3CollideHulls( manifold, capacity, hullA, hullB, transformBtoA, cache );
	if ( manifold->pointCount == 0 || ( coveredA == 0 && coveredB == 0 ) )
	{
		return;
	}

	// The outward direction of box B toward A is checked in the frame of B
	b3Vec3 normal = manifold->normal;
	bool intoA = coveredA != 0 && b3IsNormalIntoCoveredFace( coveredA, normal );
	bool intoB = coveredB != 0 && b3IsNormalIntoCoveredFace( coveredB, b3InvRotateVector( transformBtoA.q, b3Neg( normal ) ) );
	if ( intoA == false && intoB == false )
	{
		return;
	}

	b3SATCache scratch = { 0 };
	if ( intoA )
	{
		manifold->pointCount = 0;
		b3CollideHullFaces( manifold, capacity, hullA, ~coveredA & 0x3f, hullB, transformBtoA, &scratch );
		return;
	}

	// Only the covered face of B is in the way: answer on the exposed faces of B, clipping A to them
	b3LocalManifoldPoint points[B3_MAX_POINTS_PER_TRIANGLE];
	b3LocalManifold swapped = { 0 };
	swapped.points = points;
	b3Transform transformAtoB = b3InvertTransform( transformBtoA );
	b3CollideHullFaces( &swapped, B3_MAX_POINTS_PER_TRIANGLE, hullB, ~coveredB & 0x3f, hullA, transformAtoB, &scratch );

	manifold->pointCount = 0;
	if ( swapped.pointCount == 0 )
	{
		return;
	}

	manifold->normal = b3Neg( b3RotateVector( transformBtoA.q, swapped.normal ) );
	manifold->pointCount = swapped.pointCount;
	for ( int i = 0; i < swapped.pointCount; ++i )
	{
		manifold->points[i] = swapped.points[i];
		manifold->points[i].point = b3TransformPoint( transformBtoA, swapped.points[i].point );
		manifold->points[i].pair = b3FlipPair( swapped.points[i].pair );
	}
}

// Point-sampled voxel grid against voxel grid (b3World_EnableVoxelPointContacts)
//
// Each grid's exposed box faces in the region where the two can touch are sampled at nine points (corners, edge midpoints
// and centre, inset a little), and each sample is tested against the other grid's voxel occupancy. A sample whose outside
// neighbour voxel is solid in its own grid lies on an internal seam and is skipped. A sample inside a solid voxel of the other
// grid walks that grid's six axes to the nearest empty voxel (up to B3_VOXEL_SAMPLE_WALK_METERS): the shortest exit gives the
// depth and the normal. A sample whose normal is more than 60 degrees from its face's outward direction (a side face that
// sank, a face that went through a thin wall, the far face of a thin plate inside a block) or that finds no exit gives
// nothing: the touching faces' samples answer. A sample in empty space probes one speculative distance along its face's outward
// direction and gives a speculative point on the voxel it finds. Normals are axes of the grid sampled against, so the points
// cluster into at most twelve flat manifolds (six per grid frame), each reduced to four points by the shared clustering.
// A point's id is a hash of its side, box key and sample, so impulses carry over while a body rests.

// The bound on contact points from one grid pair before reduction. Over it the deepest are kept.
#define B3_MAX_VOXEL_SAMPLE_POINTS 8192

// How far a sample inside the other grid looks for the surface
#define B3_VOXEL_SAMPLE_WALK_METERS 1.0f

// A sample answers only along normals within 60 degrees of its face's outward direction: a side face that sank into the
// other grid, or a thin slab's face, would otherwise take a short exit across itself and push the wrong way
#define B3_VOXEL_SAMPLE_MIN_FACING 0.5f

typedef struct b3VoxelSampler
{
	const b3VoxelGrid* self;
	const b3VoxelGrid* other;

	// The body frame of self into the body frame of other
	b3Transform selfToOther;

	// The occupied bounds of other, grown by the speculative distance, in its frame
	b3AABB otherBounds;

	// Self is shape B of the contact: the points are flipped into the contact's convention afterwards
	bool selfIsB;

	b3LocalManifoldPoint* points;
	b3Vec3* normals;
	int capacity;
	int count;
	int shallowest;
	bool overflowed;
	int sampleCount;
} b3VoxelSampler;

static inline int b3VoxelSampleHash( int side, int key, int sample )
{
	uint32_t h = (uint32_t)key * 0x9E3779B1u;
	h ^= (uint32_t)( sample * 2 + side ) + 0x7F4A7C15u + ( h << 6 ) + ( h >> 2 );
	h *= 0x85EBCA6Bu;
	h ^= h >> 15;
	return (int)( h & 0x7fffffffu );
}

static void b3FindShallowestSample( b3VoxelSampler* sampler )
{
	sampler->shallowest = 0;
	for ( int i = 1; i < sampler->count; ++i )
	{
		if ( sampler->points[i].separation > sampler->points[sampler->shallowest].separation )
		{
			sampler->shallowest = i;
		}
	}
}

static void b3AddVoxelSample( b3VoxelSampler* sampler, b3Vec3 point, b3Vec3 normal, float separation, int id )
{
	b3LocalManifoldPoint mp = { 0 };
	mp.point = point;
	mp.separation = separation;
	mp.pair = b3FeaturePair_single;
	mp.triangleIndex = id;

	if ( sampler->count < sampler->capacity )
	{
		sampler->points[sampler->count] = mp;
		sampler->normals[sampler->count] = normal;
		sampler->count += 1;
		if ( sampler->count == sampler->capacity )
		{
			b3FindShallowestSample( sampler );
		}
		return;
	}

	sampler->overflowed = true;
	if ( separation < sampler->points[sampler->shallowest].separation )
	{
		sampler->points[sampler->shallowest] = mp;
		sampler->normals[sampler->shallowest] = normal;
		b3FindShallowestSample( sampler );
	}
}

static inline void b3VoxelCoordinates( const b3VoxelGrid* grid, float invVoxel, b3Vec3 p, int v[3] )
{
	v[0] = (int)floorf( ( p.x - grid->origin.x ) * invVoxel );
	v[1] = (int)floorf( ( p.y - grid->origin.y ) * invVoxel );
	v[2] = (int)floorf( ( p.z - grid->origin.z ) * invVoxel );
}

// One sample: p is on a face of self, just inside it, with the face's outward direction, both in the frame of self. The
// point goes in the frame of other with the normal from self to other.
static void b3SampleVoxelPoint( b3VoxelSampler* sampler, b3Vec3 p, b3Vec3 faceNormal, int id )
{
	const b3VoxelGrid* other = b3VoxelGridOccupancy( sampler->other );
	float voxel = other->cellMeters / (float)other->cellVoxels;
	float invVoxel = 1.0f / voxel;

	b3Vec3 q = b3TransformPoint( sampler->selfToOther, p );
	if ( q.x < sampler->otherBounds.lowerBound.x || q.y < sampler->otherBounds.lowerBound.y || q.z < sampler->otherBounds.lowerBound.z ||
		 q.x > sampler->otherBounds.upperBound.x || q.y > sampler->otherBounds.upperBound.y || q.z > sampler->otherBounds.upperBound.z )
	{
		return;
	}

	sampler->sampleCount += 1;
	b3Vec3 outward = b3RotateVector( sampler->selfToOther.q, faceNormal );
	const float qa[3] = { q.x, q.y, q.z };
	const float origin[3] = { other->origin.x, other->origin.y, other->origin.z };

	int v[3];
	b3VoxelCoordinates( other, invVoxel, q, v );

	if ( b3IsVoxelSolid( other, v[0], v[1], v[2], false ) )
	{
		// Inside: the nearest exit along the six axes
		int maxSteps = (int)ceilf( B3_VOXEL_SAMPLE_WALK_METERS * invVoxel );
		float best = FLT_MAX;
		int bestAxis = -1, bestSign = 0;
		for ( int axis = 0; axis < 3; ++axis )
		{
			for ( int sign = -1; sign <= 1; sign += 2 )
			{
				int w[3] = { v[0], v[1], v[2] };
				for ( int step = 1; step <= maxSteps; ++step )
				{
					float lowerBound = ( sign > 0 ? (float)( v[axis] + step ) * voxel + origin[axis] - qa[axis]
												  : qa[axis] - ( (float)( v[axis] - step + 1 ) * voxel + origin[axis] ) );
					if ( lowerBound >= best )
					{
						break;
					}

					w[axis] = v[axis] + sign * step;
					if ( b3IsVoxelSolid( other, w[0], w[1], w[2], true ) == false )
					{
						if ( lowerBound < best )
						{
							best = lowerBound;
							bestAxis = axis;
							bestSign = sign;
						}
						break;
					}
				}
			}
		}

		// Buried past the walk: no answer from this sample
		if ( bestAxis < 0 )
		{
			return;
		}

		// Moving self along the exit pushes other the opposite way: the normal from self to other is -exit. A normal that
		// points back into self's face means this face is not the one touching (a face that went through a thin wall, the
		// far face of a thin plate inside a block); the samples of the touching faces answer instead.
		float n[3] = { 0.0f, 0.0f, 0.0f };
		n[bestAxis] = -(float)bestSign;
		b3Vec3 normal = { n[0], n[1], n[2] };
		if ( b3Dot( normal, outward ) < B3_VOXEL_SAMPLE_MIN_FACING )
		{
			return;
		}
		float depth = best;

		b3Vec3 point = b3MulSub( q, 0.5f * depth, normal );
		b3AddVoxelSample( sampler, point, normal, -depth, id );
		return;
	}

	// Outside: the voxel one speculative distance out along the face
	float reach = B3_SPECULATIVE_DISTANCE + 0.05f * voxel;
	b3Vec3 probe = b3MulAdd( q, reach, outward );
	int u[3];
	b3VoxelCoordinates( other, invVoxel, probe, u );
	if ( ( u[0] == v[0] && u[1] == v[1] && u[2] == v[2] ) || b3IsVoxelSolid( other, u[0], u[1], u[2], false ) == false )
	{
		return;
	}

	// The closest point of that voxel; the normal is the axis the gap is largest along
	float gap[3];
	int axis = 0;
	for ( int a = 0; a < 3; ++a )
	{
		float low = (float)u[a] * voxel + origin[a], high = low + voxel;
		float closest = b3ClampFloat( qa[a], low, high );
		gap[a] = closest - qa[a];
		if ( fabsf( gap[a] ) > fabsf( gap[axis] ) )
		{
			axis = a;
		}
	}

	float n[3] = { 0.0f, 0.0f, 0.0f };
	n[axis] = gap[axis] >= 0.0f ? 1.0f : -1.0f;
	b3Vec3 normal = { n[0], n[1], n[2] };
	if ( b3Dot( normal, outward ) < B3_VOXEL_SAMPLE_MIN_FACING )
	{
		return;
	}

	float separation = fabsf( gap[axis] ) - 0.05f * voxel;
	if ( separation > B3_SPECULATIVE_DISTANCE )
	{
		return;
	}

	b3Vec3 point = b3MulAdd( q, 0.5f * fabsf( gap[axis] ), normal );
	b3AddVoxelSample( sampler, point, normal, separation, id );
}

static bool b3FindAnyVoxelBoxCallback( int cell, int box, b3AABB bounds, void* context )
{
	B3_UNUSED( cell, box, bounds );
	*(bool*)context = true;
	return false;
}

static bool b3SampleVoxelBoxCallback( int cell, int box, b3AABB bounds, void* context )
{
	b3VoxelSampler* sampler = context;
	const b3VoxelGrid* self = sampler->self;

	// Only a box the other grid has boxes near is sampled
	b3AABB boundsInOther = b3AABB_Transform( sampler->selfToOther, bounds );
	if ( b3AABB_Overlaps( boundsInOther, sampler->otherBounds ) == false )
	{
		return true;
	}

	boundsInOther.lowerBound = b3Sub( boundsInOther.lowerBound, ( b3Vec3 ){ B3_SPECULATIVE_DISTANCE, B3_SPECULATIVE_DISTANCE, B3_SPECULATIVE_DISTANCE } );
	boundsInOther.upperBound = b3Add( boundsInOther.upperBound, ( b3Vec3 ){ B3_SPECULATIVE_DISTANCE, B3_SPECULATIVE_DISTANCE, B3_SPECULATIVE_DISTANCE } );
	bool near = false;
	b3QueryVoxelGrid( sampler->other, boundsInOther, b3FindAnyVoxelBoxCallback, &near );
	if ( near == false )
	{
		return true;
	}

	uint8_t covered = b3GetVoxelBoxCoveredFaces( self, cell, box );
	if ( covered == 0x3f )
	{
		return true;
	}

	const b3VoxelGrid* occupancy = b3VoxelGridOccupancy( self );
	float voxel = occupancy->cellMeters / (float)occupancy->cellVoxels;
	float invVoxel = 1.0f / voxel;
	float inset = 0.05f * voxel;
	int key = b3VoxelGridKey( self, cell, box );
	const float low[3] = { bounds.lowerBound.x, bounds.lowerBound.y, bounds.lowerBound.z };
	const float high[3] = { bounds.upperBound.x, bounds.upperBound.y, bounds.upperBound.z };

	for ( int face = 0; face < 6; ++face )
	{
		if ( covered & ( 1 << face ) )
		{
			continue;
		}

		int axis = face >> 1;
		float sign = ( face & 1 ) ? 1.0f : -1.0f;
		int ua = ( axis + 1 ) % 3, va = ( axis + 2 ) % 3;
		float plane = ( face & 1 ) ? high[axis] : low[axis];
		float n[3] = { 0.0f, 0.0f, 0.0f };
		n[axis] = sign;
		b3Vec3 faceNormal = { n[0], n[1], n[2] };

		float us[3] = { low[ua] + inset, 0.5f * ( low[ua] + high[ua] ), high[ua] - inset };
		float vs[3] = { low[va] + inset, 0.5f * ( low[va] + high[va] ), high[va] - inset };
		for ( int iv = 0; iv < 3; ++iv )
		{
			for ( int iu = 0; iu < 3; ++iu )
			{
				float s[3];
				s[axis] = plane;
				s[ua] = us[iu];
				s[va] = vs[iv];

				// A sample whose outside neighbour is solid in its own grid is on an internal seam
				float outside[3] = { s[0], s[1], s[2] };
				outside[axis] += 0.5f * sign * voxel;
				int w[3];
				b3VoxelCoordinates( occupancy, invVoxel, ( b3Vec3 ){ outside[0], outside[1], outside[2] }, w );
				if ( b3IsVoxelSolid( occupancy, w[0], w[1], w[2], true ) )
				{
					continue;
				}

				s[axis] -= sign * inset;
				int id = b3VoxelSampleHash( sampler->selfIsB ? 1 : 0, key, face * 9 + iv * 3 + iu );
				b3SampleVoxelPoint( sampler, ( b3Vec3 ){ s[0], s[1], s[2] }, faceNormal, id );
			}
		}
	}
	return true;
}

static bool b3ComputeVoxelGridSampledManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
												b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, b3Arena arena )
{
	const b3VoxelGrid* gridA = shapeA->voxelGrid;
	const b3VoxelGrid* gridB = shapeB->voxelGrid;
	b3TaskContext* context = b3Array_Get( world->taskContexts, workerIndex );
	b3MeshContact* meshContact = &contact->meshContact;

	// No box pair caches on this path
	b3Array_Resize( meshContact->triangleCache, 0 );

	b3Transform transformBtoA = b3InvMulWorldTransforms( xfA, xfB );
	b3Transform transformAtoB = b3InvMulWorldTransforms( xfB, xfA );

	float margin = 2.0f * B3_SPECULATIVE_DISTANCE;
	b3AABB region = { b3Max( shapeA->aabb.lowerBound, shapeB->aabb.lowerBound ),
					  b3Min( shapeA->aabb.upperBound, shapeB->aabb.upperBound ) };
	bool overlap = region.lowerBound.x <= region.upperBound.x && region.lowerBound.y <= region.upperBound.y &&
				   region.lowerBound.z <= region.upperBound.z;

	b3LocalManifoldPoint* points = b3Bump( &arena, B3_MAX_VOXEL_SAMPLE_POINTS * (int)sizeof( b3LocalManifoldPoint ) );
	b3Vec3* normals = b3Bump( &arena, B3_MAX_VOXEL_SAMPLE_POINTS * (int)sizeof( b3Vec3 ) );
	int countA = 0;
	int totalSamples = 0;
	bool overflowed = false;

	if ( overlap )
	{
		region.lowerBound = b3Sub( region.lowerBound, ( b3Vec3 ){ margin, margin, margin } );
		region.upperBound = b3Add( region.upperBound, ( b3Vec3 ){ margin, margin, margin } );
		b3Vec3 grow = { B3_SPECULATIVE_DISTANCE, B3_SPECULATIVE_DISTANCE, B3_SPECULATIVE_DISTANCE };

		// A's faces against B, giving points in frame B with the normal from A to B
		b3VoxelSampler sampler = { 0 };
		sampler.self = gridA;
		sampler.other = gridB;
		sampler.selfToOther = transformAtoB;
		sampler.otherBounds = (b3AABB){ b3Sub( gridB->occupiedBounds.lowerBound, grow ), b3Add( gridB->occupiedBounds.upperBound, grow ) };
		sampler.selfIsB = false;
		sampler.points = points;
		sampler.normals = normals;
		sampler.capacity = B3_MAX_VOXEL_SAMPLE_POINTS;

		b3Transform gridTransformA = b3ToRelativeTransform( xfA, b3Pos_zero );
		b3Transform gridTransformB = b3ToRelativeTransform( xfB, b3Pos_zero );

		// Each side samples only where the other has boxes
		b3AABB occupiedA = b3AABB_Transform( gridTransformA, gridA->occupiedBounds );
		b3AABB occupiedB = b3AABB_Transform( gridTransformB, gridB->occupiedBounds );
		b3AABB regionForA = { b3Max( region.lowerBound, b3Sub( occupiedB.lowerBound, grow ) ), b3Min( region.upperBound, b3Add( occupiedB.upperBound, grow ) ) };
		b3AABB regionForB = { b3Max( region.lowerBound, b3Sub( occupiedA.lowerBound, grow ) ), b3Min( region.upperBound, b3Add( occupiedA.upperBound, grow ) ) };
		bool sampleA = regionForA.lowerBound.x <= regionForA.upperBound.x && regionForA.lowerBound.y <= regionForA.upperBound.y &&
					   regionForA.lowerBound.z <= regionForA.upperBound.z;
		bool sampleB = regionForB.lowerBound.x <= regionForB.upperBound.x && regionForB.lowerBound.y <= regionForB.upperBound.y &&
					   regionForB.lowerBound.z <= regionForB.upperBound.z;

		if ( sampleA )
		{
			b3QueryVoxelGrid( gridA, b3AABB_Transform( b3InvertTransform( gridTransformA ), regionForA ), b3SampleVoxelBoxCallback, &sampler );
		}
		countA = sampler.count;
		totalSamples += sampler.sampleCount;
		overflowed = sampler.overflowed;

		// B's faces against A, in frame A with the normal from B to A, flipped below. They share the buffer.
		b3VoxelSampler samplerB = { 0 };
		samplerB.self = gridB;
		samplerB.other = gridA;
		samplerB.selfToOther = transformBtoA;
		samplerB.otherBounds = (b3AABB){ b3Sub( gridA->occupiedBounds.lowerBound, grow ), b3Add( gridA->occupiedBounds.upperBound, grow ) };
		samplerB.selfIsB = true;
		samplerB.points = points + countA;
		samplerB.normals = normals + countA;
		samplerB.capacity = B3_MAX_VOXEL_SAMPLE_POINTS - countA;

		if ( samplerB.capacity > 0 && sampleB )
		{
			b3QueryVoxelGrid( gridB, b3AABB_Transform( b3InvertTransform( gridTransformB ), regionForB ), b3SampleVoxelBoxCallback, &samplerB );
		}

		for ( int i = 0; i < samplerB.count; ++i )
		{
			b3LocalManifoldPoint* mp = points + countA + i;
			mp->point = b3TransformPoint( transformAtoB, mp->point );
			normals[countA + i] = b3Neg( b3RotateVector( transformAtoB.q, normals[countA + i] ) );
		}

		countA += samplerB.count;
		totalSamples += samplerB.sampleCount;
		overflowed = overflowed || samplerB.overflowed;
	}

	if ( overflowed )
	{
		static bool s_once = false;
		if ( s_once == false )
		{
			b3Log( "WARNING: dense voxel grid contact, sample point capacity of %d reached", B3_MAX_VOXEL_SAMPLE_POINTS );
			s_once = true;
		}
	}

	int pointCount = countA;
	context->voxelGridPairCount += 1;
	context->voxelSampleCount += totalSamples;
	context->voxelPointCount += pointCount;

	b3LocalManifold** acceptedManifolds = b3Bump( &arena, ( pointCount + 1 ) * (int)sizeof( b3LocalManifold* ) );
	b3LocalManifold* manifoldBuffer = b3Bump( &arena, ( pointCount + 1 ) * (int)sizeof( b3LocalManifold ) );
	for ( int i = 0; i < pointCount; ++i )
	{
		context->voxelMinSeparation = b3MinFloat( context->voxelMinSeparation, points[i].separation );

		b3LocalManifold* accepted = manifoldBuffer + i;
		*accepted = (b3LocalManifold){ 0 };
		accepted->points = points + i;
		accepted->pointCount = 1;
		accepted->normal = normals[i];
		accepted->triangleNormal = normals[i];
		accepted->triangleIndex = points[i].triangleIndex;
		accepted->feature = b3_featureHullFace;
		acceptedManifolds[i] = accepted;
	}

	if ( b3BuildClusterManifolds( world, contact, acceptedManifolds, pointCount, xfA, xfB, 0.0f, arena ) == false )
	{
		return false;
	}

	const b3SurfaceMaterial* materialA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );
	contact->friction =
		world->frictionCallback( materialA->friction, materialA->userMaterialId, materialB->friction, materialB->userMaterialId );
	contact->restitution = world->restitutionCallback( materialA->restitution, materialA->userMaterialId, materialB->restitution,
													   materialB->userMaterialId );
	contact->rollingResistance = 0.0f;

	b3Vec3 tangentVelocityA = b3RotateVector( xfA.q, materialA->tangentVelocity );
	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );
	return true;
}

static bool b3ComputeVoxelGridPairManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
											 b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast,
											 b3Arena arena )
{
	B3_ASSERT( shapeA->type == b3_voxelGridShape && shapeB->type == b3_voxelGridShape );

	if ( world->enableVoxelPointContacts )
	{
		return b3ComputeVoxelGridSampledManifolds( world, workerIndex, contact, shapeA, xfA, shapeB, xfB, arena );
	}

	const b3VoxelGrid* gridA = shapeA->voxelGrid;
	const b3VoxelGrid* gridB = shapeB->voxelGrid;
	b3TaskContext* context = b3Array_Get( world->taskContexts, workerIndex );
	b3MeshContact* meshContact = &contact->meshContact;

	// Grid B in the frame of grid A, and grid A in the frame of grid B
	b3Transform transformBtoA = b3InvMulWorldTransforms( xfA, xfB );
	b3Transform transformAtoB = b3InvMulWorldTransforms( xfB, xfA );

	// The region where the two can touch, in the frame of each grid
	float margin = 2.0f * B3_SPECULATIVE_DISTANCE;
	b3AABB region = { b3Max( shapeA->aabb.lowerBound, shapeB->aabb.lowerBound ),
					  b3Min( shapeA->aabb.upperBound, shapeB->aabb.upperBound ) };
	bool overlap = region.lowerBound.x <= region.upperBound.x && region.lowerBound.y <= region.upperBound.y &&
				   region.lowerBound.z <= region.upperBound.z;

	// The boxes of grid A in the region, with their bounds in the frame of grid B
	b3VoxelSideBox* boxesA = NULL;
	b3VoxelSideContext sideA = { gridA, transformAtoB, NULL, 0 };
	if ( overlap )
	{
		region.lowerBound = b3Sub( region.lowerBound, ( b3Vec3 ){ margin, margin, margin } );
		region.upperBound = b3Add( region.upperBound, ( b3Vec3 ){ margin, margin, margin } );

		b3Transform gridTransformA = b3ToRelativeTransform( xfA, b3Pos_zero );
		b3AABB regionA = b3AABB_Transform( b3InvertTransform( gridTransformA ), region );

		b3QueryVoxelGrid( gridA, regionA, b3CollectVoxelSideCallback, &sideA );
		boxesA = b3Bump( &arena, sideA.count * sizeof( b3VoxelSideBox ) );
		sideA.boxes = boxesA;
		sideA.count = 0;
		b3QueryVoxelGrid( gridA, regionA, b3CollectVoxelSideCallback, &sideA );
	}

	// The pairs within speculative distance, counted first so that the memory is what the pairs need
	b3VoxelPair* pairs = NULL;
	int pairCount = 0;
	bool overflowed = false;
	if ( sideA.count > 0 )
	{
		b3VoxelPairContext pairContext = { gridB, NULL, NULL, B3_MAX_VOXEL_PAIRS, 0, 0, false };
		for ( int a = 0; a < sideA.count; ++a )
		{
			pairContext.boxA = boxesA + a;
			b3AABB query = pairContext.boxA->bounds;
			query.lowerBound = b3Sub( query.lowerBound, ( b3Vec3 ){ margin, margin, margin } );
			query.upperBound = b3Add( query.upperBound, ( b3Vec3 ){ margin, margin, margin } );
			b3QueryVoxelGrid( gridB, query, b3CollectVoxelPairCallback, &pairContext );
		}

		if ( pairContext.count > 0 )
		{
			pairs = b3Bump( &arena, pairContext.count * sizeof( b3VoxelPair ) );
			pairContext.pairs = pairs;
			pairContext.capacity = pairContext.count;
			pairContext.count = 0;
			for ( int a = 0; a < sideA.count; ++a )
			{
				pairContext.boxA = boxesA + a;
				b3AABB query = pairContext.boxA->bounds;
				query.lowerBound = b3Sub( query.lowerBound, ( b3Vec3 ){ margin, margin, margin } );
				query.upperBound = b3Add( query.upperBound, ( b3Vec3 ){ margin, margin, margin } );
				b3QueryVoxelGrid( gridB, query, b3CollectVoxelPairCallback, &pairContext );
			}

			pairCount = pairContext.count;
			overflowed = pairContext.overflowed;
		}
	}

	if ( overflowed )
	{
		static bool s_once = false;
		if ( s_once == false )
		{
			b3Log( "WARNING: dense voxel grid contact, box pair capacity of %d reached", B3_MAX_VOXEL_PAIRS );
			s_once = true;
		}
	}

	// Sort by hash to match the previous cache by merging
	for ( int i = 0; i < pairCount; ++i )
	{
		pairs[i].hash = b3VoxelPairHash( pairs[i].keyA, pairs[i].keyB );
	}

	{
#define LESS( i, j ) pairs[(int)i].hash < pairs[(int)j].hash
#define SWAP( i, j )                                                                                                             \
	do                                                                                                                           \
	{                                                                                                                            \
		b3VoxelPair tmp = pairs[(int)i];                                                                                         \
		pairs[(int)i] = pairs[(int)j];                                                                                           \
		pairs[(int)j] = tmp;                                                                                                     \
	}                                                                                                                            \
	while ( 0 )
		QSORT( pairCount, LESS, SWAP );
#undef LESS
#undef SWAP
	}

	b3ContactCache* caches = b3Bump( &arena, pairCount * sizeof( b3ContactCache ) );
	{
		int index2 = 0;
		for ( int index1 = 0; index1 < pairCount; ++index1 )
		{
			caches[index1] = (b3ContactCache){ 0 };
			while ( index2 < meshContact->triangleCache.count &&
					meshContact->triangleCache.data[index2].triangleIndex < pairs[index1].hash )
			{
				index2 += 1;
			}

			if ( index2 < meshContact->triangleCache.count && meshContact->triangleCache.data[index2].triangleIndex == pairs[index1].hash )
			{
				caches[index1] = meshContact->triangleCache.data[index2].cache;
			}
		}
	}

	b3LocalManifold** acceptedManifolds = b3Bump( &arena, pairCount * sizeof( b3LocalManifold* ) );
	int acceptedManifoldCount = 0;
	b3LocalManifold* manifoldBuffer = b3Bump( &arena, pairCount * sizeof( b3LocalManifold ) );
	b3LocalManifoldPoint* pointBuffer = b3Bump( &arena, pairCount * B3_MAX_MANIFOLD_POINTS * sizeof( b3LocalManifoldPoint ) );

	for ( int index = 0; index < pairCount; ++index )
	{
		const b3VoxelPair* pair = pairs + index;
		int cellA = b3VoxelGridKeyCell( gridA, pair->keyA ), boxIndexA = b3VoxelGridKeyBox( gridA, pair->keyA );
		int cellB = b3VoxelGridKeyCell( gridB, pair->keyB ), boxIndexB = b3VoxelGridKeyBox( gridB, pair->keyB );
		const b3VoxelGridModule* moduleA = b3GetVoxelGridCellModule( gridA, cellA );
		const b3VoxelGridModule* moduleB = b3GetVoxelGridCellModule( gridB, cellB );
		B3_ASSERT( moduleA != NULL && moduleB != NULL && boxIndexA < moduleA->boxCount && boxIndexB < moduleB->boxCount );

		b3Vec3 cornerA = b3VoxelGridCellCorner( gridA, cellA );
		b3Vec3 cornerB = b3VoxelGridCellCorner( gridB, cellB );

		// The frame of box B in the frame of box A
		b3Transform transformCellBtoCellA;
		transformCellBtoCellA.q = transformBtoA.q;
		transformCellBtoCellA.p = b3Sub( b3Add( b3RotateVector( transformBtoA.q, cornerB ), transformBtoA.p ), cornerA );

		b3ContactCache* cache = caches + index;
		if ( isFast && cache->satCache.type == b3_edgePairAxis )
		{
			cache->satCache = (b3SATCache){ 0 };
		}

		b3LocalManifoldPoint points[B3_MAX_POINTS_PER_TRIANGLE];
		b3LocalManifold manifold = { 0 };
		manifold.points = points;
		manifold.feature = b3_featureNone;

		b3CollideVoxelBoxes( &manifold, B3_MAX_POINTS_PER_TRIANGLE, &moduleA->hulls[boxIndexA].base, pair->coveredA,
							 &moduleB->hulls[boxIndexB].base, pair->coveredB, transformCellBtoCellA, &cache->satCache );
		context->satCallCount += 1;
		context->satCacheHitCount += cache->satCache.hit;

		if ( manifold.pointCount == 0 )
		{
			continue;
		}

		// Keep the points in frame B, the body frame of grid B
		b3LocalManifold* accepted = manifoldBuffer + acceptedManifoldCount;
		int pointCount = b3MinInt( manifold.pointCount, B3_MAX_MANIFOLD_POINTS );
		B3_ASSERT( manifold.pointCount <= B3_MAX_MANIFOLD_POINTS );
		*accepted = (b3LocalManifold){ 0 };
		accepted->points = pointBuffer + acceptedManifoldCount * B3_MAX_MANIFOLD_POINTS;
		accepted->pointCount = pointCount;
		accepted->normal = b3RotateVector( transformAtoB.q, manifold.normal );
		accepted->triangleNormal = accepted->normal;
		accepted->triangleIndex = pair->hash;
		accepted->feature = b3_featureHullFace;
		for ( int i = 0; i < pointCount; ++i )
		{
			accepted->points[i] = manifold.points[i];
			accepted->points[i].point = b3TransformPoint( transformAtoB, b3Add( manifold.points[i].point, cornerA ) );
		}

		acceptedManifolds[acceptedManifoldCount++] = accepted;
	}

	// Terra: the same counters as the sampled path, for comparison
	context->voxelGridPairCount += 1;
	for ( int i = 0; i < acceptedManifoldCount; ++i )
	{
		context->voxelPointCount += acceptedManifolds[i]->pointCount;
		for ( int j = 0; j < acceptedManifolds[i]->pointCount; ++j )
		{
			context->voxelMinSeparation = b3MinFloat( context->voxelMinSeparation, acceptedManifolds[i]->points[j].separation );
		}
	}

	// Save the cache, still ascending by hash
	b3Array_Resize( meshContact->triangleCache, pairCount );
	for ( int i = 0; i < pairCount; ++i )
	{
		meshContact->triangleCache.data[i] = (b3TriangleCache){ pairs[i].hash, caches[i] };
	}

	if ( b3BuildClusterManifolds( world, contact, acceptedManifolds, acceptedManifoldCount, xfA, xfB, 0.0f, arena ) == false )
	{
		return false;
	}

	// One material for each whole grid
	const b3SurfaceMaterial* materialA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );

	contact->friction =
		world->frictionCallback( materialA->friction, materialA->userMaterialId, materialB->friction, materialB->userMaterialId );
	contact->restitution = world->restitutionCallback( materialA->restitution, materialA->userMaterialId, materialB->restitution,
													   materialB->userMaterialId );
	contact->rollingResistance = 0.0f;

	b3Vec3 tangentVelocityA = b3RotateVector( xfA.q, materialA->tangentVelocity );
	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );
	return true;
}

// Mesh or height field against voxel grid
//
// The mesh is A and the voxel grid B, the order the contact table gives them. Each box of the grid that can reach the mesh is a
// convex shape B in the frame of its cell: the triangles near it are collided with it through the same code as any convex
// shape against a mesh, ghost collision handling included. A contact that would push a box through a face that solid covers is
// answered on its exposed faces. The manifolds of all the boxes are then clustered together, so a shard resting on terrain is
// a few manifolds of at most four points. A pair of triangle and box is keyed by a hash of both indices, ascending in the
// contact cache.

// The triangles one box can see, and the pairs of triangle and box in one contact. Past these the rest is left out.
#define B3_MAX_BOX_TRIANGLES 128
#define B3_MAX_BOX_POINTS_PER_TRIANGLE 8

static int b3FindVoxelPairCache( const b3TriangleCache* sorted, int count, int hash )
{
	int low = 0, high = count - 1;
	while ( low <= high )
	{
		int middle = ( low + high ) / 2;
		int value = sorted[middle].triangleIndex;
		if ( value == hash )
		{
			return middle;
		}

		if ( value < hash )
		{
			low = middle + 1;
		}
		else
		{
			high = middle - 1;
		}
	}

	return B3_NULL_INDEX;
}

bool b3ComputeMeshVoxelGridManifolds( b3World* world, int workerIndex, b3Contact* contact, const b3Shape* shapeA,
											 b3WorldTransform xfA, const b3Shape* shapeB, b3WorldTransform xfB, bool isFast,
											 b3Arena arena )
{
	B3_ASSERT( shapeA->type == b3_meshShape || shapeA->type == b3_heightShape );
	B3_ASSERT( shapeB->type == b3_voxelGridShape );

	const b3VoxelGrid* grid = shapeB->voxelGrid;
	b3TaskContext* context = b3Array_Get( world->taskContexts, workerIndex );
	b3MeshContact* meshContact = &contact->meshContact;

	// The mesh frame into the body frame of the grid, and back
	b3Transform transformAtoB = b3InvMulWorldTransforms( xfB, xfA );
	b3Transform transformBtoA = b3InvMulWorldTransforms( xfA, xfB );

	// Keep the old cache, the new one is built in the contact array
	int oldCount = meshContact->triangleCache.count;
	b3TriangleCache* oldCache = b3Bump( &arena, oldCount * sizeof( b3TriangleCache ) );
	if ( oldCount > 0 )
	{
		memcpy( oldCache, meshContact->triangleCache.data, oldCount * sizeof( b3TriangleCache ) );
	}
	b3Array_Resize( meshContact->triangleCache, 0 );

	// The boxes of the grid in the region where the two can touch, with their bounds in the mesh frame
	float margin = 2.0f * B3_SPECULATIVE_DISTANCE;
	b3AABB region = { b3Max( shapeA->aabb.lowerBound, shapeB->aabb.lowerBound ),
					  b3Min( shapeA->aabb.upperBound, shapeB->aabb.upperBound ) };
	bool overlap = region.lowerBound.x <= region.upperBound.x && region.lowerBound.y <= region.upperBound.y &&
				   region.lowerBound.z <= region.upperBound.z;

	b3VoxelSideBox* boxes = NULL;
	b3VoxelSideContext side = { grid, transformBtoA, NULL, 0 };
	if ( overlap )
	{
		region.lowerBound = b3Sub( region.lowerBound, ( b3Vec3 ){ margin, margin, margin } );
		region.upperBound = b3Add( region.upperBound, ( b3Vec3 ){ margin, margin, margin } );

		b3Transform gridTransform = b3ToRelativeTransform( xfB, b3Pos_zero );
		b3AABB gridRegion = b3AABB_Transform( b3InvertTransform( gridTransform ), region );

		b3QueryVoxelGrid( grid, gridRegion, b3CollectVoxelSideCallback, &side );
		boxes = b3Bump( &arena, side.count * sizeof( b3VoxelSideBox ) );
		side.boxes = boxes;
		side.count = 0;
		b3QueryVoxelGrid( grid, gridRegion, b3CollectVoxelSideCallback, &side );
	}

	// Buffers for the boxes, reused from one box to the next
	b3MeshConvexInput input = { 0 };
	input.context = context;
	input.meshShape = shapeA;
	input.typeB = b3_hullShape;
	input.isFast = isFast;
	input.enableSpeculative = contact->flags & b3_enableSpeculativePoints;

	int* triangleIndices = b3Bump( &arena, B3_MAX_BOX_TRIANGLES * sizeof( int ) );
	b3TriangleCache* work = b3Bump( &arena, B3_MAX_BOX_TRIANGLES * sizeof( b3TriangleCache ) );
	int pointBufferCapacity = B3_MAX_BOX_POINTS_PER_TRIANGLE * B3_MAX_BOX_TRIANGLES;
	b3LocalManifoldPoint* pointBuffer = b3Bump( &arena, pointBufferCapacity * sizeof( b3LocalManifoldPoint ) );
	b3LocalManifold* manifoldBuffer = b3Bump( &arena, B3_MAX_BOX_TRIANGLES * sizeof( b3LocalManifold ) );
	b3LocalManifold** boxAccepted = b3Bump( &arena, B3_MAX_BOX_TRIANGLES * sizeof( b3LocalManifold* ) );

	b3LocalManifold** acceptedManifolds = b3Bump( &arena, B3_MAX_VOXEL_PAIRS * sizeof( b3LocalManifold* ) );
	int acceptedManifoldCount = 0;
	int pairCount = 0;

	// Materials are averaged over the accepted points, as for any mesh contact
	const b3SurfaceMaterial* materialsA = b3GetShapeMaterials( shapeA );
	const b3SurfaceMaterial* materialB = b3GetShapeMaterials( shapeB );
	const uint8_t* materialIndices = NULL;
	if ( shapeA->materialCount > 0 )
	{
		materialIndices = shapeA->type == b3_meshShape ? b3GetMeshMaterialIndices( shapeA->mesh.data )
													   : b3GetHeightFieldMaterialIndices( shapeA->heightField );
	}
	float friction = 0.0f, restitution = 0.0f, sampleCount = 0.0f;
	b3Vec3 tangentVelocitySum = b3Vec3_zero;

	bool truncated = false;
	for ( int boxIndex = 0; boxIndex < side.count && truncated == false; ++boxIndex )
	{
		const b3VoxelSideBox* box = boxes + boxIndex;
		int cell = b3VoxelGridKeyCell( grid, box->key ), boxInCell = b3VoxelGridKeyBox( grid, box->key );
		const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
		B3_ASSERT( module != NULL && boxInCell < module->boxCount );

		// The triangles near the box
		b3AABB query = box->bounds;
		query.lowerBound = b3Sub( query.lowerBound, ( b3Vec3 ){ margin, margin, margin } );
		query.upperBound = b3Add( query.upperBound, ( b3Vec3 ){ margin, margin, margin } );

		int triangleCount;
		if ( shapeA->type == b3_meshShape )
		{
			triangleCount = b3QueryMeshTriangles( triangleIndices, B3_MAX_BOX_TRIANGLES, &shapeA->mesh, query );
		}
		else
		{
			triangleCount = b3QueryHeightFieldTriangles( triangleIndices, B3_MAX_BOX_TRIANGLES, shapeA->heightField, query );
		}

		if ( triangleCount == 0 )
		{
			continue;
		}

		if ( triangleCount == B3_MAX_BOX_TRIANGLES )
		{
			static bool s_once = false;
			if ( s_once == false )
			{
				b3Log( "WARNING: complex mesh detected, triangle buffer capacity of %d reached for a voxel box", B3_MAX_BOX_TRIANGLES );
				s_once = true;
			}
		}

		if ( pairCount + triangleCount > B3_MAX_VOXEL_PAIRS || acceptedManifoldCount + triangleCount > B3_MAX_VOXEL_PAIRS )
		{
			// Too many pairs for one contact
			truncated = true;
			break;
		}

		b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );
		int boxKey = box->key;
		for ( int i = 0; i < triangleCount; ++i )
		{
			int hash = b3VoxelPairHash( triangleIndices[i], boxKey );
			work[i].triangleIndex = triangleIndices[i];
			int found = b3FindVoxelPairCache( oldCache, oldCount, hash );
			work[i].cache = found != B3_NULL_INDEX ? oldCache[found].cache : (b3ContactCache){ 0 };
		}

		// The hull of the box is in the frame of its cell
		input.hullB = &module->hulls[boxInCell].base;
		input.transformAtoB = transformAtoB;
		input.transformAtoB.p = b3Sub( transformAtoB.p, corner );
		input.coveredB = box->covered;

		int count = b3CollideMeshTriangles( &input, work, triangleCount, boxAccepted, manifoldBuffer, pointBuffer, pointBufferCapacity,
											arena );

		for ( int i = 0; i < triangleCount; ++i )
		{
			b3TriangleCache entry = { b3VoxelPairHash( work[i].triangleIndex, boxKey ), work[i].cache };
			b3Array_Push( meshContact->triangleCache, entry );
		}
		pairCount += triangleCount;

		for ( int i = 0; i < count; ++i )
		{
			const b3LocalManifold* source = boxAccepted[i];
			int triangleIndex = source->triangleIndex;

			// Copy out, in the body frame of the grid, with the pair as the id of the manifold
			b3LocalManifold* accepted = b3Bump( &arena, sizeof( b3LocalManifold ) );
			*accepted = *source;
			accepted->points = b3Bump( &arena, source->pointCount * sizeof( b3LocalManifoldPoint ) );
			for ( int j = 0; j < source->pointCount; ++j )
			{
				accepted->points[j] = source->points[j];
				accepted->points[j].point = b3Add( source->points[j].point, corner );
			}
			accepted->triangleIndex = b3VoxelPairHash( triangleIndex, boxKey );
			acceptedManifolds[acceptedManifoldCount++] = accepted;

			int materialIndex = 0;
			if ( materialIndices != NULL )
			{
				materialIndex = shapeA->type == b3_meshShape ? materialIndices[triangleIndex] : materialIndices[triangleIndex >> 1];
				materialIndex = b3ClampInt( materialIndex, 0, shapeA->materialCount - 1 );
			}

			b3SurfaceMaterial material = materialsA[materialIndex];
			for ( int j = 0; j < source->pointCount; ++j )
			{
				friction += world->frictionCallback( material.friction, material.userMaterialId, materialB->friction,
													 materialB->userMaterialId );
				restitution += world->restitutionCallback( material.restitution, material.userMaterialId, materialB->restitution,
														   materialB->userMaterialId );
				tangentVelocitySum = b3Add( tangentVelocitySum, material.tangentVelocity );
				sampleCount += 1.0f;
			}
		}
	}

	if ( truncated )
	{
		static bool s_once = false;
		if ( s_once == false )
		{
			b3Log( "WARNING: dense voxel grid on a mesh, pair capacity of %d reached", B3_MAX_VOXEL_PAIRS );
			s_once = true;
		}
	}

	// The cache is ascending by hash
	{
		b3TriangleCache* data = meshContact->triangleCache.data;
		int n = meshContact->triangleCache.count;
#define LESS( i, j ) data[(int)i].triangleIndex < data[(int)j].triangleIndex
#define SWAP( i, j )                                                                                                             \
	do                                                                                                                           \
	{                                                                                                                            \
		b3TriangleCache tmp = data[(int)i];                                                                                      \
		data[(int)i] = data[(int)j];                                                                                             \
		data[(int)j] = tmp;                                                                                                      \
	}                                                                                                                            \
	while ( 0 )
		QSORT( n, LESS, SWAP );
#undef LESS
#undef SWAP
	}

	// A grid has no rest offset: it rests on a mesh as it rests on a plain box
	if ( b3BuildClusterManifolds( world, contact, acceptedManifolds, acceptedManifoldCount, xfA, xfB, 0.0f, arena ) == false )
	{
		return false;
	}

	if ( sampleCount > 0.0f )
	{
		float invCount = 1.0f / sampleCount;
		contact->friction = invCount * friction;
		contact->restitution = invCount * restitution;
	}
	else
	{
		contact->friction = world->frictionCallback( materialsA[0].friction, materialsA[0].userMaterialId, materialB->friction,
													 materialB->userMaterialId );
		contact->restitution = world->restitutionCallback( materialsA[0].restitution, materialsA[0].userMaterialId,
														   materialB->restitution, materialB->userMaterialId );
	}

	b3Vec3 tangentVelocityA = sampleCount > 0.0f ? b3MulSV( 1.0f / sampleCount, tangentVelocitySum ) : materialsA[0].tangentVelocity;
	tangentVelocityA = b3RotateVector( xfA.q, tangentVelocityA );
	b3Vec3 tangentVelocityB = b3RotateVector( xfB.q, materialB->tangentVelocity );
	contact->tangentVelocity = b3Sub( tangentVelocityA, tangentVelocityB );
	contact->rollingResistance = 0.0f;
	return true;
}
