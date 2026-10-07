// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

#include "test_macros.h"

#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/constants.h"
#include "box3d/math_functions.h"

#include <math.h>
#include <stdio.h>

// A 20 x 20 m floor of two triangles at y = 0. Its normal is up when facingUp, down when not.
static b3MeshData* MakeFloorMesh( bool facingUp, bool doubleSided )
{
	b3Vec3 vertices[4] = { { -10.0f, 0.0f, -10.0f }, { 10.0f, 0.0f, -10.0f }, { 10.0f, 0.0f, 10.0f }, { -10.0f, 0.0f, 10.0f } };
	int32_t up[6] = { 0, 2, 1, 0, 3, 2 };
	int32_t down[6] = { 0, 1, 2, 0, 2, 3 };

	b3MeshDef def = { 0 };
	def.vertices = vertices;
	def.indices = facingUp ? up : down;
	def.vertexCount = 4;
	def.triangleCount = 2;
	def.doubleSided = doubleSided;
	return b3CreateMesh( &def, NULL, 0 );
}

typedef struct Floor
{
	b3WorldId worldId;
	b3MeshData* mesh;
} Floor;

static bool MakeFloor( Floor* floor, bool facingUp, bool doubleSided )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	floor->worldId = b3CreateWorld( &worldDef );
	floor->mesh = MakeFloorMesh( facingUp, doubleSided );
	if ( floor->mesh == NULL )
	{
		return false;
	}

	b3BodyDef bodyDef = b3DefaultBodyDef();
	b3BodyId bodyId = b3CreateBody( floor->worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	b3CreateMeshShape( bodyId, &shapeDef, floor->mesh, b3Vec3_one );
	return true;
}

static void DestroyFloor( Floor* floor )
{
	b3DestroyWorld( floor->worldId );
	b3DestroyMesh( floor->mesh );
}

typedef struct HitResult
{
	bool hit;
	float fraction;
	b3Vec3 normal;
} HitResult;

static float CastFcn( b3ShapeId shapeId, b3Pos point, b3Vec3 normal, float fraction, uint64_t userMaterialId, int triangleIndex,
					  int childIndex, void* context )
{
	MAYBE_UNUSED( shapeId );
	MAYBE_UNUSED( point );
	MAYBE_UNUSED( userMaterialId );
	MAYBE_UNUSED( triangleIndex );
	MAYBE_UNUSED( childIndex );
	HitResult* result = context;
	result->hit = true;
	result->fraction = fraction;
	result->normal = normal;
	return fraction;
}

static HitResult Ray( b3WorldId worldId, b3Vec3 from, b3Vec3 to )
{
	HitResult result = { 0 };
	b3World_CastRay( worldId, from, b3Sub( to, from ), b3DefaultQueryFilter(), CastFcn, &result );
	return result;
}

// The ray hits a double sided mesh from either side, with the normal turned toward the ray; a one sided mesh is seen from its
// front only.
static int DoubleSidedRay( void )
{
	for ( int variant = 0; variant < 4; ++variant )
	{
		bool facingUp = ( variant & 1 ) == 0;
		bool doubleSided = ( variant & 2 ) != 0;

		Floor floor;
		ENSURE( MakeFloor( &floor, facingUp, doubleSided ) );

		HitResult above = Ray( floor.worldId, ( b3Vec3 ){ 1.0f, 5.0f, 1.0f }, ( b3Vec3 ){ 1.0f, -5.0f, 1.0f } );
		HitResult below = Ray( floor.worldId, ( b3Vec3 ){ 1.0f, -5.0f, 1.0f }, ( b3Vec3 ){ 1.0f, 5.0f, 1.0f } );

		ENSURE( above.hit == ( doubleSided || facingUp ) );
		ENSURE( below.hit == ( doubleSided || facingUp == false ) );
		if ( above.hit )
		{
			ENSURE_SMALL( above.fraction - 0.5f, 1e-4f );
			ENSURE_SMALL( above.normal.y - 1.0f, 1e-4f );
		}
		if ( below.hit )
		{
			ENSURE_SMALL( below.fraction - 0.5f, 1e-4f );
			ENSURE_SMALL( below.normal.y + 1.0f, 1e-4f );
		}

		DestroyFloor( &floor );
	}
	return 0;
}

typedef enum Dropped
{
	DROPPED_BOX,
	DROPPED_CAPSULE,
	DROPPED_SPHERE,
} Dropped;

// Drop a body onto the floor from above, or from below (the floor seen from its back, when it faces up), and report where it
// ends up.
static int Drop( Dropped dropped, bool facingUp, bool doubleSided, bool fromAbove, float* finalY )
{
	Floor floor;
	ENSURE( MakeFloor( &floor, facingUp, doubleSided ) );

	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = ( b3Pos ){ 0.3f, fromAbove ? 2.0f : -2.0f, 0.2f };
	b3BodyId bodyId = b3CreateBody( floor.worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();

	b3BoxHull box = b3MakeBoxHull( 0.5f, 0.5f, 0.5f );
	switch ( dropped )
	{
		case DROPPED_BOX:
			b3CreateHullShape( bodyId, &shapeDef, &box.base );
			break;

		case DROPPED_CAPSULE:
			b3CreateCapsuleShape( bodyId, &shapeDef, &( b3Capsule ){ { 0.0f, -0.3f, 0.0f }, { 0.0f, 0.3f, 0.0f }, 0.4f } );
			break;

		case DROPPED_SPHERE:
			b3CreateSphereShape( bodyId, &shapeDef, &( b3Sphere ){ b3Vec3_zero, 0.5f } );
			break;
	}

	// Dropped from below, the body has to rise to meet the floor: give it a push up
	if ( fromAbove == false )
	{
		b3Body_SetLinearVelocity( bodyId, ( b3Vec3 ){ 0.0f, 8.0f, 0.0f } );
		b3Body_SetGravityScale( bodyId, 0.0f );
	}

	for ( int i = 0; i < 300; ++i )
	{
		b3World_Step( floor.worldId, 1.0f / 60.0f, 4 );
	}

	*finalY = b3Body_GetPosition( bodyId ).y;
	DestroyFloor( &floor );
	return 0;
}

// A body lands on a double sided floor from either side and stops there. A one sided floor stops a body that arrives at its
// front and lets one that arrives at its back through.
static int DoubleSidedLanding( void )
{
	const Dropped kinds[3] = { DROPPED_BOX, DROPPED_CAPSULE, DROPPED_SPHERE };
	const float radii[3] = { 0.5f, 0.7f, 0.5f };

	for ( int kind = 0; kind < 3; ++kind )
	{
		for ( int variant = 0; variant < 8; ++variant )
		{
			bool facingUp = ( variant & 1 ) == 0;
			bool doubleSided = ( variant & 2 ) != 0;
			bool fromAbove = ( variant & 4 ) != 0;

			float y = 0.0f;
			ENSURE( Drop( kinds[kind], facingUp, doubleSided, fromAbove, &y ) == 0 );

			bool held = doubleSided || ( fromAbove == facingUp );
			if ( held && fromAbove )
			{
				// Rests on the floor
				ENSURE_SMALL( y - radii[kind], 0.03f );
			}
			else if ( held )
			{
				// Thrown up at the underside of the floor, it stays below
				ENSURE( y < 0.0f );
			}
			else if ( kind == 0 )
			{
				// The rounded shapes cull the back of a one sided triangle by their center, so only the box is asserted
				ENSURE( fromAbove ? y < -5.0f : y > 5.0f );
			}
		}
	}
	return 0;
}

// A character capsule: a core segment 1.1 m long with the radius 0.35 m and a 1 mm margin
static b3ShapeProxy CharacterProxy( b3Vec3 center, b3Vec3 points[2] )
{
	points[0] = b3Sub( center, ( b3Vec3 ){ 0.0f, 0.55f, 0.0f } );
	points[1] = b3Add( center, ( b3Vec3 ){ 0.0f, 0.55f, 0.0f } );
	return ( b3ShapeProxy ){ points, 2, 0.351f };
}

static int Recover( const b3MeshData* data, b3Vec3 center, b3MeshRecoverResult* results )
{
	b3Vec3 points[2];
	b3ShapeProxy proxy = CharacterProxy( center, points );
	b3Mesh mesh = { data, b3Vec3_one };
	return b3RecoverMesh( &mesh, b3Transform_identity, &proxy, results, 8 );
}

// A character sunk into a mesh is pushed out along the face it is on; the back of a one sided mesh pushes nothing
static int MeshRecovery( void )
{
	b3MeshRecoverResult results[8];

	for ( int variant = 0; variant < 4; ++variant )
	{
		bool facingUp = ( variant & 1 ) == 0;
		bool doubleSided = ( variant & 2 ) != 0;
		b3MeshData* mesh = MakeFloorMesh( facingUp, doubleSided );
		ENSURE( mesh != NULL );

		// On the floor in the middle of the diagonal, 1 mm inside the margin; and with the core through the plane
		for ( int core = 0; core < 2; ++core )
		{
			for ( int above = 0; above < 2; ++above )
			{
				float height = core == 0 ? 0.9f : 0.3f;
				b3Vec3 center = { 0.0f, above ? height : -height, 0.0f };

				int count = Recover( mesh, center, results );
				bool pushes = doubleSided || ( above == facingUp );
				ENSURE( ( count > 0 ) == pushes );
				if ( pushes == false )
				{
					continue;
				}

				float expectedDepth = core == 0 ? 0.351f - 0.35f : 0.351f + 0.25f;
				float expectedY = above ? 1.0f : -1.0f;
				for ( int i = 0; i < count; ++i )
				{
					ENSURE_SMALL( results[i].normal.y - expectedY, 1e-4f );
					ENSURE_SMALL( results[i].depth - expectedDepth, 1e-3f );
				}

				// Moving by the push leaves the mesh
				b3Vec3 moved = b3MulAdd( center, results[0].depth, results[0].normal );
				ENSURE( Recover( mesh, moved, results ) == 0 );
			}
		}

		b3DestroyMesh( mesh );
	}
	return 0;
}

static b3Vec3 CapsuleCastNormal( b3WorldId worldId, b3Vec3 center, b3Vec3 translation, float* fraction )
{
	b3Vec3 points[2];
	b3ShapeProxy proxy = CharacterProxy( center, points );
	HitResult result = { 0 };
	b3World_CastShape( worldId, b3Vec3_zero, &proxy, translation, b3DefaultQueryFilter(), CastFcn, &result );
	*fraction = result.hit ? result.fraction : 1.0f;
	return result.normal;
}

// The character motion of the plugin in miniature: recover from the mesh, cast the margin inflated capsule, stop short of the
// hit, slide the rest of the motion along the hit normal. A hit that starts inside the margin has no normal and is left to the
// recovery of the next call, so without mesh recovery the capsule falls through.
static int Walk( bool facingUp, bool doubleSided, bool useRecovery, float* finalY )
{
	Floor floor;
	ENSURE( MakeFloor( &floor, facingUp, doubleSided ) );
	b3Mesh mesh = { floor.mesh, b3Vec3_one };

	b3Vec3 position = { 0.0f, 1.5f, 0.0f };
	b3Vec3 velocity = { 2.0f, 0.0f, 0.0f };
	for ( int tick = 0; tick < 90; ++tick )
	{
		velocity.y -= 9.8f / 60.0f;
		b3Vec3 motion = b3MulSV( 1.0f / 60.0f, velocity );
		for ( int slide = 0; slide < 4 && b3LengthSquared( motion ) > 1e-10f; ++slide )
		{
			// Recovery first, as the plugin does, pushing out along each axis by the deepest triangle
			for ( int attempt = 0; attempt < 4 && useRecovery; ++attempt )
			{
				b3Vec3 points[2];
				b3ShapeProxy proxy = CharacterProxy( position, points );
				b3MeshRecoverResult results[8];
				int count = b3RecoverMesh( &mesh, b3Transform_identity, &proxy, results, 8 );
				if ( count == 0 )
				{
					break;
				}
				b3Vec3 push = b3Vec3_zero;
				for ( int i = 0; i < count; ++i )
				{
					b3Vec3 candidate = b3MulSV( results[i].depth + B3_LINEAR_SLOP, results[i].normal );
					push.x = fabsf( candidate.x ) > fabsf( push.x ) ? candidate.x : push.x;
					push.y = fabsf( candidate.y ) > fabsf( push.y ) ? candidate.y : push.y;
					push.z = fabsf( candidate.z ) > fabsf( push.z ) ? candidate.z : push.z;
				}
				position = b3Add( position, push );
			}

			float fraction;
			b3Vec3 normal = CapsuleCastNormal( floor.worldId, position, motion, &fraction );
			if ( fraction >= 1.0f || b3LengthSquared( normal ) < 0.25f )
			{
				position = b3Add( position, motion );
				break;
			}
			float length = b3Length( motion );
			float safe = b3MaxFloat( 0.0f, fraction - B3_LINEAR_SLOP / length );
			position = b3MulAdd( position, safe, motion );
			b3Vec3 remainder = b3MulSV( 1.0f - safe, motion );
			motion = b3Sub( remainder, b3MulSV( b3Dot( remainder, normal ), normal ) );
			velocity = b3Sub( velocity, b3MulSV( b3Dot( velocity, normal ), normal ) );
		}
	}

	*finalY = position.y;
	DestroyFloor( &floor );
	return 0;
}

static int CharacterWalksOnMesh( void )
{
	// Without recovery the capsule sinks into the margin and falls through: the failure this fixes
	float y = 0.0f;
	ENSURE( Walk( true, false, false, &y ) == 0 );
	ENSURE( y < -1.0f );

	for ( int variant = 0; variant < 4; ++variant )
	{
		bool facingUp = ( variant & 1 ) == 0;
		bool doubleSided = ( variant & 2 ) != 0;

		// A one sided floor facing down is no floor for a capsule above it
		if ( doubleSided == false && facingUp == false )
		{
			continue;
		}

		ENSURE( Walk( facingUp, doubleSided, true, &y ) == 0 );
		ENSURE( y > 0.88f && y < 0.92f );
	}
	return 0;
}

int MeshDoubleSidedTest( void )
{
	RUN_SUBTEST( DoubleSidedRay );
	RUN_SUBTEST( DoubleSidedLanding );
	RUN_SUBTEST( MeshRecovery );
	RUN_SUBTEST( CharacterWalksOnMesh );
	return 0;
}
