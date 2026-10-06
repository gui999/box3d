// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

#include "test_macros.h"

// Voxel grid internals are visible to the tests
#include "manifold.h"
#include "voxel_grid.h"

#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/math_functions.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CELL_METERS 2.0f
#define CELL_VOXELS 20
#define SLAB_HEIGHT 0.4f

// A flat floor of N x N cells, each cell holding the same full slab module. Optional pillar module in one cell.
typedef struct Floor
{
	b3VoxelGridModule* slab;
	b3VoxelGridModule* pillar;
	b3VoxelGrid* grid;
	int size;
} Floor;

static int PaddedIndex( const Floor* floor, int i, int j, int k )
{
	return ( ( j + 1 ) * ( floor->size + 2 ) + ( k + 1 ) ) * ( floor->size + 2 ) + ( i + 1 );
}

static bool MakeFloor( Floor* floor, int size, int emptyI, int emptyK, int pillarI )
{
	floor->size = size;

	float slabBox[6] = { 0.0f, 0.0f, 0.0f, CELL_METERS, SLAB_HEIGHT, CELL_METERS };
	floor->slab = b3CreateVoxelGridModule( slabBox, 1, CELL_METERS, CELL_VOXELS );

	// A full cell of solid, 2 m tall
	float pillarBox[6] = { 0.0f, 0.0f, 0.0f, CELL_METERS, CELL_METERS, CELL_METERS };
	floor->pillar = b3CreateVoxelGridModule( pillarBox, 1, CELL_METERS, CELL_VOXELS );
	if ( floor->slab == NULL || floor->pillar == NULL )
	{
		return false;
	}

	int paddedCount = ( size + 2 ) * 3 * ( size + 2 );
	int* padded = (int*)b3Alloc( paddedCount * sizeof( int ) );
	for ( int i = 0; i < paddedCount; ++i )
	{
		padded[i] = -1;
	}

	b3VoxelGridModule* modules[2] = { floor->slab, floor->pillar };

	b3VoxelGridDef def = { 0 };
	def.cellCountX = size;
	def.cellCountY = 1;
	def.cellCountZ = size;
	def.origin = b3Vec3_zero;
	def.cellMeters = CELL_METERS;
	def.cellVoxels = CELL_VOXELS;
	def.maxBoxesPerCell = 256;
	def.modules = modules;
	def.moduleCount = 2;

	// The padded layout has Y size 1 plus two ring layers
	for ( int k = -1; k <= size; ++k )
	{
		for ( int i = -1; i <= size; ++i )
		{
			// Ring cells carry the neighbouring module too so a floor edge is not seen as a covered seam to nothing
			int index = ( ( 0 + 1 ) * ( size + 2 ) + ( k + 1 ) ) * ( size + 2 ) + ( i + 1 );
			bool inside = 0 <= i && i < size && 0 <= k && k < size;
			if ( inside && ( i != emptyI || k != emptyK ) )
			{
				padded[index] = ( i == pillarI ) ? 1 : 0;
			}
		}
	}
	def.paddedCells = padded;

	floor->grid = b3CreateVoxelGrid( &def );
	b3Free( padded, paddedCount * sizeof( int ) );
	return floor->grid != NULL;
}

static void DestroyFloor( Floor* floor )
{
	b3ReleaseVoxelGrid( floor->grid );
	b3ReleaseVoxelGridModule( floor->slab );
	b3ReleaseVoxelGridModule( floor->pillar );
}

static int VoxelGridModuleCover( void )
{
	// A slab that spans a cell has its +-X, +-Z faces on the cell boundary: their cover depends on the neighbour
	Floor floor;
	ENSURE( MakeFloor( &floor, 4, -1, -1, -1 ) );

	const b3VoxelGridModule* slab = floor.slab;
	ENSURE( slab->boxCount == 1 );
	ENSURE( slab->boundaryLayers != NULL );
	// Top and bottom faces are on the cell boundary too? The bottom is, the top is inside the cell and open air above
	ENSURE( ( slab->boundaryFaces[0] & ( 1 << 0 ) ) != 0 );
	ENSURE( ( slab->boundaryFaces[0] & ( 1 << 2 ) ) != 0 );
	ENSURE( ( slab->boundaryFaces[0] & ( 1 << 3 ) ) == 0 );

	// Inner cells: all four side faces are covered by the neighbour slabs, the top face is exposed
	int cell = b3VoxelGridCellIndex( floor.grid, 1, 0, 1 );
	uint8_t covered = b3GetVoxelBoxCoveredFaces( floor.grid, cell, 0 );
	ENSURE( covered == ( ( 1 << 0 ) | ( 1 << 1 ) | ( 1 << 4 ) | ( 1 << 5 ) ) );

	// A corner cell keeps its outer faces exposed
	cell = b3VoxelGridCellIndex( floor.grid, 0, 0, 0 );
	covered = b3GetVoxelBoxCoveredFaces( floor.grid, cell, 0 );
	ENSURE( covered == ( ( 1 << 1 ) | ( 1 << 5 ) ) );

	DestroyFloor( &floor );
	return 0;
}

static int VoxelGridRayCast( void )
{
	Floor floor;
	// Cell (2, 2) is empty
	ENSURE( MakeFloor( &floor, 4, 2, 2, -1 ) );

	// Straight down onto the top of the slab in cell (0, 0)
	{
		b3RayCastInput input = { { 1.0f, 2.0f, 1.0f }, { 0.0f, -3.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE_SMALL( output.fraction - ( 2.0f - SLAB_HEIGHT ) / 3.0f, 1.0e-5f );
		ENSURE_SMALL( output.normal.x, 1.0e-6f );
		ENSURE_SMALL( output.normal.y - 1.0f, 1.0e-6f );
		ENSURE_SMALL( output.normal.z, 1.0e-6f );
		ENSURE_SMALL( output.point.y - SLAB_HEIGHT, 1.0e-5f );
		ENSURE( output.triangleIndex == b3VoxelGridCellIndex( floor.grid, 0, 0, 0 ) );
		ENSURE( output.childIndex == 0 );
	}

	// Into the cell (3, 1), a different cell reports its own index
	{
		b3RayCastInput input = { { 7.0f, 2.0f, 3.0f }, { 0.0f, -3.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE( output.triangleIndex == b3VoxelGridCellIndex( floor.grid, 3, 0, 1 ) );
	}

	// Through the empty cell: a miss
	{
		b3RayCastInput input = { { 5.0f, 2.0f, 5.0f }, { 0.0f, -3.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit == false );
	}

	// From the side, the entered face is the normal. Crossing several cells reaches the first box.
	{
		b3RayCastInput input = { { -1.0f, 0.2f, 1.0f }, { 3.0f, 0.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE_SMALL( output.fraction - 1.0f / 3.0f, 1.0e-5f );
		ENSURE_SMALL( output.normal.x + 1.0f, 1.0e-6f );
	}

	// Along the seam between two cells, on the slab surface plane: hits the box face in the walked column
	{
		b3RayCastInput input = { { 3.9f, 1.0f, 1.0f }, { 0.0f, -1.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE_SMALL( output.point.y - SLAB_HEIGHT, 1.0e-5f );
	}

	// Starting inside a box: a hit at the origin with no normal
	{
		b3RayCastInput input = { { 1.0f, 0.2f, 1.0f }, { 0.0f, -1.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE( output.fraction == 0.0f );
		ENSURE( output.normal.x == 0.0f && output.normal.y == 0.0f && output.normal.z == 0.0f );
	}

	// Pointing away from the grid
	{
		b3RayCastInput input = { { 1.0f, 2.0f, 1.0f }, { 0.0f, 3.0f, 0.0f }, 1.0f };
		b3CastOutput output = b3RayCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit == false );
	}

	DestroyFloor( &floor );
	return 0;
}

static b3Vec3 BoxPoints[8];

static b3ShapeProxy MakeBoxProxy( b3Vec3 center, float half )
{
	int n = 0;
	for ( int x = -1; x <= 1; x += 2 )
	{
		for ( int y = -1; y <= 1; y += 2 )
		{
			for ( int z = -1; z <= 1; z += 2 )
			{
				BoxPoints[n++] = (b3Vec3){ center.x + (float)x * half, center.y + (float)y * half, center.z + (float)z * half };
			}
		}
	}

	return (b3ShapeProxy){ BoxPoints, 8, 0.0f };
}

static int VoxelGridShapeCast( void )
{
	Floor floor;
	// The third column holds a 2 m pillar
	ENSURE( MakeFloor( &floor, 4, -1, -1, 2 ) );

	// A box just above the floor sliding across the seams between the slabs hits nothing
	{
		b3ShapeCastInput input = { 0 };
		input.proxy = MakeBoxProxy( ( b3Vec3 ){ 0.5f, SLAB_HEIGHT + 0.25f + 0.05f, 1.0f }, 0.25f );
		input.translation = ( b3Vec3 ){ 3.0f, 0.0f, 0.0f };
		input.maxFraction = 1.0f;
		b3CastOutput output = b3ShapeCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit == false );
	}

	// The same box moving down hits the floor
	{
		b3ShapeCastInput input = { 0 };
		input.proxy = MakeBoxProxy( ( b3Vec3 ){ 2.0f, SLAB_HEIGHT + 0.25f + 0.2f, 2.0f }, 0.25f );
		input.translation = ( b3Vec3 ){ 0.0f, -1.0f, 0.0f };
		input.maxFraction = 1.0f;
		b3CastOutput output = b3ShapeCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE( output.fraction > 0.0f && output.fraction < 0.5f );
		ENSURE_SMALL( output.normal.y - 1.0f, 1.0e-4f );
	}

	// Moving toward the pillar at x = 4 above the slab, a real wall, is a hit on its exposed face
	{
		b3ShapeCastInput input = { 0 };
		input.proxy = MakeBoxProxy( ( b3Vec3 ){ 2.0f, 1.0f, 1.0f }, 0.25f );
		input.translation = ( b3Vec3 ){ 4.0f, 0.0f, 0.0f };
		input.maxFraction = 1.0f;
		b3CastOutput output = b3ShapeCastVoxelGrid( floor.grid, &input );
		ENSURE( output.hit );
		ENSURE_SMALL( output.normal.x + 1.0f, 1.0e-4f );
		ENSURE( output.point.x > 3.9f && output.point.x < 4.1f );
	}

	DestroyFloor( &floor );
	return 0;
}

typedef struct World
{
	b3WorldId worldId;
	b3BodyId floorBodyId;
	b3ShapeId floorShapeId;
	Floor floor;
} World;

static void CreateWorld( World* world, int size, float friction )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	world->worldId = b3CreateWorld( &worldDef );

	MakeFloor( &world->floor, size, -1, -1, -1 );

	b3BodyDef bodyDef = b3DefaultBodyDef();
	world->floorBodyId = b3CreateBody( world->worldId, &bodyDef );

	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = friction;
	world->floorShapeId = b3CreateVoxelGridShape( world->floorBodyId, &shapeDef, world->floor.grid );
}

static void DestroyWorld( World* world )
{
	b3DestroyWorld( world->worldId );
	DestroyFloor( &world->floor );
}

static b3BodyId CreateBox( World* world, b3Vec3 position, float half, float friction )
{
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = position;
	b3BodyId bodyId = b3CreateBody( world->worldId, &bodyDef );

	b3BoxHull box = b3MakeBoxHull( half, half, half );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = friction;
	b3CreateHullShape( bodyId, &shapeDef, &box.base );
	return bodyId;
}

static int VoxelGridShapeLifetime( void )
{
	World world;
	CreateWorld( &world, 4, 0.5f );
	ENSURE( B3_IS_NON_NULL( world.floorShapeId ) );
	ENSURE( b3Shape_GetVoxelGrid( world.floorShapeId ) == world.floor.grid );
	ENSURE( world.floor.grid->refCount.value == 2 );

	// World ray cast reaches the grid through the shape dispatch
	b3RayCastInput input = { { 1.0f, 2.0f, 1.0f }, { 0.0f, -3.0f, 0.0f }, 1.0f };
	b3WorldCastOutput output = b3Shape_RayCast( world.floorShapeId, input.origin, input.translation );
	ENSURE( output.hit );
	ENSURE_SMALL( output.normal.y - 1.0f, 1.0e-5f );

	b3DestroyWorld( world.worldId );
	// The world's shape released its reference
	ENSURE( world.floor.grid->refCount.value == 1 );
	DestroyFloor( &world.floor );
	return 0;
}

static int VoxelGridRestAndSleep( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );

	// On the corner where four slabs meet
	float half = 0.5f;
	b3BodyId boxId = CreateBox( &world, ( b3Vec3 ){ 4.0f, SLAB_HEIGHT + half + 0.02f, 4.0f }, half, 0.6f );

	for ( int i = 0; i < 360; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( boxId );
	printf( "  rest: position %.4f %.4f %.4f, awake %d\n", p.x, p.y, p.z, b3Body_IsAwake( boxId ) );
	ENSURE_SMALL( p.x - 4.0f, 0.01f );
	ENSURE_SMALL( p.z - 4.0f, 0.01f );
	// No rest offset: the box rests on the grid as it rests on a plain box
	ENSURE_SMALL( p.y - ( SLAB_HEIGHT + half ), 0.002f );
	ENSURE( b3Body_IsAwake( boxId ) == false );

	DestroyWorld( &world );
	return 0;
}

// A proxy that sank into the floor is pushed out through the exposed face, never into a seam between two floor boxes
static int VoxelGridRecovery( void )
{
	// The third column holds a 2 m tall pillar wall that starts at x = 4
	Floor floor;
	ENSURE( MakeFloor( &floor, 4, -1, -1, 2 ) );
	b3Transform transform = b3Transform_identity;

	// A capsule sunk 2 cm into the slab top, centred on the seam at x = 2
	{
		b3Vec3 points[2] = { { 2.0f, SLAB_HEIGHT - 0.02f + 0.35f, 1.0f }, { 2.0f, SLAB_HEIGHT - 0.02f + 1.05f, 1.0f } };
		b3ShapeProxy proxy = { points, 2, 0.35f };
		b3Vec3 pushes[3];
		int count = b3RecoverVoxelGrid( floor.grid, transform, &proxy, pushes );
		printf( "  seam recovery: %d pushes, up %.4f, x %.4f, z %.4f\n", count, pushes[1].y, pushes[0].x, pushes[2].z );
		ENSURE( count >= 2 );
		ENSURE_SMALL( pushes[1].y - 0.02f, 1.0e-4f );
		ENSURE( pushes[0].x == 0.0f );
		ENSURE( pushes[2].z == 0.0f );
	}

	// The same capsule on the seam between four slabs at the corner (x = 2, z = 2)
	{
		b3Vec3 points[2] = { { 2.0f, SLAB_HEIGHT - 0.02f + 0.35f, 2.0f }, { 2.0f, SLAB_HEIGHT - 0.02f + 1.05f, 2.0f } };
		b3ShapeProxy proxy = { points, 2, 0.35f };
		b3Vec3 pushes[3];
		int count = b3RecoverVoxelGrid( floor.grid, transform, &proxy, pushes );
		ENSURE( count >= 4 );
		ENSURE_SMALL( pushes[1].y - 0.02f, 1.0e-4f );
		ENSURE( pushes[0].x == 0.0f );
		ENSURE( pushes[2].z == 0.0f );
	}

	// A box sunk 2 cm into the floor and 3 cm into the wall face at x = 4 leaves through both exposed faces
	{
		b3Vec3 points[8];
		MakeBoxProxy( ( b3Vec3 ){ 3.78f, SLAB_HEIGHT - 0.02f + 0.25f, 1.0f }, 0.25f );
		for ( int i = 0; i < 8; ++i )
		{
			points[i] = BoxPoints[i];
		}
		b3ShapeProxy proxy = { points, 8, 0.0f };
		b3Vec3 pushes[3];
		int count = b3RecoverVoxelGrid( floor.grid, transform, &proxy, pushes );
		printf( "  corner recovery: %d pushes, x %.4f, y %.4f, z %.4f\n", count, pushes[0].x, pushes[1].y, pushes[2].z );
		ENSURE( count >= 2 );
		ENSURE_SMALL( pushes[0].x + 0.03f, 1.0e-4f );
		ENSURE_SMALL( pushes[1].y - 0.02f, 1.0e-4f );
		ENSURE( pushes[2].z == 0.0f );
	}

	// Resting on the surface, or clear of it, nothing pushes
	{
		b3Vec3 points[2] = { { 2.0f, SLAB_HEIGHT + 0.35f, 1.0f }, { 2.0f, SLAB_HEIGHT + 1.05f, 1.0f } };
		b3ShapeProxy proxy = { points, 2, 0.35f };
		b3Vec3 pushes[3];
		ENSURE( b3RecoverVoxelGrid( floor.grid, transform, &proxy, pushes ) == 0 );
		points[0].y += 0.5f;
		points[1].y += 0.5f;
		ENSURE( b3RecoverVoxelGrid( floor.grid, transform, &proxy, pushes ) == 0 );
	}

	// A grid turned a quarter about Y: the push is answered in world space
	{
		b3Transform turned = { { 0.0f, 0.0f, 0.0f }, b3MakeQuatFromAxisAngle( b3Vec3_axisY, 0.5f * B3_PI ) };
		// The seam at grid x = 2 is the world line z = -2, the floor top stays at y = 0.4
		b3Vec3 points[2] = { { 1.0f, SLAB_HEIGHT - 0.02f + 0.35f, -2.0f }, { 1.0f, SLAB_HEIGHT - 0.02f + 1.05f, -2.0f } };
		b3ShapeProxy proxy = { points, 2, 0.35f };
		b3Vec3 pushes[3];
		int count = b3RecoverVoxelGrid( floor.grid, turned, &proxy, pushes );
		ENSURE( count >= 2 );
		ENSURE_SMALL( pushes[1].y - 0.02f, 1.0e-4f );
		ENSURE_SMALL( b3Length( pushes[0] ), 1.0e-6f );
		ENSURE_SMALL( b3Length( pushes[2] ), 1.0e-6f );
	}

	DestroyFloor( &floor );
	return 0;
}

// A box slides over the seams between slabs. Resting on a flat floor made of several boxes it must keep its direction:
// no vertical or lateral velocity spikes. The gravity is scaled up to make the contact penetrate deeper, which is when a
// ghost contact on a seam face wins the separating axis test.
static int SlideAcrossSeams( float gravity, float speed, int tickCount, float minX, float maxSpike )
{
	World world;
	CreateWorld( &world, 4, 0.01f );
	b3World_SetGravity( world.worldId, ( b3Vec3 ){ 0.0f, -gravity, 0.0f } );

	float half = 0.5f;
	b3BodyId boxId = CreateBox( &world, ( b3Vec3 ){ 1.0f, SLAB_HEIGHT + half + 0.01f, 3.5f }, half, 0.01f );

	// Settle
	for ( int i = 0; i < 30; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Body_SetLinearVelocity( boxId, ( b3Vec3 ){ speed, 0.0f, 0.0f } );

	float maxVertical = 0.0f, maxLateral = 0.0f;
	float startX = b3Body_GetPosition( boxId ).x;
	for ( int i = 0; i < tickCount; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
		b3Vec3 v = b3Body_GetLinearVelocity( boxId );
		maxVertical = b3MaxFloat( maxVertical, fabsf( v.y ) );
		maxLateral = b3MaxFloat( maxLateral, fabsf( v.z ) );
	}

	b3Pos p = b3Body_GetPosition( boxId );
	printf( "  slide g=%.0f v=%.0f: travelled %.3f m, max |vy| %.4f, max |vz| %.4f, y %.4f\n", gravity, speed,
			p.x - startX, maxVertical, maxLateral, p.y );

	// It crossed the seams (x = 2, 4, 6) without any vertical or lateral velocity spike
	ENSURE( p.x > minX );
	ENSURE( maxVertical < maxSpike );
	ENSURE( maxLateral < maxSpike );

	DestroyWorld( &world );
	return 0;
}

static int VoxelGridSlideAcrossSeams( void )
{
	return SlideAcrossSeams( 10.0f, 3.0f, 110, 6.2f, 0.05f );
}

static int VoxelGridSlideAcrossSeamsHeavy( void )
{
	return SlideAcrossSeams( 400.0f, 12.0f, 40, 6.2f, 0.05f );
}

// The normal impulse the contacts of a body report for one step: summed over every manifold point of every pair.
static float StepAndSumNormalImpulse( b3WorldId worldId, b3BodyId bodyId, float dt, int* pointCount )
{
	b3World_Step( worldId, dt, 4 );
	b3ContactData pairs[16];
	int pairCount = b3Body_GetContactData( bodyId, pairs, 16 );
	float sum = 0.0f;
	*pointCount = 0;
	for ( int i = 0; i < pairCount; ++i )
	{
		for ( int m = 0; m < pairs[i].manifoldCount; ++m )
		{
			for ( int p = 0; p < pairs[i].manifolds[m].pointCount; ++p )
			{
				sum += pairs[i].manifolds[m].points[p].appliedNormalImpulse;
				*pointCount += 1;
			}
		}
	}
	return sum;
}

static int RestingImpulse( bool onGrid )
{
	World world;
	CreateWorld( &world, 4, 0.6f );
	float surface = SLAB_HEIGHT;
	if ( onGrid == false )
	{
		// The same top plane made of one plain box hull
		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.position = ( b3Vec3 ){ 20.0f, 0.2f, 0.0f };
		b3BodyId plainId = b3CreateBody( world.worldId, &bodyDef );
		b3BoxHull slab = b3MakeBoxHull( 4.0f, 0.2f, 4.0f );
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		shapeDef.baseMaterial.friction = 0.6f;
		b3CreateHullShape( plainId, &shapeDef, &slab.base );
	}

	float half = 0.5f;
	b3Vec3 start = onGrid ? ( b3Vec3 ){ 4.0f, surface + half + 0.02f, 4.0f } : ( b3Vec3 ){ 20.0f, surface + half + 0.02f, 0.0f };
	b3BodyId boxId = CreateBox( &world, start, half, 0.6f );
	float mass = b3Body_GetMass( boxId );

	const float dt = 1.0f / 60.0f;
	int points = 0;
	float sum = 0.0f;
	for ( int i = 0; i < 30; ++i )
	{
		sum = StepAndSumNormalImpulse( world.worldId, boxId, dt, &points );
	}
	float expected = mass * 10.0f * dt;
	printf( "  resting impulse (%s): %d points, sum %.3f, expected %.3f, ratio %.3f\n", onGrid ? "grid" : "box", points, sum,
			expected, sum / expected );
	ENSURE( points > 0 );
	ENSURE( fabsf( sum - expected ) < 0.05f * expected );

	DestroyWorld( &world );
	return 0;
}

// A heavy box dropped on the floor: over the whole landing the contacts apply the momentum it arrives with plus its weight
// for the time it took, however deep the contact pushes in and out.
static int VoxelGridLandingImpulse( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );
	float half = 0.5f;
	b3BodyId boxId = CreateBox( &world, ( b3Vec3 ){ 4.0f, SLAB_HEIGHT + half + 0.1f, 4.0f }, half, 0.6f );
	b3Body_SetLinearVelocity( boxId, ( b3Vec3 ){ 0.0f, -8.0f, 0.0f } );
	float mass = b3Body_GetMass( boxId );

	const float dt = 1.0f / 60.0f;
	const int tickCount = 20;
	float applied = 0.0f, total = 0.0f;
	for ( int i = 0; i < tickCount; ++i )
	{
		b3World_Step( world.worldId, dt, 4 );
		b3ContactData pairs[16];
		int pairCount = b3Body_GetContactData( boxId, pairs, 16 );
		for ( int j = 0; j < pairCount; ++j )
		{
			for ( int m = 0; m < pairs[j].manifoldCount; ++m )
			{
				for ( int p = 0; p < pairs[j].manifolds[m].pointCount; ++p )
				{
					applied += pairs[j].manifolds[m].points[p].appliedNormalImpulse;
					total += pairs[j].manifolds[m].points[p].totalNormalImpulse;
				}
			}
		}
	}

	// The body ends at rest: all its downward momentum and the weight over the ticks went into the floor
	float expected = mass * ( 8.0f + 10.0f * dt * (float)tickCount );
	printf( "  landing impulse: applied %.0f, legacy total %.0f, expected %.0f\n", applied, total, expected );
	ENSURE( fabsf( applied - expected ) < 0.05f * expected );

	DestroyWorld( &world );
	return 0;
}

static int VoxelGridRestingImpulseGrid( void )
{
	return RestingImpulse( true );
}

static int VoxelGridRestingImpulseBox( void )
{
	return RestingImpulse( false );
}

static int VoxelGridSetCellWakes( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );

	float half = 0.5f;
	b3BodyId boxId = CreateBox( &world, ( b3Vec3 ){ 1.0f, SLAB_HEIGHT + half + 0.02f, 1.0f }, half, 0.6f );

	for ( int i = 0; i < 300; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	ENSURE( b3Body_IsAwake( boxId ) == false );
	ENSURE( b3Body_GetPosition( boxId ).y > 0.5f );

	// The cell under the box is removed, the box width fits inside the one cell
	int cells[1] = { PaddedIndex( &world.floor, 0, 0, 0 ) };
	int modules[1] = { -1 };
	b3Shape_VoxelGridSetCells( world.floorShapeId, cells, modules, 1 );

	// A cell added back later has its module appended in place too
	ENSURE( b3Body_IsAwake( boxId ) == true );

	for ( int i = 0; i < 60; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	float y = b3Body_GetPosition( boxId ).y;
	printf( "  set cell: y after 1 s %.3f\n", y );
	ENSURE( y < -1.0f );

	DestroyWorld( &world );
	return 0;
}

static int VoxelGridSetCellRestoresSupport( void )
{
	// Removing a cell next to a resting box leaves it asleep, adding a module and a cell under a falling box catches it
	World world;
	CreateWorld( &world, 4, 0.6f );

	float half = 0.25f;
	// Above cell (1, 1), starts in the air and falls onto the floor
	b3BodyId boxId = CreateBox( &world, ( b3Vec3 ){ 3.0f, 2.0f, 3.0f }, half, 0.6f );

	// Replace the cell below with a pillar module added in place: the box lands on its top at y = 2
	float pillarBox[6] = { 0.0f, 0.0f, 0.0f, CELL_METERS, 1.0f, CELL_METERS };
	b3VoxelGridModule* tall = b3CreateVoxelGridModule( pillarBox, 1, CELL_METERS, CELL_VOXELS );
	int index = b3Shape_VoxelGridAddModule( world.floorShapeId, tall );
	b3ReleaseVoxelGridModule( tall );
	ENSURE( index == 2 );

	int cells[1] = { PaddedIndex( &world.floor, 1, 0, 1 ) };
	int modules[1] = { index };
	b3Shape_VoxelGridSetCells( world.floorShapeId, cells, modules, 1 );

	for ( int i = 0; i < 240; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( boxId );
	printf( "  add module: resting y %.3f\n", p.y );
	ENSURE_SMALL( p.y - ( 1.0f + half ), 0.03f );

	DestroyWorld( &world );
	return 0;
}

// The convex routines that answer a contact on the exposed faces only. Box of half width 1 at the origin, +X face covered.
static int VoxelGridExposedFaces( void )
{
	b3BoxHull boxA = b3MakeBoxHull( 1.0f, 1.0f, 1.0f );
	int exposed = 0x3f & ~( 1 << 1 );
	b3LocalManifoldPoint points[32];
	b3LocalManifold manifold = { 0 };
	manifold.points = points;

	// A sphere centre 0.1 inside the +X face: the regular manifold pushes out through +X
	{
		b3Sphere sphere = { { 0.9f, 0.2f, 0.0f }, 0.25f };
		b3SimplexCache cache = { 0 };
		b3CollideHullAndSphere( &manifold, 32, &boxA.base, &sphere, b3Transform_identity, &cache );
		ENSURE( manifold.pointCount == 1 );
		ENSURE( manifold.normal.x > 0.99f );

		// Answered on the least penetrated exposed face, +Y
		ENSURE( b3CollideHullFacesAndSphere( &manifold, 32, &boxA.base, exposed, &sphere, b3Transform_identity ) );
		ENSURE( manifold.pointCount == 1 );
		ENSURE_SMALL( manifold.normal.y - 1.0f, 1.0e-6f );
		ENSURE_SMALL( manifold.points[0].separation - ( -0.8f - 0.25f ), 1.0e-5f );
	}

	// A capsule along Z, same place
	{
		b3Capsule capsule = { { 0.9f, 0.2f, -0.5f }, { 0.9f, 0.2f, 0.5f }, 0.25f };
		ENSURE( b3CollideHullFacesAndCapsule( &manifold, 32, &boxA.base, exposed, &capsule, b3Transform_identity ) );
		ENSURE( manifold.pointCount == 2 );
		ENSURE_SMALL( manifold.normal.y - 1.0f, 1.0e-6f );
		ENSURE_SMALL( manifold.points[0].separation - ( -0.8f - 0.25f ), 1.0e-5f );
	}

	// A box B deep in the +X face
	{
		b3BoxHull boxB = b3MakeBoxHull( 0.25f, 0.25f, 0.25f );
		b3Transform xf = { { 0.9f, 0.5f, 0.0f }, b3Quat_identity };
		b3SATCache cache = { 0 };
		ENSURE( b3CollideHullFaces( &manifold, 32, &boxA.base, exposed, &boxB.base, xf, &cache ) );
		ENSURE( manifold.pointCount == 4 );
		ENSURE_SMALL( manifold.normal.y - 1.0f, 1.0e-6f );
		ENSURE_SMALL( manifold.points[0].separation - ( 0.25f - 1.0f ), 1.0e-5f );
	}

	// Clear of every exposed face: nothing
	{
		b3Sphere sphere = { { 0.0f, 3.0f, 0.0f }, 0.25f };
		ENSURE( b3CollideHullFacesAndSphere( &manifold, 32, &boxA.base, exposed, &sphere, b3Transform_identity ) == false );
		ENSURE( manifold.pointCount == 0 );
	}

	return 0;
}

// The queries a character body makes through the plugin: casts and overlaps of a capsule and a box proxy against a static
// grid shape in a world, ray casts, and a closest-point query on the grid shape. None of them may loop or read bad memory.
typedef struct QueryCounts
{
	int hits;
	int calls;
} QueryCounts;

static float CountCast( b3ShapeId shapeId, b3Pos point, b3Vec3 normal, float fraction, uint64_t userMaterialId, int triangleIndex,
						int childIndex, void* context )
{
	B3_UNUSED( shapeId, point, userMaterialId, triangleIndex, childIndex );
	QueryCounts* counts = context;
	counts->calls += 1;
	counts->hits += 1;
	B3_UNUSED( normal );
	return fraction;
}

static bool CountOverlap( b3ShapeId shapeId, void* context )
{
	B3_UNUSED( shapeId );
	QueryCounts* counts = context;
	counts->hits += 1;
	return true;
}

static uint32_t NextRandom( uint32_t* state )
{
	uint32_t x = *state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*state = x;
	return x;
}

static float RandomRange( uint32_t* state, float low, float high )
{
	return low + ( high - low ) * ( (float)( NextRandom( state ) & 0xffffff ) / (float)0x1000000 );
}

static int VoxelGridCharacterQueries( void )
{
	// A floor, a wall, a turned wall, a railing and a crate as the game's test row has them: 5 x 1 x 4 cells of 2 m
	float floorBox[6] = { 0.0f, 0.0f, 0.0f, 2.0f, 0.2f, 2.0f };
	float wallBox[6] = { 0.0f, 0.0f, 0.8f, 2.0f, 2.0f, 1.2f };
	float turnedBox[6] = { 0.8f, 0.0f, 0.0f, 1.2f, 2.0f, 2.0f };
	float railingBoxes[12] = { 0.0f, 0.0f, 0.0f, 2.0f, 0.2f, 2.0f, 0.0f, 0.2f, 1.8f, 2.0f, 1.1f, 2.0f };
	float crateBox[6] = { 0.4f, 0.0f, 0.4f, 1.6f, 1.2f, 1.6f };
	b3VoxelGridModule* modules[5] = {
		b3CreateVoxelGridModule( floorBox, 1, 2.0f, CELL_VOXELS ), b3CreateVoxelGridModule( wallBox, 1, 2.0f, CELL_VOXELS ),
		b3CreateVoxelGridModule( turnedBox, 1, 2.0f, CELL_VOXELS ), b3CreateVoxelGridModule( railingBoxes, 2, 2.0f, CELL_VOXELS ),
		b3CreateVoxelGridModule( crateBox, 1, 2.0f, CELL_VOXELS ),
	};
	for ( int m = 0; m < 5; ++m )
	{
		ENSURE( modules[m] != NULL );
	}

	const int sx = 5, sz = 4;
	int padded[( 5 + 2 ) * 3 * ( 4 + 2 )];
	for ( int i = 0; i < ( 5 + 2 ) * 3 * ( 4 + 2 ); ++i )
	{
		padded[i] = -1;
	}

	for ( int k = 0; k < sz; ++k )
	{
		for ( int i = 0; i < sx; ++i )
		{
			int module = 0;
			if ( k == 3 )
			{
				module = 1;
			}
			else if ( i == 3 )
			{
				module = 2;
			}
			else if ( i == 4 && k == 1 )
			{
				module = 3;
			}
			else if ( i == 1 && k == 1 )
			{
				module = 4;
			}
			padded[( ( 0 + 1 ) * ( sz + 2 ) + ( k + 1 ) ) * ( sx + 2 ) + ( i + 1 )] = module;
		}
	}

	b3VoxelGridDef def = { 0 };
	def.cellCountX = sx;
	def.cellCountY = 1;
	def.cellCountZ = sz;
	def.cellMeters = 2.0f;
	def.cellVoxels = CELL_VOXELS;
	def.maxBoxesPerCell = 256;
	def.modules = modules;
	def.moduleCount = 5;
	def.paddedCells = padded;
	b3VoxelGrid* grid = b3CreateVoxelGrid( &def );
	ENSURE( grid != NULL );

	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );
	b3BodyDef bodyDef = b3DefaultBodyDef();
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	b3ShapeId shapeId = b3CreateVoxelGridShape( bodyId, &shapeDef, grid );

	// The recovery pass of a body test asks for the closest point on every shape its query overlaps
	b3Vec3 target = { 3.0f, 1.0f, 3.0f };
	b3Vec3 closest = b3Shape_GetClosestPoint( shapeId, target );
	ENSURE( b3IsValidVec3( closest ) );

	b3Vec3 capsulePoints[2] = { { 0.0f, 0.35f, 0.0f }, { 0.0f, 1.45f, 0.0f } };
	b3Vec3 boxPoints[8];
	for ( int c = 0; c < 8; ++c )
	{
		boxPoints[c] = ( b3Vec3 ){ ( c & 1 ) ? 0.3f : -0.3f, ( c & 2 ) ? 0.3f : -0.3f, ( c & 4 ) ? 0.3f : -0.3f };
	}

	b3ShapeProxy proxies[2] = { { capsulePoints, 2, 0.35f }, { boxPoints, 8, 0.0f } };
	b3QueryFilter filter = b3DefaultQueryFilter();

	// The sunk crate slides and the moves of the Godot test, then a dense random set. Zero components, axis-aligned
	// motions and starts just inside boxes are drawn on purpose.
	uint32_t state = 0x1234abcdu;
	int total = 0;
	for ( int sample = 0; sample < 40000; ++sample )
	{
		b3Vec3 start = { RandomRange( &state, -1.0f, 11.0f ), RandomRange( &state, -0.5f, 3.0f ), RandomRange( &state, -1.0f, 9.0f ) };
		b3Vec3 motion = { RandomRange( &state, -4.0f, 4.0f ), RandomRange( &state, -3.0f, 3.0f ), RandomRange( &state, -4.0f, 4.0f ) };
		uint32_t pick = NextRandom( &state );
		if ( sample < 6 )
		{
			// Seam slides
			start = ( b3Vec3 ){ 0.5f + 0.5f * (float)( sample % 3 ), 0.4995f + 0.0f, 1.0f };
			motion = ( b3Vec3 ){ 5.0f, 0.0f, 0.0f };
		}
		if ( pick & 1 )
		{
			start.x = 0.5f * floorf( start.x * 2.0f );
		}
		if ( pick & 2 )
		{
			start.y = 0.05f * floorf( start.y * 20.0f ) + ( ( pick & 4 ) ? 0.4995f - 0.5f : 0.0f );
		}
		if ( pick & 8 )
		{
			motion.x = 0.0f;
		}
		if ( pick & 16 )
		{
			motion.y = 0.0f;
		}
		if ( pick & 32 )
		{
			motion.z = 0.0f;
		}
		if ( ( pick & 0x3c0 ) == 0 )
		{
			motion = b3Vec3_zero;
		}

		for ( int p = 0; p < 2; ++p )
		{
			b3Vec3 points[8];
			for ( int i = 0; i < proxies[p].count; ++i )
			{
				points[i] = b3Add( proxies[p].points[i], start );
			}
			b3ShapeProxy proxy = { points, proxies[p].count, proxies[p].radius };

			QueryCounts counts = { 0, 0 };
			b3World_CastShape( worldId, b3Vec3_zero, &proxy, motion, filter, CountCast, &counts );
			b3World_OverlapShape( worldId, b3Vec3_zero, &proxy, filter, CountOverlap, &counts );
			total += counts.hits;
		}

		b3World_CastRay( worldId, start, motion, filter, CountCast, &( QueryCounts ){ 0, 0 } );
		b3Shape_GetClosestPoint( shapeId, start );
	}
	printf( "  character queries: %d hits\n", total );

	b3DestroyWorld( worldId );
	b3ReleaseVoxelGrid( grid );
	for ( int m = 0; m < 5; ++m )
	{
		b3ReleaseVoxelGridModule( modules[m] );
	}
	return 0;
}


// Dynamic voxel grids: shards made of 0.5 m boxes

#define SHARD_CELL 0.5f
#define SHARD_CELL_VOXELS 5

typedef struct Shard
{
	b3VoxelGridModule* module;
	b3VoxelGrid* grid;
	int sizeX, sizeY, sizeZ;
} Shard;

static int ShardPaddedIndex( const Shard* shard, int i, int j, int k )
{
	return ( ( j + 1 ) * ( shard->sizeZ + 2 ) + ( k + 1 ) ) * ( shard->sizeX + 2 ) + ( i + 1 );
}

// A block of sizeX x sizeY x sizeZ cells, each one solid 0.5 m box. The grid origin is the body origin.
static bool MakeShard( Shard* shard, int sizeX, int sizeY, int sizeZ )
{
	shard->sizeX = sizeX;
	shard->sizeY = sizeY;
	shard->sizeZ = sizeZ;

	float box[6] = { 0.0f, 0.0f, 0.0f, SHARD_CELL, SHARD_CELL, SHARD_CELL };
	shard->module = b3CreateVoxelGridModule( box, 1, SHARD_CELL, SHARD_CELL_VOXELS );
	if ( shard->module == NULL )
	{
		return false;
	}

	int paddedCount = ( sizeX + 2 ) * ( sizeY + 2 ) * ( sizeZ + 2 );
	int* padded = (int*)b3Alloc( paddedCount * sizeof( int ) );
	for ( int i = 0; i < paddedCount; ++i )
	{
		padded[i] = -1;
	}

	for ( int j = 0; j < sizeY; ++j )
	{
		for ( int k = 0; k < sizeZ; ++k )
		{
			for ( int i = 0; i < sizeX; ++i )
			{
				padded[ShardPaddedIndex( shard, i, j, k )] = 0;
			}
		}
	}

	b3VoxelGridDef def = { 0 };
	def.cellCountX = sizeX;
	def.cellCountY = sizeY;
	def.cellCountZ = sizeZ;
	def.origin = b3Vec3_zero;
	def.cellMeters = SHARD_CELL;
	def.cellVoxels = SHARD_CELL_VOXELS;
	def.maxBoxesPerCell = 1;
	def.modules = &shard->module;
	def.moduleCount = 1;
	def.paddedCells = padded;
	shard->grid = b3CreateVoxelGrid( &def );
	b3Free( padded, paddedCount * sizeof( int ) );
	return shard->grid != NULL;
}

static void DestroyShard( Shard* shard )
{
	b3ReleaseVoxelGrid( shard->grid );
	b3ReleaseVoxelGridModule( shard->module );
}

static b3BodyId CreateShardBody( b3WorldId worldId, const Shard* shard, b3Vec3 position, float friction, b3ShapeId* shapeId )
{
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = position;
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );

	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = friction;
	b3ShapeId id = b3CreateVoxelGridShape( bodyId, &shapeDef, shard->grid );
	if ( shapeId != NULL )
	{
		*shapeId = id;
	}
	return bodyId;
}

// The mass, center and inertia of a grid on a dynamic body come from its boxes
static int VoxelGridDynamicMass( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3ShapeId shapeId;
	b3BodyId bodyId = CreateShardBody( worldId, &shard, ( b3Vec3 ){ 10.0f, 5.0f, 10.0f }, 0.6f, &shapeId );
	ENSURE( B3_IS_NON_NULL( shapeId ) );

	// 4 x 2 x 4 m of the default density
	float density = b3DefaultShapeDef().density;
	float mass = density * 4.0f * 2.0f * 4.0f;
	ENSURE_SMALL( b3Body_GetMass( bodyId ) - mass, 1.0e-3f * mass );

	b3Vec3 center = b3Body_GetLocalCenter( bodyId );
	ENSURE_SMALL( center.x - 2.0f, 1.0e-4f );
	ENSURE_SMALL( center.y - 1.0f, 1.0e-4f );
	ENSURE_SMALL( center.z - 2.0f, 1.0e-4f );

	b3Matrix3 inertia = b3Body_GetLocalRotationalInertia( bodyId );
	float expectedX = mass / 12.0f * ( 2.0f * 2.0f + 4.0f * 4.0f );
	float expectedY = mass / 12.0f * ( 4.0f * 4.0f + 4.0f * 4.0f );
	ENSURE_SMALL( inertia.cx.x - expectedX, 1.0e-3f * expectedX );
	ENSURE_SMALL( inertia.cy.y - expectedY, 1.0e-3f * expectedY );
	ENSURE_SMALL( inertia.cz.z - expectedX, 1.0e-3f * expectedX );
	printf( "  shard mass %.0f kg, inertia %.0f %.0f %.0f\n", b3Body_GetMass( bodyId ), inertia.cx.x, inertia.cy.y, inertia.cz.z );

	// The AABB is that of the boxes, not of the grid's cell array (they coincide here), and it follows the body
	b3AABB aabb = b3Shape_GetAABB( shapeId );
	ENSURE( aabb.lowerBound.x < 10.0f + 0.1f && aabb.upperBound.x > 14.0f - 0.1f );

	DestroyShard( &shard );
	b3DestroyWorld( worldId );
	return 0;
}

// A shard dropped on a static hull ground rests there and falls asleep
static int VoxelGridShardOnHull( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BodyDef groundDef = b3DefaultBodyDef();
	groundDef.position = ( b3Vec3 ){ 0.0f, -0.5f, 0.0f };
	b3BodyId groundId = b3CreateBody( worldId, &groundDef );
	b3BoxHull ground = b3MakeBoxHull( 20.0f, 0.5f, 20.0f );
	b3ShapeDef groundShape = b3DefaultShapeDef();
	groundShape.baseMaterial.friction = 0.6f;
	b3CreateHullShape( groundId, &groundShape, &ground.base );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3BodyId bodyId = CreateShardBody( worldId, &shard, ( b3Vec3 ){ -2.0f, 0.3f, -2.0f }, 0.6f, NULL );

	for ( int i = 0; i < 300; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	printf( "  shard on hull: position %.4f %.4f %.4f, awake %d\n", p.x, p.y, p.z, b3Body_IsAwake( bodyId ) );
	ENSURE_SMALL( p.y, 0.003f );
	ENSURE_SMALL( p.x + 2.0f, 0.01f );
	ENSURE_SMALL( p.z + 2.0f, 0.01f );
	ENSURE( b3Body_IsAwake( bodyId ) == false );

	DestroyShard( &shard );
	b3DestroyWorld( worldId );
	return 0;
}


// A shard dropped on a static grid floor settles on it, stays put and falls asleep
static int VoxelGridShardOnGrid( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3BodyId bodyId = CreateShardBody( world.worldId, &shard, ( b3Vec3 ){ 2.0f, SLAB_HEIGHT + 0.3f, 2.0f }, 0.6f, NULL );

	for ( int i = 0; i < 360; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	b3Quat q = b3Body_GetRotation( bodyId );
	printf( "  shard on grid: position %.4f %.4f %.4f, tilt %.5f, awake %d\n", p.x, p.y, p.z, fabsf( q.v.x ) + fabsf( q.v.y ) + fabsf( q.v.z ),
			b3Body_IsAwake( bodyId ) );
	ENSURE_SMALL( p.x - 2.0f, 0.01f );
	ENSURE_SMALL( p.z - 2.0f, 0.01f );
	ENSURE_SMALL( p.y - SLAB_HEIGHT, 0.003f );
	ENSURE( b3Body_IsAwake( bodyId ) == false );

	DestroyShard( &shard );
	DestroyWorld( &world );
	return 0;
}

// A shard on a shard on a floor, the upper one offset by half a box so the boxes do not line up
static int VoxelGridShardStack( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3BodyId lowerId = CreateShardBody( world.worldId, &shard, ( b3Vec3 ){ 2.0f, SLAB_HEIGHT + 0.01f, 2.0f }, 0.6f, NULL );
	b3BodyId upperId = CreateShardBody( world.worldId, &shard, ( b3Vec3 ){ 2.25f, SLAB_HEIGHT + 2.0f + 0.02f, 2.25f }, 0.6f, NULL );

	for ( int i = 0; i < 480; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos lower = b3Body_GetPosition( lowerId );
	b3Pos upper = b3Body_GetPosition( upperId );
	printf( "  shard stack: lower %.4f %.4f %.4f, upper %.4f %.4f %.4f, awake %d %d\n", lower.x, lower.y, lower.z, upper.x, upper.y,
			upper.z, b3Body_IsAwake( lowerId ), b3Body_IsAwake( upperId ) );
	ENSURE_SMALL( lower.x - 2.0f, 0.01f );
	ENSURE_SMALL( lower.z - 2.0f, 0.01f );
	ENSURE_SMALL( lower.y - SLAB_HEIGHT, 0.004f );
	ENSURE_SMALL( upper.x - 2.25f, 0.01f );
	ENSURE_SMALL( upper.z - 2.25f, 0.01f );
	ENSURE_SMALL( upper.y - ( SLAB_HEIGHT + 2.0f ), 0.006f );
	ENSURE( b3Body_IsAwake( lowerId ) == false );
	ENSURE( b3Body_IsAwake( upperId ) == false );

	DestroyShard( &shard );
	DestroyWorld( &world );
	return 0;
}

typedef enum ShardKind
{
	SHARD_GRID,
	SHARD_BOXES,
	SHARD_SOLID
} ShardKind;

typedef struct PileResult
{
	float worstStepMs;
	float meanStepMs;
	float worstCollideMs;
	float worstSolveMs;
	int contactCount;
	int touchingManifolds;
	int touchingContacts;
	float maxDrift;
} PileResult;

// 27 shards of 8 x 4 x 8 boxes land on a grid floor, touching each other, in one of three representations
static int RunShardPile( ShardKind kind, PileResult* result )
{
	World world;
	CreateWorld( &world, 8, 0.6f );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );

	b3BoxHull solidHull = b3MakeOffsetBoxHull( 2.0f, 1.0f, 2.0f, ( b3Vec3 ){ 2.0f, 1.0f, 2.0f } );

	b3BodyId bodies[27];
	b3Vec3 starts[27];
	int count = 0;
	for ( int j = 0; j < 3; ++j )
	{
		for ( int k = 0; k < 3; ++k )
		{
			for ( int i = 0; i < 3; ++i )
			{
				b3Vec3 position = { 1.0f + (float)i * 4.01f, SLAB_HEIGHT + 0.02f + (float)j * 2.01f, 1.0f + (float)k * 4.01f };
				starts[count] = position;

				b3ShapeDef shapeDef = b3DefaultShapeDef();
				shapeDef.baseMaterial.friction = 0.6f;
				if ( kind == SHARD_GRID )
				{
					bodies[count] = CreateShardBody( world.worldId, &shard, position, 0.6f, NULL );
				}
				else
				{
					b3BodyDef bodyDef = b3DefaultBodyDef();
					bodyDef.type = b3_dynamicBody;
					bodyDef.position = position;
					bodies[count] = b3CreateBody( world.worldId, &bodyDef );
					if ( kind == SHARD_SOLID )
					{
						b3CreateHullShape( bodies[count], &shapeDef, &solidHull.base );
					}
					else
					{
						shapeDef.updateBodyMass = false;
						for ( int y = 0; y < 4; ++y )
						{
							for ( int z = 0; z < 8; ++z )
							{
								for ( int x = 0; x < 8; ++x )
								{
									b3Vec3 offset = { ( (float)x + 0.5f ) * 0.5f, ( (float)y + 0.5f ) * 0.5f, ( (float)z + 0.5f ) * 0.5f };
									b3BoxHull boxHull = b3MakeOffsetBoxHull( 0.25f, 0.25f, 0.25f, offset );
									b3CreateHullShape( bodies[count], &shapeDef, &boxHull.base );
								}
							}
						}
						b3Body_ApplyMassFromShapes( bodies[count] );
					}
				}
				count += 1;
			}
		}
	}

	memset( result, 0, sizeof( *result ) );
	const int stepCount = 90;
	float total = 0.0f;
	for ( int step = 0; step < stepCount; ++step )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
		b3Profile profile = b3World_GetProfile( world.worldId );
		total += profile.step;
		if ( profile.step > result->worstStepMs )
		{
			result->worstStepMs = profile.step;
			result->worstCollideMs = profile.collide;
			result->worstSolveMs = profile.solve;
		}

		b3Counters counters = b3World_GetCounters( world.worldId );
		result->contactCount = b3MaxInt( result->contactCount, counters.contactCount );
	}
	result->meanStepMs = total / (float)stepCount;

	// The touching contacts and manifolds at the end
	static b3ContactData data[2048];
	for ( int n = 0; n < 27; ++n )
	{
		int pairCount = b3Body_GetContactData( bodies[n], data, 2048 );
		for ( int m = 0; m < pairCount; ++m )
		{
			result->touchingContacts += 1;
			result->touchingManifolds += data[m].manifoldCount;
		}

		b3Pos p = b3Body_GetPosition( bodies[n] );
		float drift = fabsf( p.x - starts[n].x ) + fabsf( p.z - starts[n].z );
		result->maxDrift = b3MaxFloat( result->maxDrift, drift );
		ENSURE( p.y > starts[n].y - 0.1f );
	}

	DestroyShard( &shard );
	DestroyWorld( &world );
	return 0;
}

static int VoxelGridShardPile( void )
{
	PileResult grid, boxes, solid;
	ENSURE( RunShardPile( SHARD_GRID, &grid ) == 0 );
	ENSURE( RunShardPile( SHARD_BOXES, &boxes ) == 0 );
	ENSURE( RunShardPile( SHARD_SOLID, &solid ) == 0 );

	const char* names[3] = { "grid shards", "per-box shards", "solid boxes" };
	const PileResult* results[3] = { &grid, &boxes, &solid };
	for ( int i = 0; i < 3; ++i )
	{
		const PileResult* r = results[i];
		printf( "  pile %-14s: worst step %7.3f ms (collide %.3f solve %.3f), mean %.3f ms, contacts %d, touching contacts %d, manifolds %d, "
				"drift %.3f\n",
				names[i], r->worstStepMs, r->worstCollideMs, r->worstSolveMs, r->meanStepMs, r->contactCount, r->touchingContacts,
				r->touchingManifolds, r->maxDrift );
	}

	// A grid shard is one shape and its contacts are few
	ENSURE( grid.contactCount < 400 );
	ENSURE( grid.worstStepMs < 0.5f * boxes.worstStepMs );
	return 0;
}


// A shard dropped on a static mesh floor (a triangulated plane, optionally tilted) settles and falls asleep
static int ShardOnMesh( float tiltRadians, float* heightAbove )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3MeshData* mesh = b3CreateGridMesh( 16, 16, 1.0f, 0, true );
	b3BodyDef floorDef = b3DefaultBodyDef();
	floorDef.rotation = b3MakeQuatFromAxisAngle( b3Vec3_axisZ, tiltRadians );
	b3BodyId floorId = b3CreateBody( worldId, &floorDef );
	b3ShapeDef floorShape = b3DefaultShapeDef();
	floorShape.baseMaterial.friction = 0.6f;
	b3CreateMeshShape( floorId, &floorShape, mesh, b3Vec3_one );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3BodyId bodyId = CreateShardBody( worldId, &shard, ( b3Vec3 ){ -2.0f, 0.3f + 2.0f * sinf( tiltRadians ), -2.0f }, 0.6f, NULL );

	for ( int i = 0; i < 480; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	b3Quat q = b3Body_GetRotation( bodyId );
	b3Vec3 up = b3RotateVector( q, b3Vec3_axisY );
	printf( "  shard on mesh (tilt %.2f): position %.4f %.4f %.4f, up %.4f %.4f %.4f, awake %d\n", tiltRadians, p.x, p.y, p.z, up.x,
			up.y, up.z, b3Body_IsAwake( bodyId ) );
	*heightAbove = p.y;

	b3Vec3 expectedUp = b3RotateVector( b3MakeQuatFromAxisAngle( b3Vec3_axisZ, tiltRadians ), b3Vec3_axisY );
	ENSURE( b3Dot( up, expectedUp ) > 0.999f );
	ENSURE( b3Body_IsAwake( bodyId ) == false );
	ENSURE( b3Length( b3Body_GetLinearVelocity( bodyId ) ) < 0.01f );

	DestroyShard( &shard );
	b3DestroyWorld( worldId );
	b3DestroyMesh( mesh );
	return 0;
}

static int VoxelGridShardOnMesh( void )
{
	float height;
	ENSURE( ShardOnMesh( 0.0f, &height ) == 0 );
	ENSURE_SMALL( height, 0.004f );

	// On a slope the friction holds it: it stays within a box width of where it landed
	ENSURE( ShardOnMesh( 0.2f, &height ) == 0 );
	return 0;
}


// A static slab of cells, each one solid box of the given thickness, 16 x 16 cells of 0.5 m. Its top is at y = 0.
static bool MakeSlab( Shard* slab, float thickness )
{
	slab->sizeX = 16;
	slab->sizeY = 1;
	slab->sizeZ = 16;

	float box[6] = { 0.0f, 0.0f, 0.0f, SHARD_CELL, thickness, SHARD_CELL };
	slab->module = b3CreateVoxelGridModule( box, 1, SHARD_CELL, SHARD_CELL_VOXELS );
	if ( slab->module == NULL )
	{
		return false;
	}

	int paddedCount = 18 * 3 * 18;
	int* padded = (int*)b3Alloc( paddedCount * sizeof( int ) );
	for ( int i = 0; i < paddedCount; ++i )
	{
		padded[i] = -1;
	}

	for ( int k = 0; k < 16; ++k )
	{
		for ( int i = 0; i < 16; ++i )
		{
			padded[ShardPaddedIndex( slab, i, 0, k )] = 0;
		}
	}

	b3VoxelGridDef def = { 0 };
	def.cellCountX = 16;
	def.cellCountY = 1;
	def.cellCountZ = 16;
	def.origin = ( b3Vec3 ){ -4.0f, -thickness, -4.0f };
	def.cellMeters = SHARD_CELL;
	def.cellVoxels = SHARD_CELL_VOXELS;
	def.maxBoxesPerCell = 1;
	def.modules = &slab->module;
	def.moduleCount = 1;
	def.paddedCells = padded;
	slab->grid = b3CreateVoxelGrid( &def );
	b3Free( padded, paddedCount * sizeof( int ) );
	return slab->grid != NULL;
}

// A flat shard launched at 30 m/s at a 0.2 m slab moves more than the slab is thick in one step
static int FastShard( bool continuous, float speed, float startY, float* lowestY )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.enableContinuous = continuous;
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Shard slab;
	ENSURE( MakeSlab( &slab, 0.2f ) );
	b3BodyDef slabBodyDef = b3DefaultBodyDef();
	b3BodyId slabBodyId = b3CreateBody( worldId, &slabBodyDef );
	b3ShapeDef slabShapeDef = b3DefaultShapeDef();
	b3CreateVoxelGridShape( slabBodyId, &slabShapeDef, slab.grid );

	// 4 x 0.5 x 4 m above the slab
	Shard shard;
	ENSURE( MakeShard( &shard, 8, 1, 8 ) );
	b3BodyId bodyId = CreateShardBody( worldId, &shard, ( b3Vec3 ){ -2.0f, startY, -2.0f }, 0.6f, NULL );
	b3Body_SetLinearVelocity( bodyId, ( b3Vec3 ){ 0.0f, -speed, 0.0f } );

	*lowestY = FLT_MAX;
	for ( int i = 0; i < 120; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		*lowestY = b3MinFloat( *lowestY, (float)b3Body_GetPosition( bodyId ).y );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	printf( "  fast shard (continuous %d, %.0f m/s): lowest y %.3f, final %.3f %.3f %.3f, awake %d\n", continuous, speed, *lowestY, p.x,
			p.y, p.z, b3Body_IsAwake( bodyId ) );

	DestroyShard( &shard );
	b3DestroyWorld( worldId );
	DestroyShard( &slab );
	return 0;
}

static int VoxelGridFastShard( void )
{
	// At 60 m/s a step moves the shard a metre. Without continuous collision it goes through the slab.
	float lowest;
	ENSURE( FastShard( false, 60.0f, 0.4f, &lowest ) == 0 );
	ENSURE( lowest < -1.0f );

	// With it the shard lands on top and rests there, at 60 m/s and at 30 m/s
	ENSURE( FastShard( true, 60.0f, 0.4f, &lowest ) == 0 );
	ENSURE( lowest > -0.01f );
	ENSURE( FastShard( true, 30.0f, 1.0f, &lowest ) == 0 );
	ENSURE( lowest > -0.01f );
	return 0;
}

// A tilted shard dropped on a shard rolls onto a face and rests flat: the boxes of two bodies that are not aligned still touch
static int VoxelGridShardTilt( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3BodyId lowerId = CreateShardBody( world.worldId, &shard, ( b3Vec3 ){ 2.0f, SLAB_HEIGHT + 0.01f, 2.0f }, 0.6f, NULL );
	b3BodyId upperId = CreateShardBody( world.worldId, &shard, ( b3Vec3 ){ 2.0f, SLAB_HEIGHT + 2.4f, 2.0f }, 0.6f, NULL );
	b3Vec3 axis = b3Normalize( ( b3Vec3 ){ 1.0f, 0.0f, 0.6f } );
	b3Body_SetTransform( upperId, ( b3Pos ){ 2.1, SLAB_HEIGHT + 2.4, 1.9 }, b3MakeQuatFromAxisAngle( axis, 0.12f ) );

	for ( int i = 0; i < 900; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos lower = b3Body_GetPosition( lowerId );
	b3Pos upper = b3Body_GetPosition( upperId );
	b3Vec3 up = b3RotateVector( b3Body_GetRotation( upperId ), b3Vec3_axisY );
	printf( "  shard tilt: lower y %.4f, upper %.4f %.4f %.4f, up %.5f %.5f %.5f, awake %d %d\n", lower.y, upper.x, upper.y, upper.z, up.x,
			up.y, up.z, b3Body_IsAwake( lowerId ), b3Body_IsAwake( upperId ) );
	ENSURE_SMALL( lower.y - SLAB_HEIGHT, 0.004f );
	ENSURE( up.y > 0.9998f );
	ENSURE_SMALL( upper.y - ( SLAB_HEIGHT + 2.0f ), 0.01f );
	ENSURE( b3Body_IsAwake( lowerId ) == false );
	ENSURE( b3Body_IsAwake( upperId ) == false );

	DestroyShard( &shard );
	DestroyWorld( &world );
	return 0;
}

// Queries against a grid that moves: its body is turned a quarter circle about Y and moved, so the shard covers
// x 5..9, z -1..3, y 1..3 in the world
static int VoxelGridMovingQueries( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3ShapeId shapeId;
	b3BodyId bodyId = CreateShardBody( worldId, &shard, ( b3Vec3 ){ 5.0f, 1.0f, 3.0f }, 0.6f, &shapeId );
	b3Body_SetTransform( bodyId, ( b3Pos ){ 5.0, 1.0, 3.0 }, b3MakeQuatFromAxisAngle( b3Vec3_axisY, 0.5f * B3_PI ) );
	b3QueryFilter filter = b3DefaultQueryFilter();

	// A ray down onto the top, and one that would hit the shard if its turn were ignored
	b3RayResult hit = b3World_CastRayClosest( worldId, ( b3Pos ){ 7.0, 10.0, 1.0 }, ( b3Vec3 ){ 0.0f, -20.0f, 0.0f }, filter );
	ENSURE( hit.hit );
	ENSURE_SMALL( hit.point.y - 3.0, 1.0e-4 );
	ENSURE_SMALL( hit.normal.y - 1.0f, 1.0e-4f );
	hit = b3World_CastRayClosest( worldId, ( b3Pos ){ 2.0, 10.0, 2.0 }, ( b3Vec3 ){ 0.0f, -20.0f, 0.0f }, filter );
	ENSURE( hit.hit == false );

	// A ray from the side into the turned face: local +X is world -Z, so the face at local x = 0 is the world z = 3 side
	hit = b3World_CastRayClosest( worldId, ( b3Pos ){ 7.0, 2.0, 8.0 }, ( b3Vec3 ){ 0.0f, 0.0f, -10.0f }, filter );
	ENSURE( hit.hit );
	ENSURE_SMALL( hit.point.z - 3.0, 1.0e-4 );
	ENSURE_SMALL( hit.normal.z - 1.0f, 1.0e-4f );

	// A sphere overlaps inside, not where an unturned shard would be
	b3Vec3 center = { 0.0f, 0.0f, 0.0f };
	b3ShapeProxy proxy = { &center, 1, 0.3f };
	QueryCounts counts = { 0 };
	b3World_OverlapShape( worldId, ( b3Pos ){ 7.0, 2.0, 1.0 }, &proxy, filter, CountOverlap, &counts );
	ENSURE( counts.hits == 1 );
	counts.hits = 0;
	b3World_OverlapShape( worldId, ( b3Pos ){ 2.0, 2.0, 2.0 }, &proxy, filter, CountOverlap, &counts );
	ENSURE( counts.hits == 0 );

	// A sphere cast across the top lands on it
	counts = ( QueryCounts ){ 0 };
	b3World_CastShape( worldId, ( b3Pos ){ 7.0, 6.0, 1.0 }, &proxy, ( b3Vec3 ){ 0.0f, -6.0f, 0.0f }, filter, CountCast, &counts );
	ENSURE( counts.hits == 1 );

	DestroyShard( &shard );
	b3DestroyWorld( worldId );
	return 0;
}

// Removing cells of a resting shard in place: its mass follows, it wakes and settles again
static int VoxelGridShardEdit( void )
{
	World world;
	CreateWorld( &world, 4, 0.6f );

	Shard shard;
	ENSURE( MakeShard( &shard, 8, 4, 8 ) );
	b3ShapeId shapeId;
	b3BodyId bodyId = CreateShardBody( world.worldId, &shard, ( b3Vec3 ){ 2.0f, SLAB_HEIGHT + 0.01f, 2.0f }, 0.6f, &shapeId );
	float density = b3DefaultShapeDef().density;
	float boxMass = density * SHARD_CELL * SHARD_CELL * SHARD_CELL;

	for ( int i = 0; i < 240; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}
	ENSURE( b3Body_IsAwake( bodyId ) == false );
	float mass = b3Body_GetMass( bodyId );
	b3Vec3 center = b3Body_GetLocalCenter( bodyId );

	// A box at the top corner
	int cells[1] = { ShardPaddedIndex( &shard, 0, 3, 0 ) };
	int modules[1] = { -1 };
	b3Shape_VoxelGridSetCells( shapeId, cells, modules, 1 );
	float massAfterTop = b3Body_GetMass( bodyId );
	b3Vec3 centerAfterTop = b3Body_GetLocalCenter( bodyId );
	printf( "  shard edit: mass %.1f -> %.1f (a box is %.1f), center %.4f %.4f %.4f -> %.4f %.4f %.4f\n", mass, massAfterTop, boxMass,
			center.x, center.y, center.z, centerAfterTop.x, centerAfterTop.y, centerAfterTop.z );
	ENSURE_SMALL( massAfterTop - ( mass - boxMass ), 1.0e-3f * boxMass );
	ENSURE( centerAfterTop.x > center.x && centerAfterTop.y < center.y && centerAfterTop.z > center.z );
	ENSURE( b3Body_IsAwake( bodyId ) );

	// And one at the bottom corner, under the weight
	cells[0] = ShardPaddedIndex( &shard, 7, 0, 7 );
	b3Shape_VoxelGridSetCells( shapeId, cells, modules, 1 );
	ENSURE_SMALL( b3Body_GetMass( bodyId ) - ( mass - 2.0f * boxMass ), 1.0e-3f * boxMass );

	for ( int i = 0; i < 360; ++i )
	{
		b3World_Step( world.worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	printf( "  shard edit: settled at %.4f %.4f %.4f, awake %d\n", p.x, p.y, p.z, b3Body_IsAwake( bodyId ) );
	ENSURE_SMALL( p.x - 2.0f, 0.02f );
	ENSURE_SMALL( p.z - 2.0f, 0.02f );
	ENSURE_SMALL( p.y - SLAB_HEIGHT, 0.01f );
	ENSURE( b3Body_IsAwake( bodyId ) == false );

	DestroyShard( &shard );
	DestroyWorld( &world );
	return 0;
}

int VoxelGridTest( void )
{
	RUN_SUBTEST( VoxelGridModuleCover );
	RUN_SUBTEST( VoxelGridRayCast );
	RUN_SUBTEST( VoxelGridExposedFaces );
	RUN_SUBTEST( VoxelGridShapeCast );
	RUN_SUBTEST( VoxelGridShapeLifetime );
	RUN_SUBTEST( VoxelGridRestAndSleep );
	RUN_SUBTEST( VoxelGridSlideAcrossSeams );
	RUN_SUBTEST( VoxelGridSlideAcrossSeamsHeavy );
	RUN_SUBTEST( VoxelGridRestingImpulseGrid );
	RUN_SUBTEST( VoxelGridRestingImpulseBox );
	RUN_SUBTEST( VoxelGridLandingImpulse );
	RUN_SUBTEST( VoxelGridSetCellWakes );
	RUN_SUBTEST( VoxelGridSetCellRestoresSupport );
	RUN_SUBTEST( VoxelGridRecovery );
	RUN_SUBTEST( VoxelGridCharacterQueries );
	RUN_SUBTEST( VoxelGridDynamicMass );
	RUN_SUBTEST( VoxelGridShardOnHull );
	RUN_SUBTEST( VoxelGridShardOnGrid );
	RUN_SUBTEST( VoxelGridShardStack );
	RUN_SUBTEST( VoxelGridShardPile );
	RUN_SUBTEST( VoxelGridShardOnMesh );
	RUN_SUBTEST( VoxelGridShardTilt );
	RUN_SUBTEST( VoxelGridFastShard );
	RUN_SUBTEST( VoxelGridMovingQueries );
	RUN_SUBTEST( VoxelGridShardEdit );
	return 0;
}
