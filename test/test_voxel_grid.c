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
	ENSURE_SMALL( p.y - ( SLAB_HEIGHT + half ), 0.02f );
	ENSURE( b3Body_IsAwake( boxId ) == false );

	DestroyWorld( &world );
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
	RUN_SUBTEST( VoxelGridSetCellWakes );
	RUN_SUBTEST( VoxelGridSetCellRestoresSupport );
	RUN_SUBTEST( VoxelGridCharacterQueries );
	return 0;
}
