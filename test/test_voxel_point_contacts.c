// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

// Terra: voxel grid against voxel grid contacts, box pairs against point samples (b3World_EnableVoxelPointContacts), on
// building-like grids at the game's grain (2 m cells of 20 voxels): a tower block resting on and dropped onto a ground
// grid, a block sliding over the ground's cell seams, a fast shard, and two blocks landing on each other. Prints a line per
// case and mode; run on request with the filter VoxelPointBench.

#include "test_macros.h"

#include "core.h"

#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/math_functions.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct PointGrid
{
	b3VoxelGridModule* module;
	b3VoxelGrid* grid;
} PointGrid;

// sizeX x sizeY x sizeZ cells of one module, the grid origin at origin
static bool MakePointGrid( PointGrid* out, const float* boxes, int boxCount, int sizeX, int sizeY, int sizeZ, b3Vec3 origin )
{
	out->module = b3CreateVoxelGridModule( boxes, boxCount, 2.0f, 20 );
	if ( out->module == NULL )
	{
		return false;
	}

	int paddedCount = ( sizeX + 2 ) * ( sizeY + 2 ) * ( sizeZ + 2 );
	int* padded = (int*)b3Alloc( paddedCount * (int)sizeof( int ) );
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
				padded[( ( j + 1 ) * ( sizeZ + 2 ) + ( k + 1 ) ) * ( sizeX + 2 ) + ( i + 1 )] = 0;
			}
		}
	}

	b3VoxelGridDef def = { 0 };
	def.cellCountX = sizeX;
	def.cellCountY = sizeY;
	def.cellCountZ = sizeZ;
	def.origin = origin;
	def.cellMeters = 2.0f;
	def.cellVoxels = 20;
	def.maxBoxesPerCell = boxCount;
	def.modules = &out->module;
	def.moduleCount = 1;
	def.paddedCells = padded;
	out->grid = b3CreateVoxelGrid( &def );
	b3Free( padded, paddedCount * (int)sizeof( int ) );
	return out->grid != NULL;
}

static void DestroyPointGrid( PointGrid* g )
{
	b3ReleaseVoxelGrid( g->grid );
	b3ReleaseVoxelGridModule( g->module );
}

// A storey cell: a 0.2 m floor slab and four 0.2 m walls, 5 boxes
static const float s_storeyBoxes[] = {
	0.0f, 0.0f, 0.0f, 2.0f, 0.2f, 2.0f, // slab
	0.0f, 0.2f, 0.0f, 0.2f, 2.0f, 2.0f, // -x wall
	1.8f, 0.2f, 0.0f, 2.0f, 2.0f, 2.0f, // +x wall
	0.2f, 0.2f, 0.0f, 1.8f, 2.0f, 0.2f, // -z wall
	0.2f, 0.2f, 1.8f, 1.8f, 2.0f, 2.0f, // +z wall
};

// A ground cell: solid to 1 m in four 1 m columns (seams inside the cell as well as between cells)
static const float s_groundBoxes[] = {
	0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 2.0f, 1.0f, 1.0f,
	0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 1.0f, 2.0f,
};

// A solid cell, one box
static const float s_solidBoxes[] = { 0.0f, 0.0f, 0.0f, 2.0f, 2.0f, 2.0f };

typedef struct Run
{
	float meanStep, maxStep, maxCollide, maxSolve;
	float minSeparation;
	int maxPoints, maxSamples, maxSat, maxManifolds, maxPairs;
	int sleptAt;
	float maxSpeedLate;
	b3Vec3 endPosition;
	float maxAngularSpeed;
} Run;

static int CountManifolds( b3ShapeId shapeId, bool print )
{
	b3ContactData data[16];
	int count = b3Shape_GetContactData( shapeId, data, 16 );
	int manifolds = 0;
	for ( int c = 0; c < count; ++c )
	{
		manifolds += data[c].manifoldCount;
		if ( print )
		{
			for ( int m = 0; m < data[c].manifoldCount; ++m )
			{
				const b3Manifold* mf = data[c].manifolds + m;
				float least = FLT_MAX, most = 0.0f;
				for ( int p = 0; p < mf->pointCount; ++p )
				{
					least = b3MinFloat( least, mf->points[p].separation );
					most = b3MaxFloat( most, mf->points[p].totalNormalImpulse );
				}
				printf( "      manifold n=(%.2f %.2f %.2f) points=%d min_sep=%.3f max_total_impulse=%.0f\n", mf->normal.x, mf->normal.y,
						mf->normal.z, mf->pointCount, least, most );
			}
		}
	}
	return manifolds;
}

static void StepAndMeasure( b3WorldId worldId, b3BodyId bodyId, b3ShapeId shapeId, int steps, int lateFrom, Run* run )
{
	*run = (Run){ 0 };
	run->minSeparation = FLT_MAX;
	run->sleptAt = -1;
	float total = 0.0f;
	for ( int i = 0; i < steps; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		b3Profile p = b3World_GetProfile( worldId );
		b3Counters c = b3World_GetCounters( worldId );
		total += p.step;
		run->maxStep = b3MaxFloat( run->maxStep, p.step );
		run->maxCollide = b3MaxFloat( run->maxCollide, p.collide );
		run->maxSolve = b3MaxFloat( run->maxSolve, p.solve );
		if ( c.voxelGridPairCount > 0 )
		{
			run->minSeparation = b3MinFloat( run->minSeparation, c.voxelMinSeparation );
		}
		run->maxPoints = b3MaxInt( run->maxPoints, c.voxelPointCount );
		run->maxSamples = b3MaxInt( run->maxSamples, c.voxelSampleCount );
		run->maxSat = b3MaxInt( run->maxSat, c.satCallCount );
		run->maxPairs = b3MaxInt( run->maxPairs, c.voxelGridPairCount );
		run->maxManifolds = b3MaxInt( run->maxManifolds, CountManifolds( shapeId, false ) );
		bool awake = b3Body_IsAwake( bodyId );
		if ( run->sleptAt < 0 && awake == false )
		{
			run->sleptAt = i;
		}
		if ( i >= lateFrom )
		{
			run->maxSpeedLate = b3MaxFloat( run->maxSpeedLate, b3Length( b3Body_GetLinearVelocity( bodyId ) ) );
		}
		run->maxAngularSpeed = b3MaxFloat( run->maxAngularSpeed, b3Length( b3Body_GetAngularVelocity( bodyId ) ) );
	}
	run->meanStep = total / (float)steps;
	b3Pos pos = b3Body_GetPosition( bodyId );
	run->endPosition = ( b3Vec3 ){ (float)pos.x, (float)pos.y, (float)pos.z };
}

static void PrintRun( const char* name, bool points, const Run* r, float expectedY )
{
	printf( "  %-30s %-6s step mean %6.2f max %6.2f (collide %6.2f solve %6.2f) | pairs %d sat %5d samples %6d points %6d "
			"manifolds %2d | min_sep %7.3f | end y %7.3f (rest y %7.3f, sink %6.3f) late speed %6.3f max spin %5.2f slept %d\n",
			name, points ? "points" : "boxes", r->meanStep, r->maxStep, r->maxCollide, r->maxSolve, r->maxPairs, r->maxSat,
			r->maxSamples, r->maxPoints, r->maxManifolds, r->minSeparation == FLT_MAX ? 0.0f : r->minSeparation, r->endPosition.y,
			expectedY, expectedY - r->endPosition.y, r->maxSpeedLate, r->maxAngularSpeed, r->sleptAt );
}

static b3BodyId MakeStaticGround( b3WorldId worldId, PointGrid* ground, int cells )
{
	// Top at y = 0
	float half = (float)cells;
	MakePointGrid( ground, s_groundBoxes, 4, cells, 1, cells, ( b3Vec3 ){ -half, -1.0f, -half } );
	b3BodyDef bodyDef = b3DefaultBodyDef();
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = 0.6f;
	b3CreateVoxelGridShape( bodyId, &shapeDef, ground->grid );
	return bodyId;
}

// The ground as the game builds it: chunk grids of 8 x 1 x 8 cells, each its own static shape with an empty padding ring,
// so the faces at a chunk's edge are exposed. chunks x chunks of them centred on the origin, top at y = 0.
#define MAX_CHUNKS 16
typedef struct ChunkedGround
{
	PointGrid chunks[MAX_CHUNKS];
	int count;
} ChunkedGround;

static void MakeChunkedGround( b3WorldId worldId, ChunkedGround* ground, int chunks )
{
	ground->count = 0;
	float half = 8.0f * (float)chunks;
	for ( int a = 0; a < chunks; ++a )
	{
		for ( int b = 0; b < chunks; ++b )
		{
			PointGrid* g = ground->chunks + ground->count++;
			b3Vec3 origin = { -half + 16.0f * (float)a, -1.0f, -half + 16.0f * (float)b };
			MakePointGrid( g, s_groundBoxes, 4, 8, 1, 8, origin );
			b3BodyDef bodyDef = b3DefaultBodyDef();
			b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
			b3ShapeDef shapeDef = b3DefaultShapeDef();
			shapeDef.baseMaterial.friction = 0.6f;
			b3CreateVoxelGridShape( bodyId, &shapeDef, g->grid );
		}
	}
}

static void DestroyChunkedGround( ChunkedGround* ground )
{
	for ( int i = 0; i < ground->count; ++i )
	{
		DestroyPointGrid( ground->chunks + i );
	}
}

static b3BodyId MakeDynamicGrid( b3WorldId worldId, const PointGrid* g, b3Vec3 position, b3Vec3 velocity, float friction, float density,
								 b3ShapeId* shapeOut )
{
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = position;
	bodyDef.linearVelocity = velocity;
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = friction;
	shapeDef.density = density;
	*shapeOut = b3CreateVoxelGridShape( bodyId, &shapeDef, g->grid );
	return bodyId;
}

// A tower block of storey cells, cellsX x storeys x cellsZ, on the ground: rest from contact (drop 0) or dropped from a height
static int TowerCase( bool points, int cellsX, int storeys, float drop, const char* name )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );
	b3World_EnableVoxelPointContacts( worldId, points );

	PointGrid ground, tower;
	MakeStaticGround( worldId, &ground, 24 );
	float half = (float)cellsX;
	// Body origin at the block's base centre; the grid goes up from y = 0
	ENSURE( MakePointGrid( &tower, s_storeyBoxes, 5, cellsX, storeys, cellsX, ( b3Vec3 ){ -half, 0.0f, -half } ) );
	b3ShapeId shapeId;
	b3BodyId bodyId = MakeDynamicGrid( worldId, &tower, ( b3Vec3 ){ 0.3f, drop, 0.2f }, b3Vec3_zero, 0.6f, 2400.0f * 0.1f, &shapeId );

	Run run;
	StepAndMeasure( worldId, bodyId, shapeId, 600, 300, &run );
	PrintRun( name, points, &run, 0.0f );
	CountManifolds( shapeId, true );

	DestroyPointGrid( &tower );
	b3DestroyWorld( worldId );
	DestroyPointGrid( &ground );
	return 0;
}

// A 2 x 1 x 2 solid block slides at 6 m/s over the ground's seams: it must slow down without tipping or hopping
static int SlideCase( bool points )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );
	b3World_EnableVoxelPointContacts( worldId, points );

	PointGrid ground, block;
	MakeStaticGround( worldId, &ground, 24 );
	ENSURE( MakePointGrid( &block, s_solidBoxes, 1, 2, 1, 2, ( b3Vec3 ){ -2.0f, 0.0f, -2.0f } ) );
	b3ShapeId shapeId;
	b3BodyId bodyId = MakeDynamicGrid( worldId, &block, ( b3Vec3 ){ -10.0f, 0.0f, 0.1f }, ( b3Vec3 ){ 6.0f, 0.0f, 0.0f }, 0.3f, 2400.0f, &shapeId );

	Run run;
	StepAndMeasure( worldId, bodyId, shapeId, 240, 200, &run );
	PrintRun( "slide 6 m/s mu 0.3", points, &run, 0.0f );
	printf( "      slid %.2f m (expected about %.2f), z drift %.3f\n", run.endPosition.x + 10.0f, 36.0f / ( 2.0f * 0.3f * 9.8f ),
			run.endPosition.z - 0.1f );

	DestroyPointGrid( &block );
	b3DestroyWorld( worldId );
	DestroyPointGrid( &ground );
	return 0;
}

// A fast 1 x 1 x 1 cell shard hits the ground at 40 m/s, tilted: it must stay above it
static int FastCase( bool points )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );
	b3World_EnableVoxelPointContacts( worldId, points );

	PointGrid ground, shard;
	MakeStaticGround( worldId, &ground, 24 );
	ENSURE( MakePointGrid( &shard, s_storeyBoxes, 5, 1, 1, 1, ( b3Vec3 ){ -1.0f, -1.0f, -1.0f } ) );
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = ( b3Vec3 ){ 0.0f, 3.0f, 0.0f };
	bodyDef.rotation = b3MakeQuatFromAxisAngle( b3Normalize( ( b3Vec3 ){ 1.0f, 0.0f, 1.0f } ), 0.4f );
	bodyDef.linearVelocity = ( b3Vec3 ){ 2.0f, -40.0f, 0.0f };
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	b3ShapeId shapeId = b3CreateVoxelGridShape( bodyId, &shapeDef, shard.grid );

	Run run;
	StepAndMeasure( worldId, bodyId, shapeId, 300, 240, &run );
	PrintRun( "fast shard 40 m/s tilted", points, &run, 1.0f );
	printf( "      tunnelled: %s\n", run.endPosition.y < 0.0f ? "YES" : "no" );

	DestroyPointGrid( &shard );
	b3DestroyWorld( worldId );
	DestroyPointGrid( &ground );
	return 0;
}

// Two dynamic blocks: a 3 x 2 x 3 block lands on a 4 x 3 x 4 block resting on the ground
static int StackCase( bool points )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );
	b3World_EnableVoxelPointContacts( worldId, points );

	PointGrid ground, lower, upper;
	MakeStaticGround( worldId, &ground, 24 );
	ENSURE( MakePointGrid( &lower, s_storeyBoxes, 5, 4, 3, 4, ( b3Vec3 ){ -4.0f, 0.0f, -4.0f } ) );
	ENSURE( MakePointGrid( &upper, s_storeyBoxes, 5, 3, 2, 3, ( b3Vec3 ){ -3.0f, 0.0f, -3.0f } ) );
	b3ShapeId lowerShape, upperShape;
	MakeDynamicGrid( worldId, &lower, ( b3Vec3 ){ 0.0f, 0.0f, 0.0f }, b3Vec3_zero, 0.6f, 240.0f, &lowerShape );
	b3BodyId upperId = MakeDynamicGrid( worldId, &upper, ( b3Vec3 ){ 0.5f, 8.0f, 0.3f }, b3Vec3_zero, 0.6f, 240.0f, &upperShape );

	Run run;
	StepAndMeasure( worldId, upperId, upperShape, 600, 400, &run );
	PrintRun( "block lands on block", points, &run, 6.0f );

	DestroyPointGrid( &upper );
	DestroyPointGrid( &lower );
	b3DestroyWorld( worldId );
	DestroyPointGrid( &ground );
	return 0;
}

// The tower on chunked ground, straddling the chunk seams, slightly turned, full concrete density
static int SeamTowerCase( bool points, int cellsX, int storeys, float drop, float yaw, const char* name )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );
	b3World_EnableVoxelPointContacts( worldId, points );

	ChunkedGround* ground = (ChunkedGround*)b3Alloc( (int)sizeof( ChunkedGround ) );
	MakeChunkedGround( worldId, ground, 2 );
	PointGrid tower;
	float half = (float)cellsX;
	ENSURE( MakePointGrid( &tower, s_storeyBoxes, 5, cellsX, storeys, cellsX, ( b3Vec3 ){ -half, 0.0f, -half } ) );
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = ( b3Vec3 ){ 0.33f, drop, 0.27f };
	bodyDef.rotation = b3MakeQuatFromAxisAngle( ( b3Vec3 ){ 0.0f, 1.0f, 0.0f }, yaw );
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = 0.6f;
	shapeDef.density = 2400.0f;
	b3ShapeId shapeId = b3CreateVoxelGridShape( bodyId, &shapeDef, tower.grid );

	Run run;
	StepAndMeasure( worldId, bodyId, shapeId, 600, 300, &run );
	PrintRun( name, points, &run, 0.0f );
	CountManifolds( shapeId, true );

	DestroyPointGrid( &tower );
	b3DestroyWorld( worldId );
	DestroyChunkedGround( ground );
	b3Free( ground, (int)sizeof( ChunkedGround ) );
	return 0;
}

static int RunAll( bool points )
{
	ENSURE( SeamTowerCase( points, 8, 20, 0.0f, 0.0f, "seam tower 8x20x8 at rest" ) == 0 );
	ENSURE( SeamTowerCase( points, 8, 20, 0.0f, 0.03f, "seam tower turned at rest" ) == 0 );
	ENSURE( SeamTowerCase( points, 8, 20, 3.0f, 0.03f, "seam tower turned drop 3 m" ) == 0 );
	ENSURE( SeamTowerCase( points, 12, 40, 3.0f, 0.03f, "seam tower 12x40x12 drop 3 m" ) == 0 );
	ENSURE( TowerCase( points, 8, 20, 0.0f, "tower 8x20x8 at rest" ) == 0 );
	ENSURE( TowerCase( points, 8, 20, 2.0f, "tower 8x20x8 drop 2 m" ) == 0 );
	ENSURE( TowerCase( points, 12, 40, 0.0f, "tower 12x40x12 at rest" ) == 0 );
	ENSURE( SlideCase( points ) == 0 );
	ENSURE( FastCase( points ) == 0 );
	ENSURE( StackCase( points ) == 0 );
	return 0;
}

int VoxelPointBench( void )
{
	printf( "box pairs\n" );
	ENSURE( RunAll( false ) == 0 );
	printf( "point samples\n" );
	ENSURE( RunAll( true ) == 0 );
	return 0;
}
