// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

// Step cost of dynamic voxel grids in a world full of static ones: fast shards and tower-sized blocks landing, a body made
// static or destroyed while resting on many shapes, bodies with many filter joints. The tests print the worst step's profile.

#include "test_macros.h"

#include "core.h"

#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/math_functions.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct Block
{
	b3VoxelGridModule* module;
	b3VoxelGrid* grid;
} Block;

// A block of sizeX x sizeY x sizeZ cells. A cell holds boxX x boxY x boxZ boxes of boxMeters each, on the voxel lattice.
// The grid origin is the body origin.
static bool MakeBlock( Block* block, int sizeX, int sizeY, int sizeZ, float cellMeters, int cellVoxels, int boxX, int boxY,
					   int boxZ, float boxMeters, b3Vec3 origin )
{
	int boxCount = boxX * boxY * boxZ;
	float* boxes = (float*)b3Alloc( boxCount * 6 * (int)sizeof( float ) );
	int n = 0;
	for ( int j = 0; j < boxY; ++j )
	{
		for ( int k = 0; k < boxZ; ++k )
		{
			for ( int i = 0; i < boxX; ++i )
			{
				float* b = boxes + 6 * n;
				b[0] = (float)i * boxMeters;
				b[1] = (float)j * boxMeters;
				b[2] = (float)k * boxMeters;
				b[3] = b[0] + boxMeters;
				b[4] = b[1] + boxMeters;
				b[5] = b[2] + boxMeters;
				n += 1;
			}
		}
	}
	block->module = b3CreateVoxelGridModule( boxes, boxCount, cellMeters, cellVoxels );
	b3Free( boxes, boxCount * 6 * (int)sizeof( float ) );
	if ( block->module == NULL )
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
	def.cellMeters = cellMeters;
	def.cellVoxels = cellVoxels;
	def.maxBoxesPerCell = boxCount;
	def.modules = &block->module;
	def.moduleCount = 1;
	def.paddedCells = padded;
	block->grid = b3CreateVoxelGrid( &def );
	b3Free( padded, paddedCount * (int)sizeof( int ) );
	return block->grid != NULL;
}

static void DestroyBlock( Block* block )
{
	b3ReleaseVoxelGrid( block->grid );
	b3ReleaseVoxelGridModule( block->module );
}

// A floor of rows static grid shapes, each one body, rows along x. Its top is at y = 0, 2 m cells of 4 x 1 x 4 boxes of 0.5 m.
#define FLOOR_ROWS 32
#define FLOOR_CELLS 32

typedef struct Floor
{
	Block rows[FLOOR_ROWS];
	b3BodyId bodies[FLOOR_ROWS];
} Floor;

static bool MakeFloor( b3WorldId worldId, Floor* floor )
{
	for ( int r = 0; r < FLOOR_ROWS; ++r )
	{
		b3Vec3 origin = { -32.0f, 0.0f, -32.0f + 2.0f * (float)r };
		if ( MakeBlock( floor->rows + r, FLOOR_CELLS, 1, 1, 2.0f, 20, 4, 1, 4, 0.5f, origin ) == false )
		{
			return false;
		}

		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.position = ( b3Vec3 ){ 0.0f, -0.5f, 0.0f };
		floor->bodies[r] = b3CreateBody( worldId, &bodyDef );
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		shapeDef.baseMaterial.friction = 0.6f;
		b3CreateVoxelGridShape( floor->bodies[r], &shapeDef, floor->rows[r].grid );
	}
	return true;
}

static void DestroyFloor( Floor* floor )
{
	for ( int r = 0; r < FLOOR_ROWS; ++r )
	{
		DestroyBlock( floor->rows + r );
	}
}

static b3BodyId CreateGridBody( b3WorldId worldId, const Block* block, b3Vec3 position, b3Vec3 velocity, b3Vec3 spin )
{
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = position;
	bodyDef.linearVelocity = velocity;
	bodyDef.angularVelocity = spin;
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );

	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = 0.6f;
	b3CreateVoxelGridShape( bodyId, &shapeDef, block->grid );
	return bodyId;
}

typedef struct Worst
{
	b3Profile profile;
	int stepIndex;
	int contacts;
	float total;
	float maxSolve, maxTransforms, maxCollide;
	int count;
} Worst;

static void Observe( Worst* worst, b3WorldId worldId, int stepIndex )
{
	b3Profile profile = b3World_GetProfile( worldId );
	worst->total += profile.step;
	worst->count += 1;
	worst->maxSolve = b3MaxFloat( worst->maxSolve, profile.solve );
	worst->maxTransforms = b3MaxFloat( worst->maxTransforms, profile.transforms );
	worst->maxCollide = b3MaxFloat( worst->maxCollide, profile.collide );
	if ( profile.step > worst->profile.step )
	{
		worst->profile = profile;
		worst->stepIndex = stepIndex;
		worst->contacts = b3World_GetCounters( worldId ).contactCount;
	}
}

static void Report( const char* name, const Worst* worst )
{
	const b3Profile* p = &worst->profile;
	printf( "  %-34s worst step %6.2f ms at %3d (mean %5.2f) [pairs %.2f collide %.2f solve %.2f: setup %.2f constraints %.2f "
			"transforms %.2f hit %.2f refit %.2f bullets %.2f split %.2f sleep %.2f] contacts %d | max over steps: collide %.2f solve %.2f transforms %.2f\n",
			name, p->step, worst->stepIndex, worst->total / (float)( worst->count > 0 ? worst->count : 1 ), p->pairs, p->collide,
			p->solve, p->solverSetup, p->constraints, p->transforms, p->hitEvents, p->refit, p->bullets, p->splitIslands,
			p->sleepIslands, worst->contacts, worst->maxCollide, worst->maxSolve, worst->maxTransforms );
}

static float StepAll( b3WorldId worldId, int stepCount, Worst* worst, int firstIndex )
{
	for ( int i = 0; i < stepCount; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		Observe( worst, worldId, firstIndex + i );
	}
	return 0.0f;
}

// 1. A fast flat shard of about 300 boxes lands on a floor of many static grid shapes
static int FastShardOnFloor( int shardX, int shardY, int shardZ, float speed, b3Vec3 spin, const char* name )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Floor* floor = (Floor*)b3Alloc( (int)sizeof( Floor ) );
	ENSURE( MakeFloor( worldId, floor ) );

	Block shard;
	ENSURE( MakeBlock( &shard, shardX, shardY, shardZ, 0.5f, 5, 1, 1, 1, 0.5f, b3Vec3_zero ) );
	float x = -0.25f * (float)shardX;
	float z = -0.25f * (float)shardZ;
	b3BodyId bodyId = CreateGridBody( worldId, &shard, ( b3Vec3 ){ x, 0.4f + 0.5f * speed / 60.0f, z }, ( b3Vec3 ){ 0.0f, -speed, 0.0f }, spin );

	Worst worst = { 0 };
	StepAll( worldId, 60, &worst, 0 );
	Report( name, &worst );

	b3Pos p = b3Body_GetPosition( bodyId );
	ENSURE( p.y > -0.2f );

	DestroyBlock( &shard );
	b3DestroyWorld( worldId );
	DestroyFloor( floor );
	b3Free( floor, (int)sizeof( Floor ) );
	return 0;
}

static int StepFastShard( void )
{
	ENSURE( FastShardOnFloor( 10, 6, 5, 30.0f, b3Vec3_zero, "fast shard 300 boxes 30 m/s" ) == 0 );
	ENSURE( FastShardOnFloor( 10, 6, 5, 60.0f, b3Vec3_zero, "fast shard 300 boxes 60 m/s" ) == 0 );
	ENSURE( FastShardOnFloor( 10, 6, 5, 60.0f, ( b3Vec3 ){ 2.0f, 1.0f, 3.0f }, "fast shard 300 boxes 60 m/s spin" ) == 0 );
	ENSURE( FastShardOnFloor( 10, 2, 5, 60.0f, b3Vec3_zero, "fast flat shard 100 boxes 60 m/s" ) == 0 );
	return 0;
}

// 2. A tower block of 10,000 boxes lands, is made static while resting on the floor, or is destroyed
static int BigBlock( int mode, float speed, const char* name )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Floor* floor = (Floor*)b3Alloc( (int)sizeof( Floor ) );
	ENSURE( MakeFloor( worldId, floor ) );

	Block tower;
	ENSURE( MakeBlock( &tower, 20, 25, 20, 0.5f, 5, 1, 1, 1, 0.5f, b3Vec3_zero ) );

	// mode 0: falls fast, 1: lands then is made static, 2: lands then is destroyed
	b3BodyId bodyId = CreateGridBody( worldId, &tower, ( b3Vec3 ){ -5.0f, 0.4f + 0.5f * speed / 60.0f, -5.0f }, ( b3Vec3 ){ 0.0f, -speed, 0.0f }, b3Vec3_zero );

	Worst landing = { 0 };
	if ( mode == 0 )
	{
		StepAll( worldId, 60, &landing, 0 );
	}
	else
	{
		// Until the first contact, as a body is made static or destroyed when it lands
		for ( int i = 0; i < 60; ++i )
		{
			b3World_Step( worldId, 1.0f / 60.0f, 4 );
			Observe( &landing, worldId, i );
			if ( b3World_GetCounters( worldId ).contactCount > 0 )
			{
				break;
			}
		}
	}
	Report( name, &landing );

	if ( mode == 1 || mode == 2 )
	{
		uint64_t ticks = b3GetTicks();
		if ( mode == 1 )
		{
			b3Body_SetType( bodyId, b3_staticBody );
		}
		else
		{
			b3DestroyBody( bodyId );
		}
		float callMs = b3GetMilliseconds( ticks );

		Worst after = { 0 };
		StepAll( worldId, 10, &after, 0 );
		printf( "  %-34s call %6.2f ms\n", mode == 1 ? "  b3Body_SetType(static)" : "  b3DestroyBody", callMs );
		Report( mode == 1 ? "  after SetType" : "  after destroy", &after );
	}

	DestroyBlock( &tower );
	b3DestroyWorld( worldId );
	DestroyFloor( floor );
	b3Free( floor, (int)sizeof( Floor ) );
	return 0;
}

static int StepBigBlock( void )
{
	ENSURE( BigBlock( 0, 30.0f, "tower 10000 boxes 30 m/s" ) == 0 );
	ENSURE( BigBlock( 0, 60.0f, "tower 10000 boxes 60 m/s" ) == 0 );
	ENSURE( BigBlock( 0, 100.0f, "tower 10000 boxes 100 m/s" ) == 0 );
	ENSURE( BigBlock( 1, 60.0f, "tower lands, set static" ) == 0 );
	ENSURE( BigBlock( 2, 60.0f, "tower lands, destroyed" ) == 0 );
	return 0;
}

// 3. Dynamic bodies each with dozens of filter joints
static int StepFilterJoints( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Floor* floor = (Floor*)b3Alloc( (int)sizeof( Floor ) );
	ENSURE( MakeFloor( worldId, floor ) );

	Block shard;
	ENSURE( MakeBlock( &shard, 2, 2, 2, 0.5f, 5, 1, 1, 1, 0.5f, b3Vec3_zero ) );

	enum
	{
		BodyCount = 200,
		JointsPerBody = 40
	};
	b3BodyId bodies[BodyCount];
	for ( int i = 0; i < BodyCount; ++i )
	{
		float x = -20.0f + 2.0f * (float)( i % 20 );
		float z = -10.0f + 2.0f * (float)( i / 20 );
		bodies[i] = CreateGridBody( worldId, &shard, ( b3Vec3 ){ x, 3.0f + 0.01f * (float)i, z }, b3Vec3_zero, b3Vec3_zero );
	}

	uint64_t ticks = b3GetTicks();
	int jointCount = 0;
	for ( int i = 0; i < BodyCount; ++i )
	{
		for ( int n = 1; n <= JointsPerBody / 2; ++n )
		{
			int other = ( i + n ) % BodyCount;
			b3FilterJointDef def = b3DefaultFilterJointDef();
			def.base.bodyIdA = bodies[i];
			def.base.bodyIdB = bodies[other];
			b3CreateFilterJoint( worldId, &def );
			jointCount += 1;
		}
	}
	printf( "  %d filter joints created in %.2f ms\n", jointCount, b3GetMilliseconds( ticks ) );

	Worst worst = { 0 };
	StepAll( worldId, 240, &worst, 0 );
	Report( "200 bodies x 40 filter joints", &worst );

	DestroyBlock( &shard );
	b3DestroyWorld( worldId );
	DestroyFloor( floor );
	b3Free( floor, (int)sizeof( Floor ) );
	return 0;
}

// 4. A tower of floor slabs collapsing: each storey one thin dynamic grid of 400 boxes (10 x 0.5 x 10 m), falling through the
// storeys below onto a floor, with a city of static grid shapes around
#define CITY_SHAPES 900

typedef struct City
{
	Block block;
	b3BodyId bodies[CITY_SHAPES];
} City;

static bool MakeCity( b3WorldId worldId, City* city )
{
	// One module set shared by all, 8 x 6 x 8 cells of 0.5 m boxes: 384 boxes each
	if ( MakeBlock( &city->block, 8, 6, 8, 0.5f, 5, 1, 1, 1, 0.5f, b3Vec3_zero ) == false )
	{
		return false;
	}

	for ( int n = 0; n < CITY_SHAPES; ++n )
	{
		int i = n % 30;
		int k = n / 30;
		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.position = ( b3Vec3 ){ -150.0f + 10.0f * (float)i, 0.0f, -150.0f + 10.0f * (float)k };
		if ( fabsf( bodyDef.position.x ) < 40.0f && fabsf( bodyDef.position.z ) < 40.0f )
		{
			// Keep the middle open for the tower and the floor
			bodyDef.position.x += 100.0f;
		}
		city->bodies[n] = b3CreateBody( worldId, &bodyDef );
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		b3CreateVoxelGridShape( city->bodies[n], &shapeDef, city->block.grid );
	}
	return true;
}

static int CollapseSlabs( int storeys, b3Vec3 spin, bool makeStatic, const char* name )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	Floor* floor = (Floor*)b3Alloc( (int)sizeof( Floor ) );
	ENSURE( MakeFloor( worldId, floor ) );
	City* city = (City*)b3Alloc( (int)sizeof( City ) );
	ENSURE( MakeCity( worldId, city ) );

	Block slab;
	ENSURE( MakeBlock( &slab, 20, 1, 20, 0.5f, 5, 1, 1, 1, 0.5f, b3Vec3_zero ) );

	b3BodyId bodies[64];
	ENSURE( storeys <= 64 );
	for ( int s = 0; s < storeys; ++s )
	{
		b3Vec3 position = { -5.0f + 0.05f * (float)( s % 7 ), 1.0f + 3.0f * (float)s, -5.0f + 0.04f * (float)( s % 5 ) };
		bodies[s] = CreateGridBody( worldId, &slab, position, ( b3Vec3 ){ 0.0f, -2.0f, 0.0f }, spin );
	}

	Worst worst = { 0 };
	float worstSet = 0.0f;
	bool converted[64] = { 0 };
	for ( int i = 0; i < 240; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		Observe( &worst, worldId, i );

		if ( makeStatic )
		{
			// A slab that touches something and moves slowly is made static
			for ( int s = 0; s < storeys; ++s )
			{
				if ( converted[s] == false && b3Body_GetContactCapacity( bodies[s] ) > 0 )
				{
					uint64_t ticks = b3GetTicks();
					b3Body_SetType( bodies[s], b3_staticBody );
					worstSet = b3MaxFloat( worstSet, b3GetMilliseconds( ticks ) );
					converted[s] = true;
				}
			}
		}
	}
	Report( name, &worst );
	if ( makeStatic )
	{
		printf( "    worst b3Body_SetType call %.2f ms\n", worstSet );
	}

	DestroyBlock( &slab );
	b3DestroyWorld( worldId );
	DestroyBlock( &city->block );
	b3Free( city, (int)sizeof( City ) );
	DestroyFloor( floor );
	b3Free( floor, (int)sizeof( Floor ) );
	return 0;
}

static int StepCollapse( void )
{
	ENSURE( CollapseSlabs( 1, b3Vec3_zero, false, "1 slab falls from 1 m" ) == 0 );
	ENSURE( CollapseSlabs( 10, b3Vec3_zero, false, "10 slabs collapse" ) == 0 );
	ENSURE( CollapseSlabs( 40, b3Vec3_zero, false, "40 slabs collapse" ) == 0 );
	ENSURE( CollapseSlabs( 40, ( b3Vec3 ){ 0.3f, 0.2f, 0.5f }, false, "40 slabs collapse, spinning" ) == 0 );
	ENSURE( CollapseSlabs( 40, ( b3Vec3 ){ 0.3f, 0.2f, 0.5f }, true, "40 slabs collapse, made static" ) == 0 );
	return 0;
}

int VoxelGridCostTest( void )
{
	RUN_SUBTEST( StepFastShard );
	RUN_SUBTEST( StepBigBlock );
	RUN_SUBTEST( StepFilterJoints );
	RUN_SUBTEST( StepCollapse );
	return 0;
}
