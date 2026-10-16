// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

// Public-API regressions for exact terrain-owner padding. Deliberately separate
// from VoxelPointBench: no large tower, timing screen, or diagnostic dump.
#include "test_macros.h"
#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/math_functions.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

enum { HALO_CONTACT_CAPACITY = 32, HALO_LANDING_TICKS = 360 };
static const float haloStep = 1.0f / 60.0f;

typedef struct HaloCells
{
	int x, z, sizeX, sizeZ;
	uint32_t occupied;
} HaloCells;

typedef struct HaloFixture
{
	b3WorldId world;
	b3BodyId owners[4], body, fragment;
	b3ShapeId ground[4], fragmentShape;
} HaloFixture;

typedef struct HaloContacts
{
	int points, sides, persisted;
	float horizontal, positiveXImpulse, negativeXImpulse;
	bool fragmentHit, overflow;
} HaloContacts;

static bool HaloSolid( HaloCells cells, int x, int y, int z )
{
	int i = x - cells.x, k = z - cells.z;
	return y == -1 && i >= 0 && i < cells.sizeX && k >= 0 && k < cells.sizeZ &&
		( ( cells.occupied >> ( i + cells.sizeX * k ) ) & 1u ) != 0;
}

static HaloFixture HaloWorld( float gravity )
{
	HaloFixture fixture = { 0 };
	b3WorldDef def = b3DefaultWorldDef();
	def.workerCount = 1;
	def.gravity = (b3Vec3){ 0, gravity, 0 };
	fixture.world = b3CreateWorld( &def );
	b3World_EnableVoxelPointContacts( fixture.world, true );
	return fixture;
}

// Same proxy/occupancy split as the game plugin. Ring cells never enlarge the
// proxy geometry, occupied bounds, body mass, or collidable ownership.
static b3VoxelGrid* HaloGrid( int sizeX, int sizeZ, b3Vec3 origin, const int* padded,
	const float* moduleBox, bool merged )
{
	b3VoxelGridModule* module = b3CreateVoxelGridModule( moduleBox, 1, 2, 20 );
	if ( module == NULL ) return NULL;
	b3VoxelGridDef def = { 0 };
	def.cellCountX = sizeX; def.cellCountY = 1; def.cellCountZ = sizeZ;
	def.origin = origin; def.cellMeters = 2; def.cellVoxels = 20;
	def.maxBoxesPerCell = 1; def.modules = &module; def.moduleCount = 1; def.paddedCells = padded;
	b3VoxelGrid* grid = b3CreateVoxelGrid( &def );
	b3ReleaseVoxelGridModule( module );
	if ( grid == NULL || merged == false ) return grid;
	float extent = 2 * (float)( sizeX > sizeZ ? sizeX : sizeZ );
	const float box[6] = { 0, 0, 0, 2 * (float)sizeX, 2, 2 * (float)sizeZ };
	module = b3CreateVoxelGridModule( box, 1, extent, 20 );
	if ( module == NULL ) { b3ReleaseVoxelGrid( grid ); return NULL; }
	const uint8_t uncovered = 0;
	b3VoxelGridModule_SetCoveredFaces( module, &uncovered, 1 );
	int one[27];
	for ( int i = 0; i < 27; ++i ) one[i] = -1;
	one[13] = 0;
	def.cellCountX = def.cellCountZ = 1; def.cellMeters = extent;
	def.modules = &module; def.paddedCells = one; def.occupancy = grid;
	b3VoxelGrid* proxy = b3CreateVoxelGrid( &def );
	b3ReleaseVoxelGridModule( module );
	b3ReleaseVoxelGrid( grid );
	return proxy;
}

static bool HaloGround( HaloFixture* f, int slot, HaloCells cells, int x, int z, int sizeX, int sizeZ,
	bool halo, bool merged )
{
	int count = ( sizeX + 2 ) * 3 * ( sizeZ + 2 );
	int* padded = malloc( (size_t)count * sizeof( int ) );
	if ( padded == NULL ) return false;
	for ( int j = -1; j <= 1; ++j ) for ( int k = -1; k <= sizeZ; ++k ) for ( int i = -1; i <= sizeX; ++i )
	{
		bool inside = i >= 0 && i < sizeX && j == 0 && k >= 0 && k < sizeZ;
		padded[( ( j + 1 ) * ( sizeZ + 2 ) + k + 1 ) * ( sizeX + 2 ) + i + 1] =
			( inside || halo ) && HaloSolid( cells, x + i, -1 + j, z + k ) ? 0 : -1;
	}
	const float box[6] = { 0, 0, 0, 2, 2, 2 };
	b3VoxelGrid* grid = HaloGrid( sizeX, sizeZ, (b3Vec3){ 2 * (float)x, -2, 2 * (float)z }, padded, box, merged );
	free( padded );
	if ( grid == NULL ) return false;
	if ( B3_IS_NULL( f->owners[slot] ) )
	{
		b3BodyDef bodyDef = b3DefaultBodyDef();
		f->owners[slot] = b3CreateBody( f->world, &bodyDef );
	}
	// Full merged-proxy replacement also discards cached contacts on an edit.
	if ( B3_IS_NON_NULL( f->ground[slot] ) ) b3DestroyShape( f->ground[slot], true );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.baseMaterial.friction = 1; shapeDef.density = 0;
	f->ground[slot] = b3CreateVoxelGridShape( f->owners[slot], &shapeDef, grid );
	b3ReleaseVoxelGrid( grid );
	return B3_IS_NON_NULL( f->ground[slot] );
}

static b3BodyId HaloBody( HaloFixture* f, b3Pos position, b3Vec3 velocity, bool fragment, b3ShapeId* shape )
{
	b3BodyDef def = b3DefaultBodyDef();
	def.type = b3_dynamicBody; def.position = position; def.linearVelocity = velocity;
	def.linearDamping = 0; def.angularDamping = 0.1f;
	b3BodyId body = b3CreateBody( f->world, &def );
	int padded[27];
	for ( int i = 0; i < 27; ++i ) padded[i] = -1;
	padded[13] = 0;
	const float cube[6] = { 0, 0, 0, 2, 2, 2 }, chip[6] = { 0.8f, 0, 0.8f, 1.2f, 0.2f, 1.2f };
	// One-cell body geometry already is its exact collision box; no second proxy is needed.
	b3VoxelGrid* grid = HaloGrid( 1, 1, (b3Vec3){ -1, fragment ? -0.1f : -1, -1 }, padded, fragment ? chip : cube, false );
	if ( grid == NULL ) return b3_nullBodyId;
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 1900; shapeDef.baseMaterial.friction = 0.8f;
	*shape = b3CreateVoxelGridShape( body, &shapeDef, grid );
	b3ReleaseVoxelGrid( grid );
	if ( B3_IS_NULL( *shape ) ) return b3_nullBodyId;
	float mass = fragment ? 60.8f : 15200, width = fragment ? 0.4f : 2, height = fragment ? 0.2f : 2;
	b3MassData data = { 0 };
	data.mass = mass;
	data.inertia.cx.x = data.inertia.cz.z = mass * ( height * height + width * width ) / 12;
	data.inertia.cy.y = mass * 2 * width * width / 12;
	b3Body_SetMassData( body, data );
	return body;
}

static float HaloMinimumY( b3BodyId body )
{
	b3WorldTransform transform = b3Body_GetTransform( body );
	float minimum = FLT_MAX;
	for ( int x = -1; x <= 1; x += 2 ) for ( int y = -1; y <= 1; y += 2 ) for ( int z = -1; z <= 1; z += 2 )
	{
		b3Vec3 corner = b3RotateVector( transform.q, (b3Vec3){ (float)x, (float)y, (float)z } );
		minimum = fminf( minimum, (float)transform.p.y + corner.y );
	}
	return minimum;
}

static HaloContacts HaloRead( const HaloFixture* f )
{
	HaloContacts result = { 0 };
	if ( b3Body_GetContactCapacity( f->body ) > HALO_CONTACT_CAPACITY ) { result.overflow = true; return result; }
	b3ContactData pairs[HALO_CONTACT_CAPACITY];
	int count = b3Body_GetContactData( f->body, pairs, HALO_CONTACT_CAPACITY );
	for ( int c = 0; c < count; ++c )
	{
		bool isA = B3_ID_EQUALS( b3Shape_GetBody( pairs[c].shapeIdA ), f->body );
		b3ShapeId other = isA ? pairs[c].shapeIdB : pairs[c].shapeIdA;
		bool ground = false;
		for ( int s = 0; s < 4; ++s ) if ( B3_IS_NON_NULL( f->ground[s] ) && B3_ID_EQUALS( other, f->ground[s] ) ) ground = true;
		for ( int m = 0; m < pairs[c].manifoldCount; ++m )
		{
			const b3Manifold* manifold = pairs[c].manifolds + m;
			b3Vec3 normal = isA ? manifold->normal : b3Neg( manifold->normal );
			for ( int p = 0; p < manifold->pointCount; ++p )
			{
				const b3ManifoldPoint* point = manifold->points + p;
				if ( B3_IS_NON_NULL( f->fragmentShape ) && B3_ID_EQUALS( other, f->fragmentShape ) && point->appliedNormalImpulse > 0 ) result.fragmentHit = true;
				if ( ground == false ) continue;
				result.points++; result.persisted += point->persisted ? 1 : 0;
				result.horizontal = fmaxf( result.horizontal, sqrtf( normal.x * normal.x + normal.z * normal.z ) );
				result.sides += fabsf( normal.y ) < 0.5f ? 1 : 0;
				if ( normal.x > 0.95f ) result.positiveXImpulse += point->appliedNormalImpulse;
				if ( normal.x < -0.95f ) result.negativeXImpulse += point->appliedNormalImpulse;
			}
		}
	}
	return result;
}

static bool HaloLandingCase( int arm, bool falling )
{
	HaloFixture f = HaloWorld( falling ? -9.8f : 0 );
	HaloCells floor = { -2, -2, 6, 5, ( 1u << 30 ) - 1 };
	bool okay;
	if ( arm < 2 ) okay = HaloGround( &f, 0, floor, -2, -2, 2, 2, arm == 1, true ) &&
		HaloGround( &f, 1, floor, 0, -2, 4, 2, arm == 1, true ) && HaloGround( &f, 2, floor, -2, 0, 2, 3, arm == 1, true ) &&
		HaloGround( &f, 3, floor, 0, 0, 4, 3, arm == 1, true );
	else if ( arm == 2 ) okay = HaloGround( &f, 0, floor, -2, -2, 6, 5, false, true );
	else
	{
		b3BodyDef bodyDef = b3DefaultBodyDef(); bodyDef.position = (b3Pos){ 2, -1, 1 };
		f.owners[0] = b3CreateBody( f.world, &bodyDef );
		b3ShapeDef shapeDef = b3DefaultShapeDef(); shapeDef.baseMaterial.friction = 1;
		b3BoxHull box = b3MakeBoxHull( 6, 1, 5 );
		f.ground[0] = b3CreateHullShape( f.owners[0], &shapeDef, &box.base ); okay = B3_IS_NON_NULL( f.ground[0] );
	}
	if ( falling )
	{
		f.fragment = HaloBody( &f, (b3Pos){ 1, 0.15f, 1 }, b3Vec3_zero, true, &f.fragmentShape );
		okay = okay && B3_IS_NON_NULL( f.fragment );
		for ( int i = 0; i < 6; ++i ) b3World_Step( f.world, haloStep, 4 );
	}
	b3ShapeId shape;
	f.body = HaloBody( &f, (b3Pos){ 1, falling ? 2.5f : 0.98f, 1 }, (b3Vec3){ 0, falling ? -8 : 0, 0 }, false, &shape );
	okay = okay && B3_IS_NON_NULL( f.body );
	int points = 0, sides = 0; float horizontal = 0;
	if ( okay ) for ( int tick = 0; tick < ( falling ? HALO_LANDING_TICKS : 1 ); ++tick )
	{
		b3World_Step( f.world, haloStep, 4 );
		HaloContacts contacts = HaloRead( &f ); okay = okay && !contacts.overflow;
		points += contacts.points; sides += contacts.sides; horizontal = fmaxf( horizontal, contacts.horizontal );
		if ( contacts.fragmentHit )
		{
			b3DestroyBody( f.fragment ); f.fragment = b3_nullBodyId; f.fragmentShape = b3_nullShapeId;
		}
	}
	float underside = okay ? HaloMinimumY( f.body ) : FLT_MAX;
	bool valid = okay && points > 0 && fabsf( underside ) < 0.1f && ( falling ? sides == 0 : horizontal < 0.00001f ) &&
		( !falling || ( !b3Body_IsAwake( f.body ) && B3_IS_NULL( f.fragment ) ) );
	// Empty padding is the deliberately incorrect caller. Preserve its reproducer
	// rather than accepting a test that no longer distinguishes missing neighbours.
	bool passed = arm == 0 ? okay && !valid && sides > 0 && ( !falling || fabsf( underside ) >= 0.1f ) : valid;
	printf( "  terrain halo landing arm=%d fall=%d sides=%d minY=%.6f expectedFailure=%d passed=%d\n", arm, falling, sides, underside, arm == 0, passed );
	b3DestroyWorld( f.world );
	return passed;
}

static int HaloLanding( void )
{
	for ( int arm = 0; arm < 4; ++arm ) for ( int falling = 0; falling < 2; ++falling ) ENSURE( HaloLandingCase( arm, falling != 0 ) );
	return 0;
}

static int HaloExternalWall( void )
{
	HaloFixture f = HaloWorld( 0 ); HaloCells cells = { 0, 0, 1, 1, 1 };
	bool okay = HaloGround( &f, 0, cells, 0, 0, 1, 1, true, true );
	b3ShapeId shape;
	f.body = HaloBody( &f, (b3Pos){ -0.98f, -1, 1 }, (b3Vec3){ 2, 0, 0 }, false, &shape );
	okay = okay && B3_IS_NON_NULL( f.body );
	if ( okay )
	{
		b3World_Step( f.world, haloStep, 4 ); HaloContacts contact = HaloRead( &f );
		okay = !contact.overflow && contact.points > 0 && contact.positiveXImpulse > 1 &&
			fabsf( contact.positiveXImpulse - 30400 ) < 0.01f && (float)b3Body_GetPosition( f.body ).x + 1 < 0.1f;
		printf( "  terrain external wall impulse=%.6f passed=%d\n", contact.positiveXImpulse, okay );
	}
	b3DestroyWorld( f.world ); ENSURE( okay ); return 0;
}

static bool HaloCornerArm( HaloCells cells, bool partitioned, b3Pos* position, float* minimum )
{
	HaloFixture f = HaloWorld( 0 ); bool okay = true;
	if ( partitioned )
	{
		for ( int s = 0; s < 4; ++s ) if ( ( cells.occupied >> s ) & 1u )
			okay = HaloGround( &f, s, cells, cells.x + s % 2, cells.z + s / 2, 1, 1, true, true ) && okay;
	}
	else okay = HaloGround( &f, 0, cells, cells.x, cells.z, 2, 2, false, false );
	b3ShapeId shape;
	f.body = HaloBody( &f, (b3Pos){ 2 * (float)cells.x + 1, 0.98f, 2 * (float)cells.z + 1 }, b3Vec3_zero, false, &shape );
	okay = okay && B3_IS_NON_NULL( f.body );
	if ( okay )
	{
		b3World_Step( f.world, haloStep, 4 ); HaloContacts contacts = HaloRead( &f );
		*position = b3Body_GetPosition( f.body ); *minimum = HaloMinimumY( f.body );
		okay = !contacts.overflow && contacts.points > 0 && contacts.horizontal < 0.00001f && fabsf( *minimum ) < 0.1f;
		b3DestroyBody( f.body ); f.body = b3_nullBodyId;
		for ( int s = 0; s < 4; ++s )
		{
			b3RayResult ray = b3World_CastRayClosest( f.world, (b3Pos){ 2 * (float)( cells.x + s % 2 ) + 1, 1,
				2 * (float)( cells.z + s / 2 ) + 1 }, (b3Vec3){ 0, -4, 0 }, b3DefaultQueryFilter() );
			okay = ( ray.hit == ( ( ( cells.occupied >> s ) & 1u ) != 0 ) && ( !ray.hit || fabsf( (float)ray.point.y ) < 0.0001f ) ) && okay;
		}
		int hole = cells.occupied == 7 ? 3 : 1;
		f.body = HaloBody( &f, (b3Pos){ 2 * (float)( cells.x + hole % 2 ) + 1, 0.5f,
			2 * (float)( cells.z + hole / 2 ) + 1 }, b3Vec3_zero, true, &shape );
		if ( B3_IS_NON_NULL( f.body ) )
		{
			b3World_SetGravity( f.world, (b3Vec3){ 0, -9.8f, 0 } ); int blocked = 0;
			for ( int tick = 0; tick < 60; ++tick ) { b3World_Step( f.world, haloStep, 4 ); HaloContacts c = HaloRead( &f ); blocked += c.points; okay = !c.overflow && okay; }
			okay = blocked == 0 && b3Body_GetPosition( f.body ).y < -2.4f && okay;
		}
		else okay = false;
	}
	b3DestroyWorld( f.world ); return okay;
}

static int HaloCornersAndHoles( void )
{
	const HaloCells cases[2] = { { 0, 0, 2, 2, 7 }, { -2, -2, 2, 2, 9 } };
	for ( int i = 0; i < 2; ++i )
	{
		b3Pos a = { 0 }, b = { 0 }; float minimumA = 0, minimumB = 0;
		bool first = HaloCornerArm( cases[i], true, &a, &minimumA ), second = HaloCornerArm( cases[i], false, &b, &minimumB );
		float delta = b3Length( (b3Vec3){ (float)( a.x - b.x ), (float)( a.y - b.y ), (float)( a.z - b.z ) } );
		printf( "  terrain halo corner case=%d partition=%d union=%d positionDelta=%.6f minYDelta=%.6f\n", i, first, second, delta, fabsf( minimumA - minimumB ) );
		ENSURE( first && second && delta < 0.001f && fabsf( minimumA - minimumB ) < 0.001f );
	}
	return 0;
}

static void HaloReset( HaloFixture* f, b3Pos position, b3Vec3 velocity )
{
	b3Body_SetTransform( f->body, position, b3Quat_identity ); b3Body_SetLinearVelocity( f->body, velocity );
	b3Body_SetAngularVelocity( f->body, b3Vec3_zero ); b3Body_SetAwake( f->body, true );
}

static int HaloWarmReplacement( void )
{
	HaloFixture f = HaloWorld( -9.8f ); HaloCells cells = { 0, 0, 2, 2, 3 };
	bool okay = HaloGround( &f, 0, cells, 0, 0, 1, 1, true, true ) && HaloGround( &f, 1, cells, 1, 0, 1, 1, true, true );
	b3ShapeId shape; f.body = HaloBody( &f, (b3Pos){ 1, 0.98f, 1 }, b3Vec3_zero, false, &shape );
	okay = okay && B3_IS_NON_NULL( f.body );
	if ( okay )
	{
		for ( int phase = 0; phase < 3; ++phase )
		{
			if ( phase == 1 )
			{
				b3DestroyBody( f.owners[1] ); f.owners[1] = b3_nullBodyId; f.ground[1] = b3_nullShapeId; cells.occupied = 1;
				okay = HaloGround( &f, 0, cells, 0, 0, 1, 1, true, true ) && okay;
				b3World_SetGravity( f.world, b3Vec3_zero ); HaloReset( &f, (b3Pos){ 2.98f, -1, 1 }, (b3Vec3){ -2, 0, 0 } );
			}
			if ( phase == 2 )
			{
				cells.occupied = 3;
				okay = HaloGround( &f, 0, cells, 0, 0, 1, 1, true, true ) && okay;
				okay = HaloGround( &f, 1, cells, 1, 0, 1, 1, true, true ) && okay;
				b3World_SetGravity( f.world, (b3Vec3){ 0, -9.8f, 0 } ); HaloReset( &f, (b3Pos){ 1, 0.98f, 1 }, b3Vec3_zero );
			}
			int points = 0, sides = 0, persisted = 0; float impulse = 0, penetration = 0;
			for ( int tick = 0; tick < ( phase == 1 ? 30 : 50 ); ++tick )
			{
				if ( phase == 1 ) b3Body_SetLinearVelocity( f.body, (b3Vec3){ -2, 0, 0 } );
				b3World_Step( f.world, haloStep, 4 ); HaloContacts c = HaloRead( &f ); okay = !c.overflow && okay;
				points += c.points; sides += c.sides; persisted += c.persisted; impulse += c.negativeXImpulse;
				if ( phase == 1 ) penetration = fmaxf( penetration, 3 - (float)b3Body_GetPosition( f.body ).x );
			}
			bool passed = points > 0 && persisted > 0 && ( phase == 1 ? impulse > 1 && penetration < 0.1f : sides == 0 && fabsf( HaloMinimumY( f.body ) ) < 0.1f );
			printf( "  terrain halo replacement phase=%d points=%d sides=%d persisted=%d wallImpulse=%.3f penetration=%.6f passed=%d\n", phase, points, sides, persisted, impulse, penetration, passed );
			okay = passed && okay;
		}
	}
	b3DestroyWorld( f.world ); ENSURE( okay ); return 0;
}

int VoxelTerrainHaloTest( void )
{
	RUN_SUBTEST( HaloLanding );
	RUN_SUBTEST( HaloExternalWall );
	RUN_SUBTEST( HaloCornersAndHoles );
	RUN_SUBTEST( HaloWarmReplacement );
	return 0;
}
