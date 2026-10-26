// Terra: scheduler contention benchmark. Run explicitly with the filter "ContentionBench"; it is not part of the default run.
// Piles 2000 boxes, steps 300 times with N workers while M busy threads compete, per worker priority.

#include "test_macros.h"

#include "box3d/base.h"
#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/math_functions.h"

#include <stdio.h>
#include <stdlib.h>

#if defined( _WIN32 )
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>

static volatile LONG s_stopBusy;

static DWORD WINAPI BusyMain( LPVOID param )
{
	(void)param;
	volatile double x = 1.0;
	while ( s_stopBusy == 0 )
	{
		for ( int i = 0; i < 1000; ++i )
		{
			x = x * 1.0000001 + 0.5;
		}
	}
	return 0;
}

static int CompareFloat( const void* a, const void* b )
{
	float x = *(const float*)a, y = *(const float*)b;
	return ( x > y ) - ( x < y );
}

static int RunConfig( int workers, int priority, int busyCount, float* p50, float* maxMs )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, -10.0f, 0.0f };
	worldDef.workerCount = (uint32_t)workers;
	worldDef.workerPriority = priority;
	worldDef.enableSleep = false;
	b3WorldId worldId = b3CreateWorld( &worldDef );
	ENSURE( b3World_IsValid( worldId ) );

	b3BodyDef groundDef = b3DefaultBodyDef();
	groundDef.position = (b3Pos){ 0.0f, -1.0f, 0.0f };
	b3BodyId groundId = b3CreateBody( worldId, &groundDef );
	b3BoxHull groundBox = b3MakeBoxHull( 40.0f, 1.0f, 40.0f );
	b3ShapeDef groundShape = b3DefaultShapeDef();
	b3CreateHullShape( groundId, &groundShape, &groundBox.base );

	b3BoxHull box = b3MakeCubeHull( 0.5f );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	for ( int i = 0; i < 2000; ++i )
	{
		int x = i % 10, z = ( i / 10 ) % 10, y = i / 100;
		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.type = b3_dynamicBody;
		bodyDef.position = (b3Pos){ ( x - 4.5f ) * 1.05f + 0.3f * ( ( i * 7 ) % 5 - 2 ) * 0.5f, 0.6f + y * 1.05f, ( z - 4.5f ) * 1.05f + 0.3f * ( ( i * 13 ) % 5 - 2 ) * 0.5f };
		b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
		b3CreateHullShape( bodyId, &shapeDef, &box.base );
	}

	HANDLE busy[64];
	s_stopBusy = 0;
	for ( int i = 0; i < busyCount; ++i )
	{
		busy[i] = CreateThread( NULL, 0, BusyMain, NULL, 0, NULL );
	}

	float times[300];
	for ( int i = 0; i < 300; ++i )
	{
		uint64_t t = b3GetTicks();
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		times[i] = b3GetMilliseconds( t );
	}

	s_stopBusy = 1;
	for ( int i = 0; i < busyCount; ++i )
	{
		WaitForSingleObject( busy[i], INFINITE );
		CloseHandle( busy[i] );
	}
	b3DestroyWorld( worldId );

	float worst = 0.0f;
	for ( int i = 0; i < 300; ++i )
	{
		worst = times[i] > worst ? times[i] : worst;
	}
	qsort( times, 300, sizeof( float ), CompareFloat );
	*p50 = times[150];
	*maxMs = worst;
	return 0;
}

int ContentionBench( void )
{
	SYSTEM_INFO info;
	GetSystemInfo( &info );
	int hardwareThreads = (int)info.dwNumberOfProcessors;
	if ( hardwareThreads > 64 )
	{
		hardwareThreads = 64;
	}

	printf( "  hardware threads %d; 2000 boxes, 300 steps, 4 substeps; p50 / max ms\n", hardwareThreads );
	const int workerCounts[4] = { 1, 2, 4, 12 };
	for ( int busyPass = 0; busyPass < 2; ++busyPass )
	{
		int busy = busyPass == 0 ? 0 : hardwareThreads;
		printf( "  busy threads M=%d\n", busy );
		for ( int w = 0; w < 4; ++w )
		{
			for ( int priority = 0; priority < 2; ++priority )
			{
				if ( workerCounts[w] == 1 && priority == 1 )
				{
					continue;
				}
				float p50, maxMs;
				if ( RunConfig( workerCounts[w], priority, busy, &p50, &maxMs ) )
				{
					return 1;
				}
				printf( "    N=%2d %-12s p50 %6.2f  max %6.2f\n", workerCounts[w], priority == 0 ? "normal" : "above_normal", p50, maxMs );
			}
		}
	}
	return 0;
}
#else
int ContentionBench( void )
{
	return 0;
}
#endif
