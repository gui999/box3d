// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

#include "voxel_grid.h"
#include "aabb.h"

#include "platform.h"
#include "shape.h"

#include "box3d/box3d.h"
#include "box3d/constants.h"

#include <float.h>
#include <math.h>
#include <string.h>

static const int s_faceStep[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };

static size_t b3AlignUp( size_t size )
{
	return ( size + 15 ) & ~(size_t)15;
}

static int b3BitsFor( int count )
{
	int bits = 0;
	while ( ( (int64_t)1 << bits ) < (int64_t)count )
	{
		bits += 1;
	}
	return bits;
}

// Whether a neighbour's boundary layer covers a face's rectangle (u0..u1, v0..v1 in voxels)
static bool b3LayerCovers( const uint32_t* layer, int u0, int u1, int v0, int v1 )
{
	uint32_t bits = ( u1 - u0 >= 32 ? 0xffffffffu : ( ( 1u << ( u1 - u0 ) ) - 1u ) ) << u0;
	for ( int v = v0; v < v1; ++v )
	{
		if ( ( layer[v] & bits ) != bits )
		{
			return false;
		}
	}
	return true;
}

// The lattice bounds of a box in voxels (low xyz then size xyz) if it lies on the voxel lattice inside the cell.
static bool b3GetVoxelBox( const b3VoxelGridModule* module, int box, int out[6] )
{
	float voxel = module->cellMeters / (float)module->cellVoxels;
	int n = module->cellVoxels;
	const float* b = module->boxes + 6 * box;
	for ( int a = 0; a < 3; ++a )
	{
		float low = b[a] / voxel, high = b[3 + a] / voxel;
		int lowVoxel = (int)lroundf( low ), highVoxel = (int)lroundf( high );
		if ( fabsf( low - (float)lowVoxel ) > 1.0e-3f || fabsf( high - (float)highVoxel ) > 1.0e-3f || lowVoxel < 0 ||
			 highVoxel > n || highVoxel <= lowVoxel )
		{
			return false;
		}
		out[a] = lowVoxel;
		out[3 + a] = highVoxel - lowVoxel;
	}
	return true;
}

b3VoxelGridModule* b3CreateVoxelGridModule( const float* boxes, int boxCount, float cellMeters, int cellVoxels )
{
	if ( boxCount < 0 || ( boxCount > 0 && boxes == NULL ) || cellVoxels < 1 || cellVoxels > 32 || b3IsValidFloat( cellMeters ) == false ||
		 cellMeters <= 0.0f )
	{
		return NULL;
	}

	for ( int b = 0; b < boxCount; ++b )
	{
		for ( int a = 0; a < 3; ++a )
		{
			float low = boxes[6 * b + a], high = boxes[6 * b + 3 + a];
			if ( b3IsValidFloat( low ) == false || b3IsValidFloat( high ) == false || high <= low )
			{
				return NULL;
			}
		}
	}

	// One block, sections aligned
	size_t offset = b3AlignUp( sizeof( b3VoxelGridModule ) );
	size_t boxesOffset = offset;
	offset = b3AlignUp( offset + (size_t)boxCount * 6 * sizeof( float ) );
	size_t coveredOffset = offset;
	offset = b3AlignUp( offset + (size_t)boxCount );
	size_t boundaryOffset = offset;
	offset = b3AlignUp( offset + (size_t)boxCount );
	size_t layersOffset = offset;
	offset = b3AlignUp( offset + 6 * (size_t)cellVoxels * sizeof( uint32_t ) );
	size_t hullsOffset = offset;
	offset = b3AlignUp( offset + (size_t)boxCount * sizeof( b3BoxHull ) );

	uint8_t* block = (uint8_t*)b3Alloc( offset );
	memset( block, 0, offset );

	b3VoxelGridModule* module = (b3VoxelGridModule*)block;
	module->refCount.value = 1;
	module->byteCount = (int)offset;
	module->boxCount = boxCount;
	module->cellMeters = cellMeters;
	module->cellVoxels = cellVoxels;
	module->boxes = (float*)( block + boxesOffset );
	module->coveredFaces = block + coveredOffset;
	module->boundaryFaces = block + boundaryOffset;
	module->boundaryLayers = NULL;
	module->hulls = (b3BoxHull*)( block + hullsOffset );

	if ( boxCount > 0 )
	{
		memcpy( module->boxes, boxes, (size_t)boxCount * 6 * sizeof( float ) );
	}

	b3AABB bounds = { { FLT_MAX, FLT_MAX, FLT_MAX }, { -FLT_MAX, -FLT_MAX, -FLT_MAX } };
	for ( int b = 0; b < boxCount; ++b )
	{
		b3AABB box = b3GetVoxelBoxLocalBounds( module, b );
		bounds = b3AABB_Union( bounds, box );

		b3Vec3 center = b3AABB_Center( box );
		b3Vec3 extent = b3AABB_Extents( box );
		module->hulls[b] = b3MakeOffsetBoxHull( extent.x, extent.y, extent.z, center );
	}
	module->bounds = boxCount > 0 ? bounds : (b3AABB){ b3Vec3_zero, b3Vec3_zero };

	// Face cover needs every box on the voxel lattice. A module that is off it covers nothing and nothing covers it.
	int n = cellVoxels;
	int* voxelBoxes = (int*)b3Alloc( (size_t)( boxCount > 0 ? boxCount : 1 ) * 6 * sizeof( int ) );
	bool onLattice = true;
	for ( int b = 0; b < boxCount; ++b )
	{
		if ( b3GetVoxelBox( module, b, voxelBoxes + 6 * b ) == false )
		{
			onLattice = false;
			break;
		}
	}

	if ( onLattice )
	{
		module->boundaryLayers = (uint32_t*)( block + layersOffset );

		size_t solidCount = (size_t)n * n * n;
		uint8_t* solid = (uint8_t*)b3Alloc( solidCount );
		memset( solid, 0, solidCount );

#define B3_SOLID_AT( x, y, z ) ( ( (size_t)( y ) * n + ( z ) ) * n + ( x ) )

		for ( int b = 0; b < boxCount; ++b )
		{
			const int* box = voxelBoxes + 6 * b;
			for ( int y = box[1]; y < box[1] + box[4]; ++y )
			{
				for ( int z = box[2]; z < box[2] + box[5]; ++z )
				{
					for ( int x = box[0]; x < box[0] + box[3]; ++x )
					{
						solid[B3_SOLID_AT( x, y, z )] = 1;
					}
				}
			}
		}

		for ( int axis = 0; axis < 3; ++axis )
		{
			int u = ( axis + 1 ) % 3, v = ( axis + 2 ) % 3;
			for ( int side = 0; side < 2; ++side )
			{
				int layer = side == 0 ? 0 : n - 1;
				uint32_t* rows = module->boundaryLayers + ( axis * 2 + side ) * (size_t)n;
				for ( int j = 0; j < n; ++j )
				{
					for ( int i = 0; i < n; ++i )
					{
						int p[3];
						p[axis] = layer;
						p[u] = i;
						p[v] = j;
						if ( solid[B3_SOLID_AT( p[0], p[1], p[2] )] )
						{
							rows[j] |= 1u << i;
						}
					}
				}
			}
		}

		for ( int b = 0; b < boxCount; ++b )
		{
			const int* box = voxelBoxes + 6 * b;
			int low[3] = { box[0], box[1], box[2] }, size[3] = { box[3], box[4], box[5] };
			for ( int axis = 0; axis < 3; ++axis )
			{
				int u = ( axis + 1 ) % 3, v = ( axis + 2 ) % 3;
				for ( int side = 0; side < 2; ++side )
				{
					int face = axis * 2 + side;
					int layer = side == 0 ? low[axis] - 1 : low[axis] + size[axis];
					if ( layer < 0 || layer >= n )
					{
						module->boundaryFaces[b] |= (uint8_t)( 1 << face );
						continue;
					}

					bool faceCovered = true;
					for ( int j = low[v]; j < low[v] + size[v] && faceCovered; ++j )
					{
						for ( int i = low[u]; i < low[u] + size[u] && faceCovered; ++i )
						{
							int p[3];
							p[axis] = layer;
							p[u] = i;
							p[v] = j;
							faceCovered = solid[B3_SOLID_AT( p[0], p[1], p[2] )] != 0;
						}
					}

					if ( faceCovered )
					{
						module->coveredFaces[b] |= (uint8_t)( 1 << face );
					}
				}
			}
		}

#undef B3_SOLID_AT

		b3Free( solid, solidCount );
	}

	b3Free( voxelBoxes, (size_t)( boxCount > 0 ? boxCount : 1 ) * 6 * sizeof( int ) );
	return module;
}

void b3RetainVoxelGridModule( b3VoxelGridModule* module )
{
	b3AtomicFetchAddInt( &module->refCount, 1 );
}

void b3ReleaseVoxelGridModule( b3VoxelGridModule* module )
{
	if ( module == NULL )
	{
		return;
	}

	int previous = b3AtomicFetchAddInt( &module->refCount, -1 );
	B3_ASSERT( previous >= 1 );
	if ( previous == 1 )
	{
		b3Free( module, (size_t)module->byteCount );
	}
}

int b3VoxelGridModule_GetReferenceCount( const b3VoxelGridModule* module )
{
	return b3AtomicLoadInt( (b3AtomicInt*)&module->refCount );
}

int b3VoxelGridModule_GetBoxCount( const b3VoxelGridModule* module )
{
	return module->boxCount;
}

const float* b3VoxelGridModule_GetBoxes( const b3VoxelGridModule* module )
{
	return module->boxes;
}

// The bounds of the boxes of all occupied cells, inverted (lower above upper) when there are none.
static void b3RecomputeOccupiedBounds( b3VoxelGrid* grid )
{
	b3AABB bounds = { { FLT_MAX, FLT_MAX, FLT_MAX }, { -FLT_MAX, -FLT_MAX, -FLT_MAX } };
	for ( int j = 0; j < grid->size[1]; ++j )
	{
		for ( int k = 0; k < grid->size[2]; ++k )
		{
			for ( int i = 0; i < grid->size[0]; ++i )
			{
				int index = grid->paddedModules[b3VoxelGridPaddedIndex( grid, i, j, k )];
				if ( index < 0 || grid->modules[index]->boxCount == 0 )
				{
					continue;
				}

				const b3VoxelGridModule* module = grid->modules[index];
				b3Vec3 corner = { grid->origin.x + (float)i * grid->cellMeters, grid->origin.y + (float)j * grid->cellMeters,
								  grid->origin.z + (float)k * grid->cellMeters };
				bounds.lowerBound = b3Min( bounds.lowerBound, b3Add( corner, module->bounds.lowerBound ) );
				bounds.upperBound = b3Max( bounds.upperBound, b3Add( corner, module->bounds.upperBound ) );
			}
		}
	}

	// Stays inverted when there is nothing, so that adding a cell does not stretch the bounds to the origin
	grid->occupiedBounds = bounds;
}

b3VoxelGrid* b3CreateVoxelGrid( const b3VoxelGridDef* def )
{
	if ( def->cellCountX < 1 || def->cellCountY < 1 || def->cellCountZ < 1 || def->cellVoxels < 1 || def->cellVoxels > 32 ||
		 b3IsValidFloat( def->cellMeters ) == false || def->cellMeters <= 0.0f || b3IsValidVec3( def->origin ) == false )
	{
		return NULL;
	}

	int64_t cells = (int64_t)def->cellCountX * def->cellCountY * def->cellCountZ;
	int64_t paddedCells = (int64_t)( def->cellCountX + 2 ) * ( def->cellCountY + 2 ) * ( def->cellCountZ + 2 );
	if ( paddedCells > INT32_MAX / 2 )
	{
		return NULL;
	}

	int maxBoxes = def->maxBoxesPerCell > 0 ? def->maxBoxesPerCell : def->cellVoxels * def->cellVoxels * def->cellVoxels;
	int boxBits = b3BitsFor( maxBoxes );
	int cellBits = b3BitsFor( (int)cells );
	if ( cellBits + boxBits > 31 )
	{
		return NULL;
	}

	if ( def->moduleCount > 0 && def->modules == NULL )
	{
		return NULL;
	}

	for ( int i = 0; i < def->moduleCount; ++i )
	{
		if ( def->modules[i] == NULL || def->modules[i]->boxCount > ( 1 << boxBits ) || def->modules[i]->cellVoxels != def->cellVoxels )
		{
			return NULL;
		}
	}

	b3VoxelGrid* grid = (b3VoxelGrid*)b3Alloc( sizeof( b3VoxelGrid ) );
	memset( grid, 0, sizeof( b3VoxelGrid ) );
	grid->refCount.value = 1;
	grid->cellMeters = def->cellMeters;
	grid->cellVoxels = def->cellVoxels;
	grid->size[0] = def->cellCountX;
	grid->size[1] = def->cellCountY;
	grid->size[2] = def->cellCountZ;
	grid->origin = def->origin;
	grid->boxBits = boxBits;
	grid->cellBits = cellBits;
	grid->bounds.lowerBound = def->origin;
	grid->bounds.upperBound = (b3Vec3){ def->origin.x + (float)def->cellCountX * def->cellMeters,
										def->origin.y + (float)def->cellCountY * def->cellMeters,
										def->origin.z + (float)def->cellCountZ * def->cellMeters };

	grid->paddedCount = (int)paddedCells;
	grid->paddedModules = (int*)b3Alloc( (size_t)paddedCells * sizeof( int ) );
	if ( def->paddedCells != NULL )
	{
		memcpy( grid->paddedModules, def->paddedCells, (size_t)paddedCells * sizeof( int ) );
		for ( int i = 0; i < grid->paddedCount; ++i )
		{
			if ( grid->paddedModules[i] < -1 || grid->paddedModules[i] >= def->moduleCount )
			{
				b3Free( grid->paddedModules, (size_t)paddedCells * sizeof( int ) );
				b3Free( grid, sizeof( b3VoxelGrid ) );
				return NULL;
			}
		}
	}
	else
	{
		for ( int i = 0; i < grid->paddedCount; ++i )
		{
			grid->paddedModules[i] = -1;
		}
	}

	grid->moduleCapacity = def->moduleCount > 4 ? def->moduleCount : 4;
	grid->modules = (b3VoxelGridModule**)b3Alloc( (size_t)grid->moduleCapacity * sizeof( b3VoxelGridModule* ) );
	for ( int i = 0; i < def->moduleCount; ++i )
	{
		b3RetainVoxelGridModule( def->modules[i] );
		grid->modules[i] = def->modules[i];
	}
	grid->moduleCount = def->moduleCount;

	b3RecomputeOccupiedBounds( grid );
	return grid;
}

void b3RetainVoxelGrid( b3VoxelGrid* grid )
{
	b3AtomicFetchAddInt( &grid->refCount, 1 );
}

void b3ReleaseVoxelGrid( b3VoxelGrid* grid )
{
	if ( grid == NULL )
	{
		return;
	}

	int previous = b3AtomicFetchAddInt( &grid->refCount, -1 );
	B3_ASSERT( previous >= 1 );
	if ( previous == 1 )
	{
		for ( int i = 0; i < grid->moduleCount; ++i )
		{
			b3ReleaseVoxelGridModule( grid->modules[i] );
		}

		b3Free( grid->modules, (size_t)grid->moduleCapacity * sizeof( b3VoxelGridModule* ) );
		b3Free( grid->paddedModules, (size_t)grid->paddedCount * sizeof( int ) );
		b3Free( grid, sizeof( b3VoxelGrid ) );
	}
}

int b3VoxelGrid_AddModule( b3VoxelGrid* grid, b3VoxelGridModule* module )
{
	if ( module->boxCount > ( 1 << grid->boxBits ) || module->cellVoxels != grid->cellVoxels )
	{
		return -1;
	}

	if ( grid->moduleCount == grid->moduleCapacity )
	{
		int capacity = 2 * grid->moduleCapacity;
		b3VoxelGridModule** modules = (b3VoxelGridModule**)b3Alloc( (size_t)capacity * sizeof( b3VoxelGridModule* ) );
		memcpy( modules, grid->modules, (size_t)grid->moduleCount * sizeof( b3VoxelGridModule* ) );
		b3Free( grid->modules, (size_t)grid->moduleCapacity * sizeof( b3VoxelGridModule* ) );
		grid->modules = modules;
		grid->moduleCapacity = capacity;
	}

	b3RetainVoxelGridModule( module );
	grid->modules[grid->moduleCount] = module;
	grid->moduleCount += 1;
	return grid->moduleCount - 1;
}

void b3VoxelGrid_SetCells( b3VoxelGrid* grid, const int* paddedCells, const int* modules, int count )
{
	// The occupied bounds grow with the cells added and are recomputed when a box that touched them goes
	b3AABB occupied = grid->occupiedBounds;
	bool recompute = false;
	int sx = grid->size[0] + 2, sz = grid->size[2] + 2;
	for ( int n = 0; n < count; ++n )
	{
		int padded = paddedCells[n];
		B3_ASSERT( 0 <= padded && padded < grid->paddedCount );
		B3_ASSERT( -1 <= modules[n] && modules[n] < grid->moduleCount );

		int oldIndex = grid->paddedModules[padded];
		int newIndex = modules[n];
		grid->paddedModules[padded] = newIndex;

		// The padding ring is not part of the shape
		int i = padded % sx - 1;
		int k = ( padded / sx ) % sz - 1;
		int j = padded / ( sx * sz ) - 1;
		if ( i < 0 || i >= grid->size[0] || j < 0 || j >= grid->size[1] || k < 0 || k >= grid->size[2] || oldIndex == newIndex )
		{
			continue;
		}

		b3Vec3 corner = { grid->origin.x + (float)i * grid->cellMeters, grid->origin.y + (float)j * grid->cellMeters,
						  grid->origin.z + (float)k * grid->cellMeters };

		if ( oldIndex >= 0 && grid->modules[oldIndex]->boxCount > 0 && recompute == false )
		{
			const b3AABB* old = &grid->modules[oldIndex]->bounds;
			const float tolerance = 1.0e-4f;
			b3Vec3 lower = b3Add( corner, old->lowerBound );
			b3Vec3 upper = b3Add( corner, old->upperBound );
			recompute = lower.x <= occupied.lowerBound.x + tolerance || lower.y <= occupied.lowerBound.y + tolerance ||
						lower.z <= occupied.lowerBound.z + tolerance || upper.x >= occupied.upperBound.x - tolerance ||
						upper.y >= occupied.upperBound.y - tolerance || upper.z >= occupied.upperBound.z - tolerance;
		}

		if ( newIndex >= 0 && grid->modules[newIndex]->boxCount > 0 )
		{
			const b3AABB* added = &grid->modules[newIndex]->bounds;
			occupied.lowerBound = b3Min( occupied.lowerBound, b3Add( corner, added->lowerBound ) );
			occupied.upperBound = b3Max( occupied.upperBound, b3Add( corner, added->upperBound ) );
		}
	}

	if ( recompute )
	{
		b3RecomputeOccupiedBounds( grid );
	}
	else
	{
		grid->occupiedBounds = occupied;
	}
}

int b3VoxelGrid_GetPaddedCellCount( const b3VoxelGrid* grid )
{
	return grid->paddedCount;
}

int b3VoxelGrid_GetModuleCount( const b3VoxelGrid* grid )
{
	return grid->moduleCount;
}

int b3VoxelGrid_GetBoxCount( const b3VoxelGrid* grid )
{
	int count = 0;
	int cells = grid->size[0] * grid->size[1] * grid->size[2];
	for ( int cell = 0; cell < cells; ++cell )
	{
		const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
		if ( module != NULL )
		{
			count += module->boxCount;
		}
	}
	return count;
}

b3AABB b3VoxelGrid_GetBounds( const b3VoxelGrid* grid )
{
	return grid->bounds;
}

b3AABB b3ComputeVoxelGridAABB( const b3VoxelGrid* grid, b3Transform transform )
{
	return b3AABB_Transform( transform, grid->bounds );
}

static inline bool b3IsOccupiedBoundsEmpty( const b3VoxelGrid* grid )
{
	return grid->occupiedBounds.lowerBound.x > grid->occupiedBounds.upperBound.x;
}

b3AABB b3ComputeVoxelGridOccupiedAABB( const b3VoxelGrid* grid, b3Transform transform )
{
	if ( b3IsOccupiedBoundsEmpty( grid ) )
	{
		b3Vec3 point = b3TransformPoint( transform, grid->origin );
		return (b3AABB){ point, point };
	}

	return b3AABB_Transform( transform, grid->occupiedBounds );
}

b3ShapeExtent b3ComputeVoxelGridExtent( const b3VoxelGrid* grid, b3Vec3 localCenter )
{
	b3ShapeExtent extent = { B3_HUGE, b3Vec3_zero };
	if ( b3IsOccupiedBoundsEmpty( grid ) )
	{
		return extent;
	}

	const b3AABB* bounds = &grid->occupiedBounds;
	b3Vec3 size = b3Sub( bounds->upperBound, bounds->lowerBound );

	// A body this thick cannot pass through anything thinner than its own motion before a contact sees it
	extent.minExtent = 0.5f * b3MinFloat( size.x, b3MinFloat( size.y, size.z ) );

	b3Vec3 farthest = b3FarthestPointOnAABB( *bounds, localCenter );
	extent.maxExtent = b3Abs( b3Sub( farthest, localCenter ) );
	return extent;
}

bool b3IsVoxelBoxEnclosed( const b3VoxelGrid* grid, int cell, int box )
{
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	if ( module->boundaryFaces[box] == 0 && module->coveredFaces[box] != 0x3f )
	{
		// Nothing outside the cell can cover a face that is inside it
		return false;
	}

	return b3GetVoxelBoxCoveredFaces( grid, cell, box ) == 0x3f;
}

// The boxes of a grid are disjoint, so mass, center and inertia are sums over them.
b3MassData b3ComputeVoxelGridMass( const b3VoxelGrid* grid, float density )
{
	b3MassData massData = { 0 };
	int cellCount = grid->size[0] * grid->size[1] * grid->size[2];

	// Pass one: mass and center
	float mass = 0.0f;
	b3Vec3 center = b3Vec3_zero;
	for ( int cell = 0; cell < cellCount; ++cell )
	{
		const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
		if ( module == NULL )
		{
			continue;
		}

		b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );
		for ( int box = 0; box < module->boxCount; ++box )
		{
			const float* b = module->boxes + 6 * box;
			float boxMass = density * ( b[3] - b[0] ) * ( b[4] - b[1] ) * ( b[5] - b[2] );
			b3Vec3 boxCenter = { corner.x + 0.5f * ( b[0] + b[3] ), corner.y + 0.5f * ( b[1] + b[4] ),
								 corner.z + 0.5f * ( b[2] + b[5] ) };
			mass += boxMass;
			center = b3MulAdd( center, boxMass, boxCenter );
		}
	}

	if ( mass <= 0.0f )
	{
		return massData;
	}

	center = b3MulSV( 1.0f / mass, center );

	// Pass two: inertia about the center. Each box about its own center, then shifted.
	b3Matrix3 inertia = b3Mat3_zero;
	for ( int cell = 0; cell < cellCount; ++cell )
	{
		const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
		if ( module == NULL )
		{
			continue;
		}

		b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );
		for ( int box = 0; box < module->boxCount; ++box )
		{
			const float* b = module->boxes + 6 * box;
			float dx = b[3] - b[0], dy = b[4] - b[1], dz = b[5] - b[2];
			float boxMass = density * dx * dy * dz;
			b3Vec3 boxCenter = { corner.x + 0.5f * ( b[0] + b[3] ), corner.y + 0.5f * ( b[1] + b[4] ),
								 corner.z + 0.5f * ( b[2] + b[5] ) };

			float k = boxMass / 12.0f;
			b3Matrix3 own = b3Mat3_zero;
			own.cx.x = k * ( dy * dy + dz * dz );
			own.cy.y = k * ( dx * dx + dz * dz );
			own.cz.z = k * ( dx * dx + dy * dy );

			inertia = b3AddMM( inertia, b3AddMM( own, b3Steiner( boxMass, b3Sub( center, boxCenter ) ) ) );
		}
	}

	massData.mass = mass;
	massData.center = center;
	massData.inertia = inertia;
	return massData;
}

bool b3QueryVoxelGrid( const b3VoxelGrid* grid, b3AABB bounds, b3VoxelBoxFcn* fcn, void* context )
{
	// Clamp before converting to avoid overflowing the int conversion
	int lo[3], hi[3];
	for ( int a = 0; a < 3; ++a )
	{
		float low = ( ( &bounds.lowerBound.x )[a] - ( &grid->origin.x )[a] ) / grid->cellMeters;
		float high = ( ( &bounds.upperBound.x )[a] - ( &grid->origin.x )[a] ) / grid->cellMeters;
		float size = (float)grid->size[a];
		low = b3ClampFloat( low, -1.0f, size );
		high = b3ClampFloat( high, -1.0f, size );
		lo[a] = b3MaxInt( 0, (int)floorf( low ) );
		hi[a] = b3MinInt( grid->size[a] - 1, (int)floorf( high ) );
	}

	for ( int j = lo[1]; j <= hi[1]; ++j )
	{
		for ( int k = lo[2]; k <= hi[2]; ++k )
		{
			for ( int i = lo[0]; i <= hi[0]; ++i )
			{
				int index = grid->paddedModules[b3VoxelGridPaddedIndex( grid, i, j, k )];
				if ( index < 0 )
				{
					continue;
				}

				const b3VoxelGridModule* module = grid->modules[index];
				int cell = b3VoxelGridCellIndex( grid, i, j, k );
				b3Vec3 corner = { grid->origin.x + (float)i * grid->cellMeters, grid->origin.y + (float)j * grid->cellMeters,
								  grid->origin.z + (float)k * grid->cellMeters };

				b3AABB local = { b3Sub( bounds.lowerBound, corner ), b3Sub( bounds.upperBound, corner ) };
				if ( module->boxCount == 0 || b3AABB_Overlaps( module->bounds, local ) == false )
				{
					continue;
				}

				const float* boxes = module->boxes;
				for ( int box = 0; box < module->boxCount; ++box )
				{
					const float* b = boxes + 6 * box;
					if ( b[3] < local.lowerBound.x || b[0] > local.upperBound.x || b[4] < local.lowerBound.y ||
						 b[1] > local.upperBound.y || b[5] < local.lowerBound.z || b[2] > local.upperBound.z )
					{
						continue;
					}

					b3AABB boxBounds = { { corner.x + b[0], corner.y + b[1], corner.z + b[2] },
										 { corner.x + b[3], corner.y + b[4], corner.z + b[5] } };
					if ( fcn( cell, box, boxBounds, context ) == false )
					{
						return false;
					}
				}
			}
		}
	}

	return true;
}

uint8_t b3GetVoxelBoxCoveredFaces( const b3VoxelGrid* grid, int cell, int box )
{
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	B3_ASSERT( module != NULL );
	if ( module->boundaryLayers == NULL )
	{
		return 0;
	}

	uint8_t covered = module->coveredFaces[box];
	uint8_t boundary = module->boundaryFaces[box];
	if ( boundary == 0 )
	{
		return covered;
	}

	int i, j, k;
	b3VoxelGridCellCoordinates( grid, cell, &i, &j, &k );

	int voxelBox[6];
	bool onLattice = b3GetVoxelBox( module, box, voxelBox );
	B3_ASSERT( onLattice );
	B3_UNUSED( onLattice );

	int low[3] = { voxelBox[0], voxelBox[1], voxelBox[2] }, size[3] = { voxelBox[3], voxelBox[4], voxelBox[5] };
	int n = grid->cellVoxels;
	for ( int face = 0; face < 6; ++face )
	{
		if ( ( boundary & ( 1 << face ) ) == 0 )
		{
			continue;
		}

		const int* step = s_faceStep[face];
		int neighbour = grid->paddedModules[b3VoxelGridPaddedIndex( grid, i + step[0], j + step[1], k + step[2] )];
		if ( neighbour < 0 )
		{
			continue;
		}

		const b3VoxelGridModule* other = grid->modules[neighbour];
		if ( other->boundaryLayers == NULL )
		{
			continue;
		}

		int axis = face / 2, u = ( axis + 1 ) % 3, v = ( axis + 2 ) % 3;

		// The neighbour's layer on the opposite side
		if ( b3LayerCovers( other->boundaryLayers + (size_t)( face ^ 1 ) * n, low[u], low[u] + size[u], low[v], low[v] + size[v] ) )
		{
			covered |= (uint8_t)( 1 << face );
		}
	}

	return covered;
}

bool b3IsNormalIntoCoveredFace( uint8_t coveredFaces, b3Vec3 normal )
{
	const float* n = &normal.x;
	for ( int a = 0; a < 3; ++a )
	{
		float c = n[a];
		if ( fabsf( c ) >= 1.0e-4f && ( coveredFaces & ( 1 << ( a * 2 + ( c > 0.0f ? 1 : 0 ) ) ) ) )
		{
			return true;
		}
	}
	return false;
}

// Ray

// A ray against a box. Returns false on a miss. A start inside reports fraction zero and face -1.
static bool b3RayCastVoxelBox( b3Vec3 origin, b3Vec3 direction, b3AABB box, float maxFraction, float* fraction, int* face )
{
	float enter = -FLT_MAX;
	float exit = FLT_MAX;
	int enterFace = -1;
	const float* o = &origin.x;
	const float* d = &direction.x;
	const float* lo = &box.lowerBound.x;
	const float* hi = &box.upperBound.x;
	for ( int a = 0; a < 3; ++a )
	{
		if ( d[a] == 0.0f )
		{
			if ( o[a] < lo[a] || o[a] > hi[a] )
			{
				return false;
			}
			continue;
		}

		float t1 = ( lo[a] - o[a] ) / d[a];
		float t2 = ( hi[a] - o[a] ) / d[a];
		int towardsMin = d[a] > 0.0f;
		if ( t1 > t2 )
		{
			float temp = t1;
			t1 = t2;
			t2 = temp;
		}

		if ( t1 > enter )
		{
			enter = t1;
			// Moving towards +: enters through the minimum face
			enterFace = a * 2 + ( towardsMin ? 0 : 1 );
		}

		if ( t2 < exit )
		{
			exit = t2;
		}
	}

	if ( enter > exit || exit <= 0.0f )
	{
		return false;
	}

	if ( enter < 0.0f )
	{
		// Starts inside
		*fraction = 0.0f;
		*face = -1;
		return true;
	}

	if ( enter > maxFraction )
	{
		return false;
	}

	*fraction = enter;
	*face = enterFace;
	return true;
}

typedef struct b3VoxelRayContext
{
	b3Vec3 origin;
	b3Vec3 direction;
	float best;
	bool found;
	int cell;
	int box;
	int face;
} b3VoxelRayContext;

b3CastOutput b3RayCastVoxelGrid( const b3VoxelGrid* grid, const b3RayCastInput* input )
{
	b3CastOutput output = { 0 };
	output.triangleIndex = B3_NULL_INDEX;
	output.childIndex = B3_NULL_INDEX;

	b3Vec3 origin = input->origin;
	b3Vec3 direction = input->translation;
	float maxFraction = input->maxFraction;

	// Clip the ray to the grid
	float tmin = 0.0f, tmax = maxFraction;
	{
		const float* o = &origin.x;
		const float* d = &direction.x;
		const float* lo = &grid->bounds.lowerBound.x;
		const float* hi = &grid->bounds.upperBound.x;
		for ( int a = 0; a < 3; ++a )
		{
			if ( d[a] == 0.0f )
			{
				if ( o[a] < lo[a] || o[a] > hi[a] )
				{
					return output;
				}
				continue;
			}

			float t1 = ( lo[a] - o[a] ) / d[a];
			float t2 = ( hi[a] - o[a] ) / d[a];
			if ( t1 > t2 )
			{
				float temp = t1;
				t1 = t2;
				t2 = temp;
			}

			tmin = b3MaxFloat( tmin, t1 );
			tmax = b3MinFloat( tmax, t2 );
			if ( tmin > tmax )
			{
				return output;
			}
		}
	}

	b3VoxelRayContext ray = { origin, direction, FLT_MAX, false, 0, 0, 0 };

	b3Vec3 start = b3MulAdd( origin, tmin, direction );
	const float cellMeters = grid->cellMeters;

	// A ray along a cell boundary (no motion on that axis) touches the cells on both sides of it: a box face lying on the
	// boundary belongs to the other cell, so both columns are walked.
	int low[3], high[3];
	for ( int a = 0; a < 3; ++a )
	{
		float c = ( ( &start.x )[a] - ( &grid->origin.x )[a] ) / cellMeters;
		bool along = ( &direction.x )[a] == 0.0f;
		low[a] = b3ClampInt( (int)floorf( along ? c - 1.0e-4f : c ), 0, grid->size[a] - 1 );
		high[a] = b3ClampInt( (int)floorf( along ? c + 1.0e-4f : c ), 0, grid->size[a] - 1 );
	}

	for ( int ci = low[0]; ci <= high[0]; ++ci )
	{
		for ( int cj = low[1]; cj <= high[1]; ++cj )
		{
			for ( int ck = low[2]; ck <= high[2]; ++ck )
			{
				int cell[3] = { ci, cj, ck };
				int step[3];
				float next[3], delta[3];
				for ( int a = 0; a < 3; ++a )
				{
					float d = ( &direction.x )[a];
					float o = ( &origin.x )[a];
					float g = ( &grid->origin.x )[a];
					if ( d > 0.0f )
					{
						step[a] = 1;
						delta[a] = cellMeters / d;
						next[a] = ( g + (float)( cell[a] + 1 ) * cellMeters - o ) / d;
					}
					else if ( d < 0.0f )
					{
						step[a] = -1;
						delta[a] = -cellMeters / d;
						next[a] = ( g + (float)cell[a] * cellMeters - o ) / d;
					}
					else
					{
						step[a] = 0;
						delta[a] = FLT_MAX;
						next[a] = FLT_MAX;
					}
				}

				// A walk crosses at most one boundary per cell on each axis
				int maxSteps = grid->size[0] + grid->size[1] + grid->size[2] + 3;
				for ( int walked = 0; walked < maxSteps; ++walked )
				{
					int a = next[0] < next[1] ? ( next[0] < next[2] ? 0 : 2 ) : ( next[1] < next[2] ? 1 : 2 );
					float exit = b3MinFloat( next[a], tmax );

					int moduleIndex = grid->paddedModules[b3VoxelGridPaddedIndex( grid, cell[0], cell[1], cell[2] )];
					if ( moduleIndex >= 0 )
					{
						const b3VoxelGridModule* module = grid->modules[moduleIndex];
						int cellIndex = b3VoxelGridCellIndex( grid, cell[0], cell[1], cell[2] );
						b3Vec3 corner = b3VoxelGridCellCorner( grid, cellIndex );
						b3Vec3 localOrigin = b3Sub( origin, corner );
						for ( int box = 0; box < module->boxCount; ++box )
						{
							float fraction;
							int face;
							float limit = ray.found ? ray.best : maxFraction;
							if ( b3RayCastVoxelBox( localOrigin, direction, b3GetVoxelBoxLocalBounds( module, box ), limit, &fraction,
													&face ) &&
								 ( ray.found == false || fraction < ray.best ) )
							{
								ray.found = true;
								ray.best = fraction;
								ray.cell = cellIndex;
								ray.box = box;
								ray.face = face;
							}
						}
					}

					if ( ray.found && ray.best <= exit )
					{
						break;
					}

					if ( next[a] > tmax )
					{
						break;
					}

					cell[a] += step[a];
					if ( cell[a] < 0 || cell[a] >= grid->size[a] )
					{
						break;
					}

					next[a] += delta[a];
				}
			}
		}
	}

	if ( ray.found )
	{
		output.hit = true;
		output.fraction = ray.best;
		output.point = b3MulAdd( origin, ray.best, direction );
		output.triangleIndex = ray.cell;
		output.childIndex = ray.box;
		if ( ray.face >= 0 )
		{
			b3Vec3 normal = b3Vec3_zero;
			( &normal.x )[ray.face / 2] = ( ray.face & 1 ) ? 1.0f : -1.0f;
			output.normal = normal;
		}
		else
		{
			output.point = origin;
		}
	}

	return output;
}

// Shape cast and overlap

typedef struct b3VoxelCastContext
{
	const b3VoxelGrid* grid;
	b3ShapeCastInput input;
	b3CastOutput output;
	bool found;
} b3VoxelCastContext;

static bool b3VoxelShapeCastFcn( int cell, int box, b3AABB boxBounds, void* context )
{
	B3_UNUSED( boxBounds );
	b3VoxelCastContext* castContext = context;
	const b3VoxelGrid* grid = castContext->grid;
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );

	// Work in the cell frame, where the hull lives
	b3Vec3 points[B3_MAX_SHAPE_CAST_POINTS];
	int count = b3MinInt( castContext->input.proxy.count, B3_MAX_SHAPE_CAST_POINTS );
	for ( int i = 0; i < count; ++i )
	{
		points[i] = b3Sub( castContext->input.proxy.points[i], corner );
	}

	const b3HullData* hull = &module->hulls[box].base;

	b3ShapeCastPairInput pairInput = { 0 };
	pairInput.proxyA = (b3ShapeProxy){ b3GetHullPoints( hull ), hull->vertexCount, 0.0f };
	pairInput.proxyB = (b3ShapeProxy){ points, count, castContext->input.proxy.radius };
	pairInput.transform = b3Transform_identity;
	pairInput.translationB = castContext->input.translation;
	pairInput.maxFraction = castContext->found ? castContext->output.fraction : castContext->input.maxFraction;
	pairInput.canEncroach = castContext->input.canEncroach;

	b3CastOutput output = b3ShapeCast( &pairInput );
	if ( output.hit == false )
	{
		return true;
	}

	if ( castContext->found && output.fraction >= castContext->output.fraction )
	{
		return true;
	}

	// A hit whose normal passes through a face that solid covers is a hit on the seam between boxes, not on the surface
	uint8_t covered = b3GetVoxelBoxCoveredFaces( grid, cell, box );
	if ( covered != 0 && b3IsNormalIntoCoveredFace( covered, output.normal ) )
	{
		return true;
	}

	output.point = b3Add( output.point, corner );
	output.triangleIndex = cell;
	output.childIndex = box;
	castContext->output = output;
	castContext->found = true;
	return true;
}

b3CastOutput b3ShapeCastVoxelGrid( const b3VoxelGrid* grid, const b3ShapeCastInput* input )
{
	b3AABB bounds = b3MakeAABB( input->proxy.points, input->proxy.count, input->proxy.radius );
	b3Vec3 delta = b3MulSV( input->maxFraction, input->translation );
	b3AABB moved = { b3Add( bounds.lowerBound, delta ), b3Add( bounds.upperBound, delta ) };
	b3AABB swept = b3AABB_Union( bounds, moved );

	b3VoxelCastContext context = { 0 };
	context.grid = grid;
	context.input = *input;
	context.output.triangleIndex = B3_NULL_INDEX;
	context.output.childIndex = B3_NULL_INDEX;

	if ( b3AABB_Overlaps( swept, grid->bounds ) == false )
	{
		return context.output;
	}

	b3QueryVoxelGrid( grid, swept, b3VoxelShapeCastFcn, &context );
	return context.output;
}

typedef struct b3VoxelOverlapContext
{
	const b3VoxelGrid* grid;
	const b3ShapeProxy* proxy;
	bool overlap;
} b3VoxelOverlapContext;

static bool b3VoxelOverlapFcn( int cell, int box, b3AABB boxBounds, void* context )
{
	B3_UNUSED( boxBounds );
	b3VoxelOverlapContext* overlapContext = context;
	const b3VoxelGrid* grid = overlapContext->grid;
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );

	b3Vec3 points[B3_MAX_SHAPE_CAST_POINTS];
	int count = b3MinInt( overlapContext->proxy->count, B3_MAX_SHAPE_CAST_POINTS );
	for ( int i = 0; i < count; ++i )
	{
		points[i] = b3Sub( overlapContext->proxy->points[i], corner );
	}

	const b3HullData* hull = &module->hulls[box].base;

	b3DistanceInput input;
	input.proxyA = (b3ShapeProxy){ b3GetHullPoints( hull ), hull->vertexCount, 0.0f };
	input.proxyB = (b3ShapeProxy){ points, count, overlapContext->proxy->radius };
	input.transform = b3Transform_identity;
	input.useRadii = true;

	b3SimplexCache cache = { 0 };
	b3DistanceOutput output = b3ShapeDistance( &input, &cache, NULL, 0 );
	if ( output.distance < B3_OVERLAP_SLOP )
	{
		overlapContext->overlap = true;
		return false;
	}

	return true;
}

bool b3OverlapVoxelGrid( const b3VoxelGrid* grid, b3Transform shapeTransform, const b3ShapeProxy* proxy )
{
	// The proxy is in world space, bring it into the grid frame
	b3Vec3 localPoints[B3_MAX_SHAPE_CAST_POINTS];
	b3ShapeProxy localProxy = b3MakeLocalProxy( proxy, shapeTransform, localPoints );

	b3AABB bounds = b3ComputeProxyAABB( &localProxy );
	if ( b3AABB_Overlaps( bounds, grid->bounds ) == false )
	{
		return false;
	}

	b3VoxelOverlapContext context = { grid, &localProxy, false };
	b3QueryVoxelGrid( grid, bounds, b3VoxelOverlapFcn, &context );
	return context.overlap;
}

// Recovery of a proxy that sank into the grid

typedef struct b3VoxelRecoverContext
{
	const b3VoxelGrid* grid;
	const b3ShapeProxy* proxy;
	b3AABB query;
	// The deepest signed push along each grid axis
	float push[3];
	int count;
} b3VoxelRecoverContext;

static bool b3VoxelRecoverFcn( int cell, int box, b3AABB boxBounds, void* context )
{
	b3VoxelRecoverContext* recover = context;
	const b3VoxelGrid* grid = recover->grid;
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );

	// The query box is a loose bound of the proxy, so confirm the proxy touches the box itself
	b3Vec3 points[B3_MAX_SHAPE_CAST_POINTS];
	int count = b3MinInt( recover->proxy->count, B3_MAX_SHAPE_CAST_POINTS );
	for ( int i = 0; i < count; ++i )
	{
		points[i] = b3Sub( recover->proxy->points[i], corner );
	}

	const b3HullData* hull = &module->hulls[box].base;
	b3DistanceInput input;
	input.proxyA = (b3ShapeProxy){ b3GetHullPoints( hull ), hull->vertexCount, 0.0f };
	input.proxyB = (b3ShapeProxy){ points, count, recover->proxy->radius };
	input.transform = b3Transform_identity;
	input.useRadii = true;
	b3SimplexCache cache = { 0 };
	if ( b3ShapeDistance( &input, &cache, NULL, 0 ).distance >= B3_OVERLAP_SLOP )
	{
		return true;
	}

	// Only a face that solid does not cover can push the proxy out: the faces on a seam between two boxes cannot, so a proxy
	// across a seam goes straight up or out through the exposed face and never into the seam.
	uint8_t covered = b3GetVoxelBoxCoveredFaces( grid, cell, box );

	const float* lowerQuery = &recover->query.lowerBound.x;
	const float* upperQuery = &recover->query.upperBound.x;
	const float* lowerBox = &boxBounds.lowerBound.x;
	const float* upperBox = &boxBounds.upperBound.x;

	// Boxes that only touch the query, or miss it by the proxy's rounded shape, have nothing to push out of
	for ( int axis = 0; axis < 3; ++axis )
	{
		float depth = b3MinFloat( upperQuery[axis], upperBox[axis] ) - b3MaxFloat( lowerQuery[axis], lowerBox[axis] );
		if ( depth <= 1.0e-5f )
		{
			return true;
		}
	}

	int bestAxis = -1;
	float bestPush = 0.0f;
	for ( int axis = 0; axis < 3; ++axis )
	{
		// Out through the + face moves the proxy up past the box's upper bound, through the - face down past its lower
		if ( ( covered & ( 1 << ( 2 * axis + 1 ) ) ) == 0 )
		{
			float distance = upperBox[axis] - lowerQuery[axis];
			if ( distance > 0.0f && ( bestAxis < 0 || distance < fabsf( bestPush ) ) )
			{
				bestAxis = axis;
				bestPush = distance;
			}
		}

		if ( ( covered & ( 1 << ( 2 * axis ) ) ) == 0 )
		{
			float distance = upperQuery[axis] - lowerBox[axis];
			if ( distance > 0.0f && ( bestAxis < 0 || distance < fabsf( bestPush ) ) )
			{
				bestAxis = axis;
				bestPush = -distance;
			}
		}
	}

	if ( bestAxis < 0 || fabsf( bestPush ) <= 1.0e-5f )
	{
		return true;
	}

	recover->count += 1;
	if ( fabsf( bestPush ) > fabsf( recover->push[bestAxis] ) )
	{
		recover->push[bestAxis] = bestPush;
	}

	return true;
}

int b3RecoverVoxelGrid( const b3VoxelGrid* grid, b3Transform shapeTransform, const b3ShapeProxy* proxy, b3Vec3 pushes[3] )
{
	pushes[0] = pushes[1] = pushes[2] = b3Vec3_zero;

	b3Vec3 localPoints[B3_MAX_SHAPE_CAST_POINTS];
	b3ShapeProxy localProxy = b3MakeLocalProxy( proxy, shapeTransform, localPoints );

	b3AABB bounds = b3ComputeProxyAABB( &localProxy );
	if ( b3AABB_Overlaps( bounds, grid->bounds ) == false )
	{
		return 0;
	}

	b3VoxelRecoverContext context = { grid, &localProxy, bounds, { 0.0f, 0.0f, 0.0f }, 0 };
	b3QueryVoxelGrid( grid, bounds, b3VoxelRecoverFcn, &context );

	// Each push is along a grid axis, answer it in the world frame
	b3Vec3 axes[3] = { b3RotateVector( shapeTransform.q, b3Vec3_axisX ), b3RotateVector( shapeTransform.q, b3Vec3_axisY ),
					   b3RotateVector( shapeTransform.q, b3Vec3_axisZ ) };
	for ( int axis = 0; axis < 3; ++axis )
	{
		pushes[axis] = b3MulSV( context.push[axis], axes[axis] );
	}

	return context.count;
}

// Contacts of a proxy with a grid, one per box

typedef struct b3VoxelContactContext
{
	const b3VoxelGrid* grid;
	const b3ShapeProxy* proxy;
	b3VoxelContact* results;
	int capacity;
	int count;
} b3VoxelContactContext;

static bool b3VoxelContactFcn( int cell, int box, b3AABB boxBounds, void* context )
{
	b3VoxelContactContext* contactContext = context;
	const b3VoxelGrid* grid = contactContext->grid;
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );
	const b3ShapeProxy* proxy = contactContext->proxy;

	b3Vec3 points[B3_MAX_SHAPE_CAST_POINTS];
	int count = b3MinInt( proxy->count, B3_MAX_SHAPE_CAST_POINTS );
	for ( int i = 0; i < count; ++i )
	{
		points[i] = b3Sub( proxy->points[i], corner );
	}

	// The distance between the cores, no radii: the contact normal and depth come from it
	const b3HullData* hull = &module->hulls[box].base;
	b3DistanceInput input;
	input.proxyA = (b3ShapeProxy){ b3GetHullPoints( hull ), hull->vertexCount, 0.0f };
	input.proxyB = (b3ShapeProxy){ points, count, 0.0f };
	input.transform = b3Transform_identity;
	input.useRadii = false;
	b3SimplexCache cache = { 0 };
	b3DistanceOutput output = b3ShapeDistance( &input, &cache, NULL, 0 );
	if ( output.distance >= proxy->radius )
	{
		return true;
	}

	// A face that solid covers is the seam between boxes: contacts through it belong to the neighbour
	uint8_t covered = b3GetVoxelBoxCoveredFaces( grid, cell, box );

	b3Vec3 normal;
	float depth;
	b3Vec3 point;
	if ( output.distance > 1.0e-5f )
	{
		// The core is outside the box. The true closest feature gives the normal: a face, or the diagonal of an edge or corner
		normal = b3MulSV( 1.0f / output.distance, b3Sub( output.pointB, output.pointA ) );
		depth = proxy->radius - output.distance;
		point = output.pointA;
		if ( covered != 0 && b3IsNormalIntoCoveredFace( covered, normal ) )
		{
			return true;
		}
	}
	else
	{
		// The core is inside the box: leave through the nearest exposed face
		b3Vec3 coreLower = points[0];
		b3Vec3 coreUpper = points[0];
		for ( int i = 1; i < count; ++i )
		{
			coreLower = b3Min( coreLower, points[i] );
			coreUpper = b3Max( coreUpper, points[i] );
		}
		const float* lowerCore = &coreLower.x;
		const float* upperCore = &coreUpper.x;
		b3Vec3 lowerBox = b3Sub( boxBounds.lowerBound, corner );
		b3Vec3 upperBox = b3Sub( boxBounds.upperBound, corner );
		const float* lowerBoxF = &lowerBox.x;
		const float* upperBoxF = &upperBox.x;

		int bestAxis = -1;
		float bestPush = 0.0f;
		for ( int axis = 0; axis < 3; ++axis )
		{
			if ( ( covered & ( 1 << ( 2 * axis + 1 ) ) ) == 0 )
			{
				float distance = upperBoxF[axis] - lowerCore[axis] + proxy->radius;
				if ( distance > 0.0f && ( bestAxis < 0 || distance < fabsf( bestPush ) ) )
				{
					bestAxis = axis;
					bestPush = distance;
				}
			}

			if ( ( covered & ( 1 << ( 2 * axis ) ) ) == 0 )
			{
				float distance = upperCore[axis] - lowerBoxF[axis] + proxy->radius;
				if ( distance > 0.0f && ( bestAxis < 0 || distance < fabsf( bestPush ) ) )
				{
					bestAxis = axis;
					bestPush = -distance;
				}
			}
		}

		if ( bestAxis < 0 )
		{
			return true;
		}

		normal = b3Vec3_zero;
		( &normal.x )[bestAxis] = bestPush > 0.0f ? 1.0f : -1.0f;
		depth = fabsf( bestPush );
		// The point of the box nearest the middle of the core
		b3Vec3 center = b3Vec3_zero;
		for ( int i = 0; i < count; ++i )
		{
			center = b3Add( center, b3MulSV( 1.0f / (float)count, points[i] ) );
		}
		point = b3Clamp( center, lowerBox, upperBox );
	}

	if ( depth <= 1.0e-5f )
	{
		return true;
	}

	b3VoxelContact contact = { normal, depth, b3Add( point, corner ) };
	if ( contactContext->count < contactContext->capacity )
	{
		contactContext->results[contactContext->count++] = contact;
	}
	else
	{
		// Keep the deepest
		int shallowest = 0;
		for ( int i = 1; i < contactContext->capacity; ++i )
		{
			if ( contactContext->results[i].depth < contactContext->results[shallowest].depth )
			{
				shallowest = i;
			}
		}
		if ( contact.depth > contactContext->results[shallowest].depth )
		{
			contactContext->results[shallowest] = contact;
		}
	}

	return true;
}

int b3CollideVoxelGrid( const b3VoxelGrid* grid, b3Transform shapeTransform, const b3ShapeProxy* proxy, b3VoxelContact* results, int capacity )
{
	b3Vec3 localPoints[B3_MAX_SHAPE_CAST_POINTS];
	b3ShapeProxy localProxy = b3MakeLocalProxy( proxy, shapeTransform, localPoints );

	b3AABB bounds = b3ComputeProxyAABB( &localProxy );
	if ( capacity <= 0 || b3AABB_Overlaps( bounds, grid->bounds ) == false )
	{
		return 0;
	}

	b3VoxelContactContext context = { grid, &localProxy, results, capacity, 0 };
	b3QueryVoxelGrid( grid, bounds, b3VoxelContactFcn, &context );

	// Answer in the world frame
	for ( int i = 0; i < context.count; ++i )
	{
		results[i].normal = b3RotateVector( shapeTransform.q, results[i].normal );
		results[i].point = b3TransformPoint( shapeTransform, results[i].point );
	}

	return context.count;
}

// Character mover

typedef struct b3VoxelMoverContext
{
	const b3VoxelGrid* grid;
	const b3Capsule* mover;
	b3PlaneResult* results;
	int capacity;
	int count;
} b3VoxelMoverContext;

static bool b3VoxelMoverFcn( int cell, int box, b3AABB boxBounds, void* context )
{
	B3_UNUSED( boxBounds );
	b3VoxelMoverContext* moverContext = context;
	const b3VoxelGrid* grid = moverContext->grid;
	const b3VoxelGridModule* module = b3GetVoxelGridCellModule( grid, cell );
	b3Vec3 corner = b3VoxelGridCellCorner( grid, cell );

	b3Capsule localMover = *moverContext->mover;
	localMover.center1 = b3Sub( localMover.center1, corner );
	localMover.center2 = b3Sub( localMover.center2, corner );

	b3PlaneResult result;
	if ( b3CollideMoverAndHull( &result, &module->hulls[box].base, &localMover ) == 0 )
	{
		return true;
	}

	// A plane through a face that solid covers is the seam between boxes
	uint8_t covered = b3GetVoxelBoxCoveredFaces( grid, cell, box );
	if ( covered != 0 && b3IsNormalIntoCoveredFace( covered, result.plane.normal ) )
	{
		return true;
	}

	result.point = b3Add( result.point, corner );
	moverContext->results[moverContext->count++] = result;
	return moverContext->count < moverContext->capacity;
}

int b3CollideMoverAndVoxelGrid( b3PlaneResult* results, int capacity, const b3VoxelGrid* grid, const b3Capsule* mover )
{
	b3Vec3 centers[2] = { mover->center1, mover->center2 };
	b3AABB bounds = b3MakeAABB( centers, 2, mover->radius + B3_SPECULATIVE_DISTANCE );
	if ( b3AABB_Overlaps( bounds, grid->bounds ) == false )
	{
		return 0;
	}

	b3VoxelMoverContext context = { grid, mover, results, capacity, 0 };
	b3QueryVoxelGrid( grid, bounds, b3VoxelMoverFcn, &context );
	return context.count;
}
