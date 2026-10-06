// SPDX-FileCopyrightText: 2026 Terra Prime
// SPDX-License-Identifier: MIT

#pragma once

#include "core.h"
#include "math_internal.h"

#include "box3d/collision.h"
#include "box3d/types.h"

// A voxel grid is a static shape made of cells. Each cell names a module (or none) and a module is a set of solid
// axis-aligned boxes in the cell's frame. There is no tree and no triangle: queries walk the cells their bounds
// overlap and collide with each box through the convex code on a prebuilt box hull.
//
// A box face that solid voxels cover whole, in its own cell or in a neighbour cell (the padding ring included), cannot
// be touched from outside. A contact that would push through such a face is answered on the box's exposed faces
// instead, so a body sliding over several boxes never snags on the seams between them.
//
// Faces are numbered like the faces of b3BoxHull: 0 -X, 1 +X, 2 -Y, 3 +Y, 4 -Z, 5 +Z.

// A module's boxes and what is derived from them. Immutable once created, shared by any number of grids.
struct b3VoxelGridModule
{
	b3AtomicInt refCount;
	int byteCount;
	int boxCount;
	float cellMeters;
	int cellVoxels;

	// The bounds of all boxes in the cell frame
	b3AABB bounds;

	// boxCount * 6 floats: minimum xyz then maximum xyz
	float* boxes;

	// One prebuilt hull per box, in the cell frame. b3BoxHull is relocatable.
	b3BoxHull* hulls;

	// Per box, the faces its own cell's voxels cover whole, and the faces on the cell boundary (their cover depends on
	// the neighbour cell).
	uint8_t* coveredFaces;
	uint8_t* boundaryFaces;

	// Per side, the solid voxels of the boundary layer: cellVoxels rows of bits, one row per v, bit u. For side
	// axis * 2 + side, u = axis + 1 and v = axis + 2 modulo 3. NULL when the boxes are not on the voxel lattice.
	uint32_t* boundaryLayers;
};

struct b3VoxelGrid
{
	b3AtomicInt refCount;
	float cellMeters;
	int cellVoxels;
	int size[3];
	b3Vec3 origin;

	// Sub-shape keys are ( cell << boxBits ) | box, cell = ( j * sizeZ + k ) * sizeX + i
	int boxBits;
	int cellBits;

	// The whole grid, occupied or not, so a cell update never moves it.
	b3AABB bounds;

	// The bounds of the boxes of the occupied cells, in the grid frame. A grid on a moving body uses it for its AABB, which
	// moves every step and wants to be tight. Kept up to date by b3VoxelGrid_SetCells.
	b3AABB occupiedBounds;

	int paddedCount;
	int* paddedModules;

	int moduleCount;
	int moduleCapacity;
	b3VoxelGridModule** modules;
};

static inline int b3VoxelGridPaddedIndex( const b3VoxelGrid* grid, int i, int j, int k )
{
	return ( ( j + 1 ) * ( grid->size[2] + 2 ) + ( k + 1 ) ) * ( grid->size[0] + 2 ) + ( i + 1 );
}

static inline int b3VoxelGridCellIndex( const b3VoxelGrid* grid, int i, int j, int k )
{
	return ( j * grid->size[2] + k ) * grid->size[0] + i;
}

static inline void b3VoxelGridCellCoordinates( const b3VoxelGrid* grid, int cell, int* i, int* j, int* k )
{
	*i = cell % grid->size[0];
	*k = ( cell / grid->size[0] ) % grid->size[2];
	*j = cell / ( grid->size[0] * grid->size[2] );
}

static inline b3Vec3 b3VoxelGridCellCorner( const b3VoxelGrid* grid, int cell )
{
	int i, j, k;
	b3VoxelGridCellCoordinates( grid, cell, &i, &j, &k );
	return (b3Vec3){ grid->origin.x + (float)i * grid->cellMeters, grid->origin.y + (float)j * grid->cellMeters,
					 grid->origin.z + (float)k * grid->cellMeters };
}

static inline int b3VoxelGridKey( const b3VoxelGrid* grid, int cell, int box )
{
	return ( cell << grid->boxBits ) | box;
}

static inline int b3VoxelGridKeyCell( const b3VoxelGrid* grid, int key )
{
	return key >> grid->boxBits;
}

static inline int b3VoxelGridKeyBox( const b3VoxelGrid* grid, int key )
{
	return key & ( ( 1 << grid->boxBits ) - 1 );
}

// The module of a cell or NULL when it is empty
static inline const b3VoxelGridModule* b3GetVoxelGridCellModule( const b3VoxelGrid* grid, int cell )
{
	int i, j, k;
	b3VoxelGridCellCoordinates( grid, cell, &i, &j, &k );
	int index = grid->paddedModules[b3VoxelGridPaddedIndex( grid, i, j, k )];
	return index >= 0 ? grid->modules[index] : NULL;
}

// A box's bounds in the cell frame
static inline b3AABB b3GetVoxelBoxLocalBounds( const b3VoxelGridModule* module, int box )
{
	const float* b = module->boxes + 6 * box;
	return (b3AABB){ { b[0], b[1], b[2] }, { b[3], b[4], b[5] } };
}

// Visit the boxes whose bounds overlap a box in the grid frame. bounds is the box's bounds in the grid frame. Return false
// to stop. Returns false if stopped.
typedef bool b3VoxelBoxFcn( int cell, int box, b3AABB bounds, void* context );
bool b3QueryVoxelGrid( const b3VoxelGrid* grid, b3AABB bounds, b3VoxelBoxFcn* fcn, void* context );

// The faces of a box that solid voxels cover whole, its own cell's and the neighbour cells'.
uint8_t b3GetVoxelBoxCoveredFaces( const b3VoxelGrid* grid, int cell, int box );

// Does the normal (the box's outward direction toward the other shape) have a component into a covered face?
bool b3IsNormalIntoCoveredFace( uint8_t coveredFaces, b3Vec3 normal );

// Mass of the solid boxes of a grid: exact for a union of disjoint boxes. The inertia is about the returned center.
b3MassData b3ComputeVoxelGridMass( const b3VoxelGrid* grid, float density );

// The AABB of the occupied boxes under a transform
b3AABB b3ComputeVoxelGridOccupiedAABB( const b3VoxelGrid* grid, b3Transform transform );

// Collision extents of the occupied boxes about a local center: half the thinnest side of the occupied bounds and the farthest
// corner on each axis. A grid with no boxes is not a collision hazard (huge minimum extent).
b3ShapeExtent b3ComputeVoxelGridExtent( const b3VoxelGrid* grid, b3Vec3 localCenter );

// Is a box covered on every face, so that nothing can touch it?
bool b3IsVoxelBoxEnclosed( const b3VoxelGrid* grid, int cell, int box );
