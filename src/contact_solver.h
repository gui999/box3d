// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "math_internal.h"
#include "solver.h"

typedef struct b3ManifoldConstraintPoint
{
	b3Vec3 rA, rB;
	float baseSeparation;
	float relativeVelocity;
	float normalImpulse;
	float totalNormalImpulse;
	// Net impulse applied over the step: each sub-step's accumulated impulse after relax, plus restitution
	float appliedNormalImpulse;
	float normalMass;
	float leverArm;
	// Terra: positive restitution demand awaiting the shared allowance clamp.
	float restitutionImpulse;
} b3ManifoldConstraintPoint;

typedef struct b3ManifoldConstraint
{
	b3ManifoldConstraintPoint points[4];
	int pointCount;
	b3Vec3 normal;
	b3Vec3 tangent1;
	b3Vec3 tangent2;
	// Friction centers
	b3Vec3 centerA, centerB;
	float twistMass;
	float twistImpulse;
	b3Matrix2 tangentMass;
	b3Vec2 frictionImpulse;
	b3Vec3 rollingImpulse;
	float tangentVelocity1;
	float tangentVelocity2;
} b3ManifoldConstraint;

typedef struct b3ContactConstraint
{
	b3ManifoldConstraint* constraints;
	struct b3Contact* contact;
	int indexA;
	int indexB;
	float invMassA, invMassB;
	b3Matrix3 invIA, invIB;
	b3Softness softness;
	b3Matrix3 rollingMass;
	float friction;
	float restitution;
	float rollingResistance;
	// Terra: the impulse limit over the step (FLT_MAX for none) and the impulse the finished sub-steps applied
	float maxNormalImpulse;
	float spentNormalImpulse;
	// Terra: transient links rebuilt after scalar prepare. Only the leader executes this body pair.
	struct b3ContactConstraint* impulseGroup;
	struct b3ContactConstraint* nextImpulseMember;
	struct b3ContactConstraint* impulseGroupTail;
	float sharedMaxNormalImpulse;
	float sharedSpentNormalImpulse;
	int impulseTargetBodyId;
	int impulseOpposingBodyId;
	int manifoldCount;
} b3ContactConstraint;

// Terra: step-local hash entry. Finalization copies exhaustion here before prepared constraints
// are released, so both ordinary and later bullet CCD can query the pair without dangling pointers.
typedef struct b3ContactImpulseGroupSlot
{
	b3ContactConstraint* leader;
	int targetBodyId;
	int opposingBodyId;
	bool exhausted;
} b3ContactImpulseGroupSlot;

int b3GetWideContactConstraintByteCount( void );

void b3PrepareContactImpulseGroups( b3StepContext* context );
void b3FinalizeContactImpulseGroups( b3StepContext* context );
bool b3ContactImpulseGroupAllowsPassThrough( const b3StepContext* context, int targetBodyId, int opposingBodyId );
void b3FreeContactImpulseGroups( b3StepContext* context );

// Overflow contacts don't fit into the constraint graph coloring
void b3PrepareContacts_Overflow( b3StepContext* context );
void b3WarmStartContacts_Overflow( b3StepContext* context );
void b3SolveContacts_Overflow( b3StepContext* context, bool useBias );
void b3ApplyRestitution_Overflow( b3StepContext* context );
void b3StoreImpulses_Overflow( b3StepContext* context );

void b3PrepareContacts_Mesh( b3SolverBlock block, b3StepContext* context );
void b3WarmStartContacts_Mesh( b3SolverBlock block, b3StepContext* context );
void b3SolveContacts_Mesh( b3SolverBlock block, b3StepContext* context, bool useBias );
void b3ApplyRestitution_Mesh( b3SolverBlock block, b3StepContext* context );
void b3StoreImpulses_Mesh( b3SolverBlock block, b3StepContext* context, int workerIndex );

void b3PrepareContacts_Convex( b3SolverBlock block, b3StepContext* context );
void b3WarmStartContacts_Convex( b3SolverBlock block, b3StepContext* context );
void b3SolveContacts_Convex( b3SolverBlock block, b3StepContext* context, bool useBias );
void b3ApplyRestitution_Convex( b3SolverBlock block, b3StepContext* context );
void b3StoreImpulses_Convex( b3SolverBlock block, b3StepContext* context, int workerIndex );
