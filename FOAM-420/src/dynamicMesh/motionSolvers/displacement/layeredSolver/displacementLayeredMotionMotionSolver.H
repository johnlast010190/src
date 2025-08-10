/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Class
    Foam::displacementLayeredMotionMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an (multi-block) extruded fvMesh. Gets given the
    structure of the mesh blocks and boundary conditions on these blocks.

    The displacementLayeredMotionCoeffs subdict of dynamicMeshDict specifies
    per region (=cellZone) the boundary conditions on two opposing patches
    (=faceZones). It then interpolates the boundary values using topological
    walking so requires the cellZone to be a layered mesh.

    The boundary conditions on faceZones are currently:

    follow: the faceZone follows the overall mesh displacement.
            Use this for faceZones on boundary faces (so it uses the
            proper boundary conditions on the pointDisplacement).

    uniformFollow: like 'follow' but takes the average value of
            a specified 'patch' (which is not necessarily colocated)

    fixedValue: fixed value.

    timeVaryingUniformFixedValue: table-driven fixed value.

    slip: the second of a pair of faceZones follows the tangential movement
          specified by the first faceZone. (= removes the normal component).

SourceFiles
    displacementLayeredMotionMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementLayeredMotionMotionSolver_H
#define displacementLayeredMotionMotionSolver_H

#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "containers/Lists/PackedList/PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations

/*---------------------------------------------------------------------------*\
             Class displacementLayeredMotionMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementLayeredMotionMotionSolver
:
    public displacementMotionSolver
{
    // Private Member Functions

        void calcZoneMask
        (
            const label cellZoneI,
            PackedBoolList& isZonePoint,
            PackedBoolList& isZoneEdge
        ) const;

        void walkStructured
        (
            const label cellZoneI,
            const PackedBoolList& isZonePoint,
            const PackedBoolList& isZoneEdge,
            const labelList& seedPoints,
            const vectorField& seedData,
            scalarField& distance,
            vectorField& data
        ) const;

        tmp<vectorField> faceZoneEvaluate
        (
            const faceZone& fz,
            const labelList& meshPoints,
            const dictionary& dict,
            const PtrList<pointVectorField>& patchDisp,
            const label patchi
        ) const;

        void cellZoneSolve
        (
            const label cellZoneI,
            const dictionary& zoneDict
        );


        //- Disallow default bitwise copy construct
        displacementLayeredMotionMotionSolver
        (
            const displacementLayeredMotionMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementLayeredMotionMotionSolver&);


public:

    //- Runtime type information
    TypeName("displacementLayeredMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        displacementLayeredMotionMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Construct from polyMesh, dictionary and displacement and
        //  zero-time fields
        displacementLayeredMotionMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const pointVectorField&,
            const pointIOField&
        );


    //- Destructor
    ~displacementLayeredMotionMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
