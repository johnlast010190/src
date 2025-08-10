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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::displacementMotionSolverMeshMover

Description
    Quality-based under-relaxation wrapped around generic
    displacementMotionSolver.

    Example of use in layer settings in snappyHexMeshDict:
    \verbatim
    meshShrinker displacementMotionSolver;
    solver displacementLaplacian;
    displacementLaplacianCoeffs
    {
        diffusivity quadratic inverseDistance 1(wall);
    }
    \endverbatim

SourceFiles
    displacementMotionSolverMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMotionSolverMeshMover_H
#define displacementMotionSolverMeshMover_H

#include "externalDisplacementMeshMover/externalDisplacementMeshMover.H"
#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "motionSmoother/motionSmootherAlgo.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class displacementMotionSolverMeshMover Declaration
\*---------------------------------------------------------------------------*/

class displacementMotionSolverMeshMover
:
    public externalDisplacementMeshMover
{
    // Private Data

        //- Mesh motion solver
        autoPtr<displacementMotionSolver> solverPtr_;

        //- IDs of fixedValue patches that we can modify
        const labelList adaptPatchIDs_;

        //- Combined indirect fixedValue patches that we can modify
        autoPtr<indirectPrimitivePatch> adaptPatchPtr_;

        //- Scale factor for displacemen
        pointScalarField scale_;

        //- Old point field
        pointField oldPoints_;

        //- Mesh mover algorithm
        motionSmootherAlgo meshMover_;

        //- Field smoothing
        fieldSmoother fieldSmoother_;


    // Private Member Functions

        //- Apply the mesh mover algorithm
        bool moveMesh
        (
            const dictionary& moveDict,
            const label nAllowableErrors,
            labelList& checkFaces
        );


public:

    //- Runtime type information
    TypeName("displacementMotionSolver");


    // Constructors

        //- Construct from a polyMesh and an IOdictionary
        displacementMotionSolverMeshMover
        (
            const dictionary& dict,
            const List<labelPair>& baffles,
            pointVectorField& pointDisplacemen
        );


    //- Destructor
    virtual ~displacementMotionSolverMeshMover();


    // Member Functions

        //- Move mesh using current pointDisplacement boundary values.
        //  Return true if succesful (errors on checkFaces less than
        //  allowable). Updates pointDisplacement.
        virtual bool move
        (
            const dictionary&,
            const label nAllowableErrors,
            const PackedList<1>& isConvexEdgePoint,
            const PackedList<1>& isConcaveEdgePoint,
            labelList& checkFaces
        );

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&)
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
