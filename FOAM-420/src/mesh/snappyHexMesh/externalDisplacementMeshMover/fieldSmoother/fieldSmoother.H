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
    (c) 2015 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::fieldSmoother

Description
    Utility functions for mesh motion solvers

SourceFiles
    fieldSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef fieldSmoother_H
#define fieldSmoother_H

#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldSmoother Declaration
\*---------------------------------------------------------------------------*/

class fieldSmoother
{
    // Private data

        //- Reference to the poly mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fieldSmoother(const fieldSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const fieldSmoother&);


public:

    // Run-time type information
    TypeName("fieldSmoother");


    // Constructors

        //- Construct from a polyMesh
        fieldSmoother(const polyMesh&);


    //- Destructor
    virtual ~fieldSmoother();


    // Member Functions

        //- Smooth interior normals
        void smoothNormals
        (
            const label nIter,
            const PackedBoolList& isMeshMasterPoint,
            const PackedBoolList& isMeshMasterEdge,
            const labelList& fixedPoints,
            pointVectorField& normals
        ) const;

        //- Smooth patch normals
        void smoothPatchNormals
        (
            const label nIter,
            const PackedBoolList& isPatchMasterPoint,
            const PackedBoolList& isPatchMasterEdge,
            const boolList& isExternalEdge,
            const boolList& isExternalPoint,
            const boolList& isCornerPoint,
            const labelList& meshEdges,
            const indirectPrimitivePatch& adaptPatch,
            pointField& normals
        ) const;

        //- Smooth a scalar field towards, but not beyond, a minimum value
        template <class Type>
        void minSmoothField
        (
            const label nIter,
            const PackedBoolList& isPatchMasterPoint,
            const PackedBoolList& isPatchMasterEdge,
             const labelList& meshEdges,
            const indirectPrimitivePatch& adaptPatch,
            const scalarField& fieldMin,
            Field<Type>& field
        ) const;

        //- Smooth and then un-smooth a displacement
        void smoothLambdaMuDisplacement
        (
            const label nIter,
            const PackedBoolList& isMeshMasterPoint,
            const PackedBoolList& isMeshMasterEdge,
            const PackedBoolList& isSmoothable,
            vectorField& displacement
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmootherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
