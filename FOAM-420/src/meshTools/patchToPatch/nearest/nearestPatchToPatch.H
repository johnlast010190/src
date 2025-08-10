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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::patchToPatches::nearest

Description
    Class to generate patchToPatch coupling geometry. Couples a face to the
    single nearest opposite face only.

SourceFiles
    nearestPatchToPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nearestPatchToPatch_H
#define nearestPatchToPatch_H

#include "patchToPatch/patchToPatch/patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{

/*---------------------------------------------------------------------------*\
                           Class nearest Declaration
\*---------------------------------------------------------------------------*/

class nearest
:
    public patchToPatch
{
    // Private Member Data

        //- For each source face, the distance to its coupled target face
        List<scalar> srcDistances_;

        //- For each target face, the distance to its coupled source face
        List<scalar> tgtDistances_;


    // Private Member Functions

        //- Get the bound box for a source face
        virtual treeBoundBox srcBox
        (
            const face& srcFace,
            const pointField& srcPoints,
            const vectorField& srcPointNormals
        ) const override;

        //- Intersect two faces
        bool intersectFaces
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const label facei,
            const label otherFacei,
            DynamicList<label>& faceOtherFaces,
            scalar& faceDistance
        ) const;

        //- Intersect two faces
        virtual bool intersectFaces
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const label srcFacei,
            const label tgtFacei
        ) override;

        //- Initialisation
        virtual void initialise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch
        ) override;

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes.
        virtual void rDistributeTgt
        (
            const primitiveOldTimePatch& tgtPatch,
            const mapDistribute& tgtMap
        ) override;


public:

    //- Runtime type information
    TypeName("nearest");


    // Constructors

        //- Construct from components
        nearest(const bool reverse);


    //- Destructor
    ~nearest();


    // Member Functions

        //- For each source face, the coupled target weights
        virtual tmpNrc<List<DynamicList<scalar>>> srcWeights
        (
            const primitivePatch& srcPatch
        ) const override;

        //- For each target face, the coupled source weights
        virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights
        (
            const primitivePatch& tgtPatch
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
