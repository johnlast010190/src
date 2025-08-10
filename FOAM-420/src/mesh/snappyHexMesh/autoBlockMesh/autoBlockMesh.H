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
    (c) 2019 Esi Ltd.

Application
    Foam::autoBlockMesh

Description
    Automatically generate a starting blockMesh for foamHexMesh

\*---------------------------------------------------------------------------*/

#ifndef autoBlockMesh_H
#define autoBlockMesh_H

#include "refinementSurfaces/refinementSurfaces.H"
#include "coordinate/systems/coordinateSystem.H"
#include "primitives/Vector/labelVector/labelVector.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class autoBlockMesh Declaration
\*---------------------------------------------------------------------------*/

class autoBlockMesh
{

    // Private Member Functions

    fvMesh& mesh_;

    const refinementSurfaces& surfaces_;

    coordinateFrame* coordFramePtr_;
    autoPtr<coordinateSystem> coordSysPtr_;

    const Tuple2<scalar, label> blockData_;

    const Switch primitiveBounding_;

    const Switch blockIso_;

    const point blockOrigin_;

public:

    //- Runtime type information
    ClassName("autoBlockMesh");

    // Constructors

        //- Construct from refinement surfaces
        autoBlockMesh
        (
            const refinementSurfaces& surfaces,
            const dictionary& meshDict,
            fvMesh& mesh
        );

    // Member Functions

        // Helper

            //Add an empty patch
            void addEmptyPatch
            (
                const word& patchName,
                const dictionary& patchInfo
            );

            //Read coordinate rotation dictonary or set to default
            dictionary getCoordinateDict
            (
                const dictionary& meshDict
            );

            //- point label offset for a particular i,j,k position
            inline label pointLabel
            (
                const label i,
                const label j,
                const label k,
                const labelVector n
            ) const;

            //- cell label offset for a particular i,j,k position
            inline label cellLabel
            (
                const label i,
                const label j,
                const label k,
                const labelVector n
            ) const;

        // Access

        // create the block mesh

            void createBlock();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
