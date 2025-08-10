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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::orientFaceZone

Description
    Corrects orientation of faceZone.

    - correct in parallel - excludes coupled faceZones from walk
    - correct for non-manifold faceZones - restarts walk

SourceFiles
    orientFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef orientFaceZone_H
#define orientFaceZone_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class orientFaceZone Declaration
\*---------------------------------------------------------------------------*/

class orientFaceZone
{
    // Private data

        //- polyMesh reference
        polyMesh& mesh_;

        //- faceZone name
        const word& fzName_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        orientFaceZone(const orientFaceZone&);

        //- Disallow default bitwise assignment
        void operator=(const orientFaceZone&);


public:

    //- Runtime type information
    ClassName("orientFaceZone");


    // Constructors

        //- Construct from cells to cut and pattern of cuts
        orientFaceZone
        (
            polyMesh&,
            const word&
        );


    //- Destructor
    ~orientFaceZone(){};


    // Member Functions

        // Access


        // Write

            //- Create and write mesh surface
            void orientByPoint
            (
                const point&
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
