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
    (c) 2017 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.
    (c) 2020 Esi Ltd.

Class
    Foam::methToSurface

Description
    Create and export surface from mesh patches and zones

SourceFiles
    meshToSurface.C

\*---------------------------------------------------------------------------*/

#ifndef meshToSurface_H
#define meshToSurface_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class meshToSurface Declaration
\*---------------------------------------------------------------------------*/

class meshToSurface
{
    // Private data

        //- polyMesh reference
        const polyMesh& mesh_;

        //- Hash set with the patches
        const labelHashSet& patchHash_;

        //- Hash set with the zones
        const labelHashSet& zoneHash_;

        //- Keep surface order consistent with the mesh patches
        const bool consistent_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshToSurface(const meshToSurface&);

        //- Disallow default bitwise assignment
        void operator=(const meshToSurface&);


public:

    //- Runtime type information
    ClassName("meshToSurface");


    // Constructors

        //- Construct from cells to cut and pattern of cuts
        meshToSurface
        (
            const polyMesh&,
            const labelHashSet&,
            const labelHashSet&,
            const bool&
        );


    //- Destructor
    ~meshToSurface(){};


    // Member Functions

        // Access


        // Write

            //- Create and write mesh surface
            void write(const fileName&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
