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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    triSurfModifier

Description
    A class allowing non-const access to the member of the surface mesh

SourceFiles
    triSurfModifier.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfModifier_H
#define triSurfModifier_H

#include "utilities/meshes/triSurf/triSurf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class triSurfModifier Declaration
\*---------------------------------------------------------------------------*/

class triSurfModifier
{
    // Private data
        //- reference to the surface mesh
        triSurf& surface_;

    // Disallow bitwise assignment
        //- Default construct
        triSurfModifier();

        //- Disallow default bitwise assignment
        void operator=(const triSurfModifier&);

public:

    // Constructors

        //- Construct from surface mesh
        triSurfModifier(triSurf& surface);

    // Destructor

        ~triSurfModifier();


    // Member Functions
        //- non-const access to points
        inline pointField& pointsAccess();

        //- access to facets
        inline LongList<labelledTri>& facetsAccess();

        //- non-const access to feature edges
        inline edgeLongList& featureEdgesAccess();

        //- access to patches
        inline geometricSurfacePatchList& patchesAccess();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/triSurf/triSurfModifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
