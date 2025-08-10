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
    symmetryPlaneOptimisation

Description
    Smoothing of symmetry planes in the mesh such that all points
    are in the plane.

SourceFiles
    symmetryPlaneOptimisation.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneOptimisation_H
#define symmetryPlaneOptimisation_H

#include "utilities/containers/DynList/DynList.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "meshes/boundBox/boundBox.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/primitiveShapes/plane/plane.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class symmetryPlaneOptimisation Declaration
\*---------------------------------------------------------------------------*/

class symmetryPlaneOptimisation
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- symmetry planes in the mesh
        std::map<label, plane> symmetryPlanes_;

    // Private member functions
        //- detect symmetry planes
        void detectSymmetryPlanes();

        //- point-planes addressing
        void pointInPlanes(VRWGraph&) const;

public:

    // Constructors

        //- Construct from mesh
        symmetryPlaneOptimisation(polyMeshGen& mesh);


    // Destructor

        ~symmetryPlaneOptimisation();

    // Member Functions
        //- move vertices to the symmetry planes
        void optimizeSymmetryPlanes();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
