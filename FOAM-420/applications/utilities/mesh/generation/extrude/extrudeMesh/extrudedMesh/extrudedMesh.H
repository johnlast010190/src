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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::extrudedMesh

Description

SourceFiles
    extrudedMesh.C
    extrudedMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extrudedMesh_H
#define extrudedMesh_H

#include "meshes/polyMesh/polyMesh.H"
#include "extrudeModel/extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class extrudedMesh Declaration
\*---------------------------------------------------------------------------*/

class extrudedMesh
:
    public polyMesh
{
    // Private data

        const extrudeModel& model_;


    // Private Member Functions

        //- Do edge and face use points in same order?
        static bool sameOrder(const face&, const edge&);

        //- Construct and return the extruded mesh points
        template<class FaceList, class PointField>
        Xfer<pointField> extrudedPoints
        (
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );

        //- Construct and return the extruded mesh faces
        template<class FaceList, class PointField>
        Xfer<faceList> extrudedFaces
        (
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );

        //- Construct and return the extruded mesh cells
        template<class FaceList, class PointField>
        Xfer<cellList> extrudedCells
        (
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );


        //- Disallow default bitwise copy construct
        extrudedMesh(const extrudedMesh&);

        //- Disallow default bitwise assignment
        void operator=(const extrudedMesh&);


public:

    // Constructors

        //- Construct from the primitivePatch to extrude
        template<class FaceList, class PointField>
        extrudedMesh
        (
            const IOobject&,
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extrudedMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
