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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::fvMeshMapper

Description
    Class holds all the necessary information for mapping fields associated
    with fvMesh.
SourceFiles
    fvMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMapper_H
#define fvMeshMapper_H

#include "meshes/polyMesh/mapPolyMesh/faceMapper/faceMapper.H"
#include "meshes/polyMesh/mapPolyMesh/cellMapper/cellMapper.H"
#include "fvMesh/fvMeshMapper/fvSurfaceMapper.H"
#include "fvMesh/fvMeshMapper/fvBoundaryMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                        Class fvMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class fvMeshMapper
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Face mapper
        faceMapper faceMap_;

        //- Cell mapper
        cellMapper cellMap_;

        //- Surface mapper (needs to be shortened for internal faces only)
        fvSurfaceMapper surfaceMap_;

        //- Boundary mapper
        fvBoundaryMeshMapper boundaryMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvMeshMapper(const fvMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const fvMeshMapper&);


public:

    // Constructors

        //- Construct from fvMesh
        fvMeshMapper(const fvMesh& mesh, const mapPolyMesh& mpm)
        :
            mesh_(mesh),
            faceMap_(mpm),
            cellMap_(mpm),
            surfaceMap_(mesh, faceMap_),
            boundaryMap_(mesh, faceMap_)
        {}


    // Member Functions

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to objectRegistry storing fields. Can be
        //  removed once fields stored on pointMesh.
        const objectRegistry& thisDb() const
        {
            return mesh_;
        }

        //- Return volume mapper
        const morphFieldMapper& volMap() const
        {
            return cellMap_;
        }

        //- Return surface mapper
        const fvSurfaceMapper& surfaceMap() const
        {
            return surfaceMap_;
        }

        //- Return boundary mapper
        const fvBoundaryMeshMapper& boundaryMap() const
        {
            return boundaryMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
