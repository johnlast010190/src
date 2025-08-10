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
    fpmaMesh

Description

SourceFiles
    fpmaMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fpmaMesh_H
#define fpmaMesh_H

#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fpmaMesh Declaration
\*---------------------------------------------------------------------------*/

class fpmaMesh
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fpmaMesh(const fpmaMesh&);

        //- Disallow default bitwise assignment
        void operator=(const fpmaMesh&);

        void writePoints(OFstream& fpmaGeometryFile) const;

        void writeFaces(OFstream& fpmaGeometryFile) const;

        void writeCells(OFstream& fpmaGeometryFile) const;

        void writeSubsets(OFstream& fpmaGeometryFile) const;

public:

    // Public data

        const polyMeshGen& mesh_;

    // Constructors

        //- Construct from polyMeshGen
        fpmaMesh(const polyMeshGen& mesh);


    // Destructor

        ~fpmaMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

            void write(OFstream& fpmaGeometryFile) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
