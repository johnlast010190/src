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

Class
    Foam::IOmapDistributePolyMesh

Description
    IOmapDistributePolyMesh is derived from mapDistributePolyMesh and
    IOobject to give the mapDistributePolyMesh
    automatic IO functionality via the objectRegistry.

SourceFiles
    IOmapDistributePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef IOmapDistributePolyMesh_H
#define IOmapDistributePolyMesh_H

#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IOmapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class IOmapDistributePolyMesh
:
    public regIOobject,
    public mapDistributePolyMesh
{

public:

    //- Runtime type information
    TypeName("mapDistributePolyMesh");

    // Constructors

        //- Construct given an IOobject
        IOmapDistributePolyMesh(const IOobject&);

        //- Construct given an IOobject and mapDistributePolyMesh
        IOmapDistributePolyMesh(const IOobject&, const mapDistributePolyMesh&);


    //- Destructor
    virtual ~IOmapDistributePolyMesh();


    // Member functions

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
