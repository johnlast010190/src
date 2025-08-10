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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2015 Esi Ltd.

Class
    Foam::lduRegisteredPrimitiveMesh

Description
    lduPrimitiveMesh which stores an object registry (needed by GAMG).

SourceFiles
    lduRegisteredPrimitiveMesh.C

\*---------------------------------------------------------------------------*/

#ifndef lduRegisteredPrimitiveMesh_H
#define lduRegisteredPrimitiveMesh_H

#include "meshes/lduMesh/lduPrimitiveMesh.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduRegisteredPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduRegisteredPrimitiveMesh
:
    public lduPrimitiveMesh
{
    // Private data

        //- Object registry returned by thisDb
        const objectRegistry& obr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduRegisteredPrimitiveMesh(const lduRegisteredPrimitiveMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduRegisteredPrimitiveMesh&);


public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("lduRegisteredPrimitiveMesh");

    // Constructors

        //- Construct from components and re-use storage.
        lduRegisteredPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            PtrList<const lduInterface>& primitiveInterfaces,
            const lduSchedule& ps,
            const label comm,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~lduRegisteredPrimitiveMesh()
    {}


    // Member Functions

        // Access

            //- Return the object registry
            virtual const objectRegistry& thisDb() const
            {
                return obr_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
