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
    (c) 2010-2016 Esi Ltd.

Class
    Foam::mappedLduInterface

Description
    An abstract base class for region coupled interfaces

SourceFiles
    mappedLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mappedLduInterface_H
#define mappedLduInterface_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "AMIInterpolation/AMIInterpolation/AMIPatchToPatchInterpolation.H"
#include "meshes/polyMesh/polyMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mappedLduInterface Declaration
\*---------------------------------------------------------------------------*/

class mappedLduInterface
{

public:

    //- Runtime type information
    TypeName("mappedLduInterface");


    // Constructors

        //- Construct null
        mappedLduInterface()
        {}


    //- Destructor
    virtual ~mappedLduInterface();


    // Member Functions

        // Access

            //- Whether this can indeed act as an implicitly coupled interface
            virtual bool isInterface() const = 0;

            //- Return AMI
            virtual const AMIPatchToPatchInterpolation& AMI() const = 0;

            //- Return direct map
            virtual const mapDistribute& map() const = 0;

            //- Whether we are using AMI or direct map
            virtual bool usingAMI() const = 0;

            virtual bool owner() const = 0;

            virtual const word& patchName() const = 0;

            virtual const word& regionName() const = 0;

            virtual const word& nbrPatchName() const = 0;

            virtual const word& nbrRegionName() const = 0;

            //- Return neighbour patch
            virtual const mappedLduInterface& nbrPatch() const = 0;

            //- Whether neighbour mesh exists
            virtual bool nbrMeshFound() const = 0;

            //- Return nbrMesh
            virtual const polyMesh& nbrMesh() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
