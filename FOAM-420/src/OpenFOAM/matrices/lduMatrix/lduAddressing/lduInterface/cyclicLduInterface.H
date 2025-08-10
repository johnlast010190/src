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
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicLduInterface

Description
    An abstract base class for cyclic coupled interfaces.

SourceFiles
    cyclicLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicLduInterface_H
#define cyclicLduInterface_H

#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "primitives/transform/transformer/transformer.H"
#include "fields/Fields/primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicLduInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicLduInterface
{

public:

    //- Runtime type information
    TypeName("cyclicLduInterface");

    // Constructors

    //- Destructor
    virtual ~cyclicLduInterface();


    // Member Functions

        // Access

            //- Return neighbour
            virtual label nbrPatchID() const = 0;

            virtual bool owner() const = 0;

            //- Return processor number
            virtual const cyclicLduInterface& nbrPatch() const = 0;

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
