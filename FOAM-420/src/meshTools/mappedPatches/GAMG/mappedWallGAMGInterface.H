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
    (c) 2019 Esi Ltd.

Class
    Foam::mappedWallGAMGInterface

Description
    GAMG agglomerated coupled wall interface

SourceFiles
    mappedWallGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mappedWallGAMGInterface_H
#define mappedWallGAMGInterface_H

#include "mappedPatches/GAMG/mappedGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mappedGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class mappedWallGAMGInterface
:
    virtual public mappedGAMGInterface
{
public:

    //- Runtime type information
    TypeName("mappedWall");

    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        mappedWallGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Disallow default bitwise copy construct
        mappedWallGAMGInterface(const mappedWallGAMGInterface&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const mappedWallGAMGInterface&) = delete;


    //- Destructor
    virtual ~mappedWallGAMGInterface();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
