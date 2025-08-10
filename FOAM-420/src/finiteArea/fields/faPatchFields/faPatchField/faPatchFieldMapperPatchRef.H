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
    (c) 1991-2005 OpenCFD Ltd.

Class
    faPatchFieldMapperPatchRef

Description

\*---------------------------------------------------------------------------*/

#ifndef faPatchFieldMapperPatchRef_H
#define faPatchFieldMapperPatchRef_H

#include "fields/faPatchFields/faPatchField/faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class faPatchFieldMapperPatchRef Declaration
\*---------------------------------------------------------------------------*/

class faPatchFieldMapperPatchRef
:
    public faPatchFieldMapper
{
    // Private data

        const faPatch& sourcePatch_;
        const faPatch& targetPatch_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faPatchFieldMapperPatchRef(const faPatchFieldMapperPatchRef&);

        //- Disallow default bitwise assignment
        void operator=(const faPatchFieldMapperPatchRef&);


public:

    // Constructors

        //- Construct from components
        faPatchFieldMapperPatchRef(const faPatch& source, const faPatch& target)
        :
            sourcePatch_(source),
            targetPatch_(target)
        {}


    // Member functions

        const faPatch& sourcePatch() const
        {
            return sourcePatch_;
        }

        const faPatch& targetPatch() const
        {
            return targetPatch_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
