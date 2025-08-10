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
    (c) 2022 Esi Ltd.

Class
    Foam::distributedFvPatchFieldMapper

Description
    A fvPatchFieldMapper with distributed functionality.

\*---------------------------------------------------------------------------*/

#ifndef distributedFvPatchFieldMapper_H
#define distributedFvPatchFieldMapper_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributedFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedFvPatchFieldMapper
:
    public generalFvPatchFieldMapper
{
protected:

    // Protected Data

        const mapDistributeBase* distMapPtr_;

        bool hasUnmapped_;


public:

    // Constructors

        //- Construct given addressing
        distributedFvPatchFieldMapper(const mapDistributeBase* distMapPtr)
        :
            distMapPtr_(distMapPtr),
            hasUnmapped_(false)
        {}


    //- Destructor
    virtual ~distributedFvPatchFieldMapper()
    {}


    // Member Functions

        //- Is the mapping distributed
        virtual bool distributed() const = 0;


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, , = 0)

        //- Map a label field
        DECLARE_DISTRIBUTED_MAPPER_OPERATOR(label, , = 0)

        //- Specialisations for VectorN types
        forAllVectorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, = 0)

        forAllTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, = 0)

        forAllDiagTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, = 0)

        forAllSphericalTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, = 0)
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
