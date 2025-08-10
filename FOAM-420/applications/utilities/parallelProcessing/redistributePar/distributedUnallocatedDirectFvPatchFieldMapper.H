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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::distributedUnallocatedDirectFvPatchFieldMapper

Description
    FieldMapper with direct mapping from remote quantities.

\*---------------------------------------------------------------------------*/

#ifndef distributedUnallocatedDirectFvPatchFieldMapper_H
#define distributedUnallocatedDirectFvPatchFieldMapper_H

#include "fields/fvPatchFields/fvPatchField/fieldMappers/distributedFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class distributedUnallocatedDirectFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedUnallocatedDirectFvPatchFieldMapper
:
    public distributedFvPatchFieldMapper
{
    // Private Data

        const labelUList& directAddressing_;


    // Private Member Functions

        label size() const
        {
            return
            (
                notNull(directAddressing_)
              ? directAddressing_.size()
              : distMapPtr_->constructSize()
            );
        }

        template<class Type>
        void map
        (
            Field<Type>& f,
            const Field<Type>& mapF,
            const bool applyFlip
        ) const;

        template<class Type>
        tmp<Field<Type>> map
        (
            const Field<Type>& f,
            const bool applyFlip
        ) const;


public:

    // Constructors

        //- Construct given addressing
        distributedUnallocatedDirectFvPatchFieldMapper
        (
            const labelUList& directAddressing,
            const mapDistributeBase& distMap
        );


    //- Destructor
    virtual ~distributedUnallocatedDirectFvPatchFieldMapper()
    {}


    // Member Functions

        //- Is the mapping direct
        virtual bool direct() const override
        {
            return true;
        }

        //- No indirect mapping
        virtual bool indirect() const override
        {
            return false;
        }

        //- Is the mapping distributed
        virtual bool distributed() const override
        {
            return true;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return hasUnmapped_;
        }

        //- Return the direct map addressing
        virtual const labelUList& directAddressing() const override
        {
            return directAddressing_;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, , override)

        //- Map a label field
        DECLARE_DISTRIBUTED_MAPPER_OPERATOR(label, , override)

        //- Specialisations for VectorN types
        forAllVectorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, override)

        forAllTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, override)

        forAllDiagTensorNTypes(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, override)

        forAllSphericalTensorNTypes
        (
            DECLARE_DISTRIBUTED_MAPPER_OPERATOR,
            override
        )
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distributedUnallocatedDirectFvPatchFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
