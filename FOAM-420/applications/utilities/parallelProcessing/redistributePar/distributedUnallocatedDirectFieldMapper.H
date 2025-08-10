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
    Foam::distributedUnallocatedDirectFieldMapper

Description
    FieldMapper with direct mapping from remote quantities.

\*---------------------------------------------------------------------------*/

#ifndef distributedUnallocatedDirectFieldMapper_H
#define distributedUnallocatedDirectFieldMapper_H

#include "fields/Fields/fieldMappers/directFieldMapper/directFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class distributedUnallocatedDirectFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedUnallocatedDirectFieldMapper
:
    public directFieldMapper
{
    // Private Data

        const labelUList& addressing_;

        const mapDistributeBase& distMap_;

        bool hasUnmapped_;


    // Private Member Functions

        label size() const
        {
            return
            (
                notNull(addressing_)
              ? addressing_.size()
              : distMap_.constructSize()
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
        distributedUnallocatedDirectFieldMapper
        (
            const labelUList& addressing,
            const mapDistributeBase& distMap
        );


    //- Destructor
    virtual ~distributedUnallocatedDirectFieldMapper()
    {}


    // Member Functions

        //- Is the mapping distributed
        bool distributed() const
        {
            return true;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        //- Return the direct map addressing
        virtual const labelUList& addressing() const
        {
            return addressing_;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DECLARE_DISTRIBUTED_MAPPER_OPERATOR, ,)
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distributedUnallocatedDirectFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
