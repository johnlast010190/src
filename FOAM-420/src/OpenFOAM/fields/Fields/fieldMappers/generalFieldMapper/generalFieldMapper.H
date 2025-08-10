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
    (c) 2019-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::generalFieldMapper

Description
    General field mapper supporting both direct, indirect and weighted mapping.

\*---------------------------------------------------------------------------*/

#ifndef generalFieldMapper_H
#define generalFieldMapper_H

#include "fields/Fields/fieldMappers/fieldMapper/fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class generalFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class generalFieldMapper
:
    virtual public fieldMapper
{
    // Private Member Functions

        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& f) const;


public:

    // Constructors

        //- Null constructor
        generalFieldMapper()
        {}


    //- Destructor
    virtual ~generalFieldMapper()
    {}


    // Member Functions

        //- Is the mapping direct
        virtual bool direct() const = 0;

        //- Is the mapping indirect
        virtual bool indirect() const override = 0;

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override = 0;

        //- Return the direct map addressing
        virtual const labelUList& directAddressing() const;

        //- Return the indirect map addressing
        virtual const labelUList& indirectAddressing() const;

        //- Return the interpolated map addressing
        virtual const labelListList& addressing() const;

        //- Return interpolation weights
        virtual const scalarListList& weights() const;


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DECLARE_FIELD_MAPPER_OPERATOR, , override)

        //- Map a label field
        DECLARE_FIELD_MAPPER_OPERATOR(label, , override)

        //- Specialisations for VectorN types
        forAllVectorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, override)

        forAllTensorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, override)

        forAllDiagTensorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, override)

        forAllSphericalTensorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, override)
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/Fields/fieldMappers/generalFieldMapper/generalFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
