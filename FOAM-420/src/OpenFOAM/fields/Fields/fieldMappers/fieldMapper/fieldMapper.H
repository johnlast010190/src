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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::fieldMapper

Description
    Abstract base class for field mapping

\*---------------------------------------------------------------------------*/

#ifndef fieldMapper_H
#define fieldMapper_H

#include "fields/Fields/Field/Field.H"
#include "fields/Fields/fieldMappers/fieldMapper/fieldMapperM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldMapper Declaration
\*---------------------------------------------------------------------------*/

class fieldMapper
{
public:

    // Constructors

        //- Null constructor
        fieldMapper()
        {}


    //- Destructor
    virtual ~fieldMapper()
    {}


    // Member Functions

        //- Is the mapping indirect
        virtual bool indirect() const = 0;

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const = 0;


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DECLARE_FIELD_MAPPER_OPERATOR, , = 0)

        //- Map a label field
        DECLARE_FIELD_MAPPER_OPERATOR(label, , = 0)

        //- Specialisations for VectorN types
        forAllVectorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, = 0)

        forAllTensorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, = 0)

        forAllDiagTensorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, = 0)

        forAllSphericalTensorNTypes(DECLARE_FIELD_MAPPER_OPERATOR, = 0)

        //- Map a temporary field
        template<class Type>
        void operator()(Field<Type>& f, const tmp<Field<Type>>& tmapF) const;

        //- Map a temporary field
        template<class Type>
        tmp<Field<Type>> operator()(const tmp<Field<Type>>& tmapF) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/Fields/fieldMappers/fieldMapper/fieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
