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
    (c) 2016 Esi Ltd.

Class
    coupledFaePatchField

Description

SourceFiles
    coupledFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFaePatchField_H
#define coupledFaePatchField_H

#include "fields/faePatchFields/faePatchField/faePatchField.H"
#include "faMesh/faPatches/basic/coupled/coupledFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coupledFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledFaePatchField
:
    public faePatchField<Type>
{

public:

    //- Runtime type information
    TypeName(coupledFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        coupledFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct from patch and internal field and patch field
        coupledFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        coupledFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given coupledFaePatchField onto a new patch
        coupledFaePatchField
        (
            const coupledFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        coupledFaePatchField
        (
            const coupledFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone() const = 0;

        //- Construct as copy setting internal field reference
        coupledFaePatchField
        (
            const coupledFaePatchField<Type>&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, faEdgeMesh>&
        ) const = 0;


    // Member functions

        // Access

            //- Return true if this patch field is derived from
            //  coupledFaePatchField<Type>.
            virtual bool coupled() const
            {
                return true;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fields/faePatchFields/basic/coupled/coupledFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
