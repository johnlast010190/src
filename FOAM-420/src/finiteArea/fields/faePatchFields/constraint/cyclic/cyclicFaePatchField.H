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
    Foam::cyclicFaePatchField

Description
    Foam::cyclicFaePatchField

SourceFiles
    cyclicFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFaePatchField_H
#define cyclicFaePatchField_H

#include "fields/faePatchFields/basic/coupled/coupledFaePatchField.H"
#include "faMesh/faPatches/constraint/cyclic/cyclicFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFaePatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicFaePatchField
:
    public coupledFaePatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicFaPatch& cyclicPatch_;


public:

    //- Runtime type information
    TypeName(cyclicFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicFaePatchField onto a new patch
        cyclicFaePatchField
        (
            const cyclicFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        cyclicFaePatchField
        (
            const cyclicFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return tmp<faePatchField<Type>>
            (
                new cyclicFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicFaePatchField
        (
            const cyclicFaePatchField<Type>&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, faEdgeMesh>& iF
        ) const
        {
            return tmp<faePatchField<Type>>
            (
                new cyclicFaePatchField<Type>(*this, iF)
            );
        }

// Member functions

    // Access

        virtual word constraintType()
        {
            return cyclicFaPatch::typeName_();
        }


};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fields/faePatchFields/constraint/cyclic/cyclicFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
