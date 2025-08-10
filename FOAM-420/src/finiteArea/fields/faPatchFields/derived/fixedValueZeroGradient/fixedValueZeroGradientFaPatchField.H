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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016-2017 Wikki Ltd
    (c) 2018 Esi Ltd.

Class
    Foam::fixedValueZeroGradientFaPatchField

Description
    This boundary condition provides a fixed value condition and a zero gradient
    condition.

See also
    Foam::fixedValueFaPatchField

SourceFiles
    fixedValueZeroGradientFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueZeroGradientFaPatchField_H
#define fixedValueZeroGradientFaPatchField_H

#include "fields/faPatchFields/basic/fixedValue/fixedValueFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedValueZeroGradientFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValueZeroGradientFaPatchField
:
    public fixedValueFaPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("fixedValueZeroGradient");


    // Constructors

        //- Construct from patch and internal field
        fixedValueZeroGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueZeroGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        fixedValueZeroGradientFaPatchField
        (
            const fixedValueZeroGradientFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueZeroGradientFaPatchField
        (
            const fixedValueZeroGradientFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return tmp<faPatchField<Type>>
            (
                new fixedValueZeroGradientFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedValueZeroGradientFaPatchField
        (
            const fixedValueZeroGradientFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type>>
            (
                new fixedValueZeroGradientFaPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~fixedValueZeroGradientFaPatchField()
    {}



    // Member functions

        //- Manipulate matrix
        virtual void manipulateMatrix(faMatrix<Type>& matrix);

        //- Manipulate matrix with given weights
        virtual void manipulateMatrix
        (
            faMatrix<Type>& matrix,
            const Field<Type>& weights
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/faPatchFields/derived/fixedValueZeroGradient/fixedValueZeroGradientFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
