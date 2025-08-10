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
    (c) 2018 Esi Ltd.

Class
    Foam::fixedValueZoneFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition constrains the solution of the equations inside
    cellZone called "inactive"+faceZoneName, where the faceZoneName is the
    faceZone name that the GIB is linked to.

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            fixedValueZone;
        value           uniform 0;
        cellZoneValue   0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    fixedValueZoneFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueZoneFvPatchField_H
#define fixedValueZoneFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedValueZoneFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValueZoneFvPatchField
:
    public fixedValueFvPatchField<Type>
{

private:

    //- Private members
    Type cellZoneValue_;


public:

    //- Runtime type information
    TypeName("fixedValueZone");


    // Constructors

        //- Construct from patch and internal field
        fixedValueZoneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueZoneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueZoneFvPatchField
        //  onto a new patch
        fixedValueZoneFvPatchField
        (
            const fixedValueZoneFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueZoneFvPatchField
        (
            const fixedValueZoneFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fixedValueZoneFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedValueZoneFvPatchField
        (
            const fixedValueZoneFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fixedValueZoneFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const Type& cellZoneValue() const
            {
                return cellZoneValue_;
            }

            virtual Type& cellZoneValue()
            {
                return cellZoneValue_;
            }

        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<Type>&
            );


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
