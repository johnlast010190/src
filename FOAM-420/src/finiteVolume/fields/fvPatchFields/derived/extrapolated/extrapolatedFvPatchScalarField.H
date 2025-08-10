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
    (c) 2010-2012 Esi Ltd.

Class
    extrapolatedFvPatchScalarField

Description

    Uses cell centred gradient to extrapolate boundary value

SourceFiles
    extrapolatedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedFvPatchScalarField_H
#define extrapolatedFvPatchScalarField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class extrapolatedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extrapolatedFvPatchScalarField
:
    public fvPatchField<scalar>
{
private:

        //- patch intrnal gradient
        vectorField patchInternalGradient_;

        //- vectors from cell centres to face centres
        tmp<vectorField> patchCellToFace_;

    // Private member functions

        //- calculate patch internal cell centres
        tmp<vectorField> calcPatchCellToFace();

        //- update patch internal gradient
        void calcPatchInternalGradient();

protected:

        static labelList fieldGradCalcTime_;
        static wordList fieldGradCalcName_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("extrapolated");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        extrapolatedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Field<scalar>&
        );

        //- Construct by mapping given extrapolatedFvPatchScalarField
        //  onto a new patch
        extrapolatedFvPatchScalarField
        (
            const extrapolatedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct as copy
        extrapolatedFvPatchScalarField(const extrapolatedFvPatchScalarField&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new extrapolatedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extrapolatedFvPatchScalarField
        (
            const extrapolatedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (const DimensionedField<scalar, volMesh>& iF) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new extrapolatedFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        // Access


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<scalar>> snGrad() const;

            // update gradients
            virtual void updateCoeffs();

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes
                commsType=Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
