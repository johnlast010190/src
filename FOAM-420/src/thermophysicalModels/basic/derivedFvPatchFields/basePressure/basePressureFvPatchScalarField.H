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
    Foam::basePressureFvPatchScalarField

Description
    Dirichlet BC for pressure which depends on U

    p = p0 + (C1 + C2*v_i n_i) v_i + C3*(v_i t_i)^2
    or
    p = p0 + (C1 + rho*C2*v_i n_i) v_i + rho*C3*(v_i t_i)^2

    This BC is used specifically for block solver because it provides off
    diagonal matrix contributions resulted from in grad(p), div(U), laplacian(p)
    The BC can be used in many different ways. Most common usages are
    - for totalPressure with C0 = 0 and C2=-0.5 C3=0.5
    - For resistive outlets where the C0 C1 C2 depends on the resistive model
    (Darcy-Forchheimer -- D is C1 and C2, C3 if F).

See also
    Foam::fixedValueFvPatchField
    Foam::pressureVelocityFvPatchVectorField
    Foam::porousPressureFvPatchScalarField
    Foam::totalPressurePressureFvPatchScalarField

SourceFiles
    basePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef basePressureFvPatchScalarField_H
#define basePressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class basePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class basePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:

    // Private data

        //- Total pressure
        scalarField p0_;

        //- Name of the velocity field
        word UName_;

        //- Name of the density field. Needed for compressible formulation
        word rhoName_;

        //- Ability to add more pressure in derived class
        virtual void computeOtherSources(){}


protected:

    // Protected Functions

        // is coupled
        bool isCoupledSolver() const;

public:

    //- Runtime type information
    TypeName("basePressure");


    // Constructors

        //- Construct from patch and internal field
        basePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        basePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given basePressureFvPatchScalarField
        //  onto a new patch
        basePressureFvPatchScalarField
        (
            const basePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        basePressureFvPatchScalarField
        (
            const basePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new basePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        basePressureFvPatchScalarField
        (
            const basePressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new basePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return p0
            virtual const scalarField& p0() const
            {
                return p0_;
            }

            //- Return the name of the velocity field
            virtual const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the density field
            //  to allow adjustment
            virtual word& UName()
            {
                return UName_;
            }

            //- Return the name of the density field
            virtual const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            virtual word& rhoName()
            {
                return rhoName_;
            }

            //- Return flux values at patch to check in-out flow
            virtual const scalarField& phiP() const;


        // Coefficients - these are overidden by the derived BCs

            //- C0 Field
            virtual tmp<scalarField> C0Field() const
            {
                return tmp<scalarField>(p0_);
            }

            //- C1 Field
            virtual tmp<scalarField> C1Field() const
            {
                return tmp<scalarField>(new scalarField(this->size(), Zero));
            }

            //- C2 Field
            virtual tmp<scalarField> C2Field() const
            {
                return tmp<scalarField>(new scalarField(this->size(), Zero));
            }

            //- C3 Field
            virtual tmp<scalarField> C3Field() const
            {
                return tmp<scalarField>(new scalarField(this->size(), Zero));
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientBoundaryCoeffs() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //-
            virtual void boundaryRelaxMatrix(fvBlockMatrix<vector>&) const
            {};


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
