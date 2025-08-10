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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::resistivePressureFvPatchScalarField

Description
    Foam::resistivePressureFvPatchScalarField

    Resistance pressure inlet/outlet boundary incorporating
    a Darcy-Forchheimer type loss coefficient. Incorporates the
    effect of various types of downstream geometries and other resistances
    on the pressure at the outlet. Uses an iterative solver to predict
    change in velocity due to change in pressure.

    The pressure at the boundary is generally described as

    p = p0 + (C1   + 0.5 * C2 * mag(Un))*Un + 0.5*Ct*(Ut)^2

    where velocity Un can be found from

    Un_new = Un_old - 1/A * (grad(p) - grad(p_old))

    here

    p0 is the known downstream pressure
    C1 and C2 are loss coefficients for dp/dx*L = rho*C1*U + 0.5*rho*C2*U^2
    Ct is a coefficient between 0 and 1 depending on whether tangential velocity
    should be accounted for. Ct=0 --> full slip, Ct=1 --> no-slip
    Un is the surface normal velocity - positive outward and negative in
    Ut tangential velocity magnitude


    \heading Patch usage

    \table
        Property     | Description             | Required   | Default value
        value        | Restart value           | yes        |
        p0           | the known downstream pressure | yes  |
        C1           | linear loss coefficients | no        | 0
        C2           | quadratic loss coefficients | no     | 0
        Ct           | tangential loss coefficients | no    | 0

        U            | Name of velocity field   | no        | U
        phi          | Name of flux field       | no        | phi
        rho          | Name of density field    | no        | rho
        rAU          | Name of matrix diagonal field | no   | (1|A(U))
    \endtable



    Example of the boundary condition specification:
    \verbatim
        myPatch
        {
            type            resistivePressure;
            p0                uniform 0;
            C1              1.0;
            C2              0.0;
            Ct              0.0;
            U               U;
            phi             phi;
            rho             rho; // Not needed if phi volumetric
            rAU             rAU

        }

    \endverbatim

Note


SourceFiles
    resistivePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef resistivePressureFvPatchScalarField_H
#define resistivePressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class resistivePressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class resistivePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- name of the matrix diagonal inverse field
        word rAUName_;

        //- Far side pressure
        scalarField p0_;

        //- Linear loss coefficient
        scalar C1_;

        //- Quadratic loss coefficient
        scalar C2_;

        //- Tangential loss coefficient.
        scalar Ct_;

public:

    //- Runtime type information
    TypeName("resistivePressure");


    // Constructors

        //- Construct from patch and internal field
        resistivePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        resistivePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given resistivePressureFvPatchScalarField
        //  onto a new patch
        resistivePressureFvPatchScalarField
        (
            const resistivePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        resistivePressureFvPatchScalarField
        (
            const resistivePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new resistivePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        resistivePressureFvPatchScalarField
        (
            const resistivePressureFvPatchScalarField&,
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
                new resistivePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the far field pressure
            const scalarField& p0() const
            {
                return p0_;
            }

            //- Return reference to the far field pressure to allow adjustment
            scalarField& p0()
            {
                return p0_;
            }

            //- Return the linear loss coefficient
            scalar C1() const
            {
                return C1_;
            }

            //- Return reference to the linear loss coefficient
            scalar& C1()
            {
                return C1_;
            }

            //- Return the quadratic loss coefficient
            scalar C2() const
            {
                return C2_;
            }

            //- Return reference to the tangential loss coefficient
            scalar& C2()
            {
                return C2_;
            }

            //- Return the tangential loss coefficient
            scalar Ct() const
            {
                return Ct_;
            }

            //- Return reference to the quadratic loss coefficient
            scalar& Ct()
            {
                return Ct_;
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

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  using the given patch velocity field
            virtual void updateCoeffs(const vectorField& Up);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
