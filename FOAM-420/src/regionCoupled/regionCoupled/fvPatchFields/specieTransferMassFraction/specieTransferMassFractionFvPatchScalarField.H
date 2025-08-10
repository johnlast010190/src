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
    (c) 2023 Esi Ltd.
    (c) 2019-2023 OpenFOAM Foundation

Class
    Foam::specieTransferMassFractionFvPatchScalarField

Description
    Abstract base class for specie-transferring mass fraction boundary
    conditions.

    Derived classes compute the actual species flux. Those fluxes are summed up
    by the specieTransferVelocityFvPatchVectorField velocity condition to
    generate the net mass transfer into or out of the patch. This mass-fraction
    condition then generates a corrective diffusive flux to ensure that the
    correct amounts of the transferring species are transported through the
    patch.

SourceFiles
    specieTransferMassFractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef specieTransferMassFractionFvPatchScalarField_H
#define specieTransferMassFractionFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;

/*---------------------------------------------------------------------------*\
        Class specieTransferMassFractionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class specieTransferMassFractionFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
public:

        //- Enumeration for the property driving the transfer
        enum property
        {
            massFraction,
            moleFraction,
            molarConcentration,
            partialPressure
        };

        //- Property type names
        static const NamedEnum<property, 4> propertyNames_;


private:

    // Private Data

        //- Name of the flux field
        const word phiName_;

        //- Name of the velocity field
        const word UName_;

        //- Flux of this species
        mutable scalarField phiYp_;

        //- Current time index
        mutable label timeIndex_;

        //- Clip mass fractions?
        Switch clipMassFractions_;


protected:

    // Protected Data

        //- Transfer coefficient
        const scalar c_;

        //- Property driving the transfer
        const property property_;


public:


    // Constructors

        //- Construct from patch and internal field
        specieTransferMassFractionFvPatchScalarField
        (
            const fvPatch& ptf,
            const DimensionedField<scalar, volMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        specieTransferMassFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  specieTransferMassFractionFvPatchScalarField
        //  onto a new patch
        specieTransferMassFractionFvPatchScalarField
        (
            const specieTransferMassFractionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor setting internal field reference
        specieTransferMassFractionFvPatchScalarField
        (
            const specieTransferMassFractionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions


        // Mapping functions

            //- Return the flux of this species
            virtual tmp<scalarField> calcPhiYp() const = 0;

            //- Return the flux of this species
            virtual const scalarField& phiYp() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) override;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;


    // Member Operators

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
