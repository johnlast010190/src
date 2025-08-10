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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

Class
    Foam::regionCoupledEnergyFluxFvPatchScalarField

Description
    Energy region coupled implicit boundary condition based on energy equation
    flux equality. Supports boundary sources, and optional thin thermal
    layer resistance. Can be used for T, h or e variables.
    For further details, see regionCoupledFluxFvPatchScalarField.

SourceFiles
    regionCoupledEnergyFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledEnergyFluxFvPatchScalarField_H
#define regionCoupledEnergyFluxFvPatchScalarField_H

#include "regionCoupledFluxFvPatchScalarField.H"
#include "fvMesh/fvPatches/derived/mapped/mappedWallFvPatch.H"
#include "matrices/LduMatrix/LduMatrix/LduInterfaceField/LduInterfaceField.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class regionCoupledEnergyFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledEnergyFluxFvPatchScalarField
:
    public regionCoupledFluxFvPatchScalarField
{
private:

    // Private data

    //- My thermo
    mutable const basicThermo* thermoPtr_;

    //- Multiplier and offset for the field such that
    // (fieldValue/pMultiplier_-pOffset_) is matched at the interface
    scalarField pMultiplier_;
    scalarField pOffset_;
    //- Multiplier for the flux such that fluxScale_*flux is matched at the
    //  at the interface
    scalarField fluxScale_;

    // Private functions

        //- Set thermo pointer
        void setThermo() const;

        //- Detect neighbour field name based on info from thermo
        void setNeighbourFieldName() const;

protected:

    // Protected functions

        //- Overridable functions to scale and offset such that
        //  fieldValue = fieldScale*(matchedValue+fieldOffset)
        //  where matchedValue is continuous across the interface
        virtual tmp<scalarField> fieldScale(const tmp<scalarField>& f) const override
        {
            return f*pMultiplier_;
        }

        virtual tmp<scalarField>
        invFieldScale(const tmp<scalarField>& f) const override
        {
            return f/pMultiplier_;
        }

        virtual tmp<scalarField>
        fieldOffset(const tmp<scalarField>& f) const override
        {
            return f + pOffset_;
        }

        virtual tmp<scalarField>
        negFieldOffset(const tmp<scalarField>& f) const override
        {
            return f - pOffset_;
        }


public:

    //- Runtime type information
    TypeName("regionCoupledEnergyFlux");


    // Constructors

    //- Construct from patch and internal field
    regionCoupledEnergyFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    regionCoupledEnergyFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given regionCoupledEnergyFluxFvPatchScalarField
    // onto a new patch
    regionCoupledEnergyFluxFvPatchScalarField
    (
        const regionCoupledEnergyFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    regionCoupledEnergyFluxFvPatchScalarField
    (
        const regionCoupledEnergyFluxFvPatchScalarField&
    );

    //- Construct and return a clone
    tmp<fvPatchField<scalar>> clone() const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new regionCoupledEnergyFluxFvPatchScalarField(*this)
            );
    }

    //- Construct as copy setting internal field reference
    regionCoupledEnergyFluxFvPatchScalarField
    (
        const regionCoupledEnergyFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new regionCoupledEnergyFluxFvPatchScalarField(*this, iF)
            );
    }


    //- Destructor
    virtual ~regionCoupledEnergyFluxFvPatchScalarField()
    {}


    //Access functions

        //- Return name of other side field that is coupled
        virtual word coupledField() const override
        {
            if (neighbourFieldName_ == word::null)
            {
                setNeighbourFieldName();
            }
            return regionCoupledFluxFvPatchScalarField::coupledField();
        }


    // Evaluation functions

        //- Is this patch on the h/e field rather than temperature?
        bool ishe() const
        {
            setThermo();
            return internalField().name() == thermoPtr_->he().name();
        }

        void updateCoeffs() override;

        virtual tmp<scalarField>
        boundarySources(const scalarField& pf, scalarField& df) const override;

    //- Write
    virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
