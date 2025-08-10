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
    (c) 2010-2020 Esi Ltd.

Class
    Foam::coupledWallEnergyFvPatchScalarField

Description
    Energy region coupled implicit boundary condition supporting optional
    radiation and thin thermal layer resistance through the thicknessLayers
    and kappaLayers entries. Can be used for T, h or e variables.
    The fvPatch is treated as uncoupled from the delta point of view, as
    all regions are not always available for interpolations etc.
    In the mesh the fvPatch is an interface and is incorporated
    into the matrix implicitly.
    The thermal conductivity, \c kappa, can either be retrieved from the
    mesh database using the \c lookup option, or from a \c solidThermo
    or \c fluidThermo thermophysical package.
    \heading Patch usage
    \table
        Property     | Description             | Required    | Default value
        kappaMethod  | thermal conductivity option | yes     |
        kappaName    | name of thermal conductivity field | no  |
        Qr         | name of the radiative flux in this region or "none" | no | Qr
        QrNbr      | name of the radiative flux in the nbr region or "none" | no | Qr
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivites per layer [W/m/K] | no |
    \endtable
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            compressible::coupledWallEnergy;
        kappaMethod     lookup;
        KappaName       kappa;
        Qr              Qr; // or none. Name of Qr field on local region
        QrNbr           Qr; // or none. Name of Qr field on neighbour region
        kappaLayers     (10 0.1)
        thicknessLayers (0.1 0.01);
        value           uniform 300;
    }
    \endverbatim
    Needs to be on underlying mappedWallFvPatch.
    Note: kappa : heat conduction at patch. Gets supplied how to
            lookup/calculate
     kappaMethod:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'fluidThermo' : use fluidThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use solidThermo kappa()
    - 'directionalSolidThermo' directionalKappa()

SourceFiles
    coupledWallEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledWallEnergyFvPatchScalarField_H
#define coupledWallEnergyFvPatchScalarField_H

#include "fvMesh/fvPatches/derived/mapped/mappedWallFvPatch.H"
#include "matrices/LduMatrix/LduMatrix/LduInterfaceField/LduInterfaceField.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.H"
#include "primitives/enums/NamedEnum.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
             Class coupledWallEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledWallEnergyFvPatchScalarField
:
    public regionCoupledFvPatchField<scalar>,
    public boundaryKappa
{
private:

    // Private data

    //- AutoPtr to nbr thermo
    mutable const basicThermo* nbrThermoPtr_;

    //- AutoPtr to my thermo
    mutable const basicThermo* thermoPtr_;

    //- Owner cell weight
    scalarField lWeights_;

    //- Explicit correction to face value (for radiation and humidity)
    scalarField fCorr_;

    //- Previous value of fCorr (at last call to evaluate)
    scalarField fCorrEval_;

    //- Store computed values for convenience
    scalarField pCpEffByNbrCpEff_;

    //- Conductivity of layers
    scalarList kappaLayers_;

    //- Thickness of layers
    scalarList thicknessLayers_;

    //- Inverse of contact conductance per unit area (computed from above)
    scalar contactInsulance_;

    //- Radiation field name or "none"
    word QrName_;

    //- Neighbouring radiation field name or "none"
    word QrNbrName_;

    //- net patch heat source (used for damping)
    scalarField Qsum_;

    //- heat source under relaxation factor
    scalar qURF_;


    // Private functions

    //- Local weight for this coupled field
    tmp<scalarField> weights() const;

    //- Return nbr temperature internal field
    tmp<scalarField> patchNeighbourTemperatureField() const;

    //- Return local temperature internal field
    tmp<scalarField> patchInternalTemperatureField() const;

    //- Set thermo pointers
    void setThermos() const;

public:

    //- Runtime type information
    TypeName("compressible::coupledWallEnergy");


    // Constructors

    //- Construct from patch and internal field
    coupledWallEnergyFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    coupledWallEnergyFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given coupledWallEnergyFvPatchScalarField
    // onto a new patch
    coupledWallEnergyFvPatchScalarField
    (
        const coupledWallEnergyFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    coupledWallEnergyFvPatchScalarField
    (
        const coupledWallEnergyFvPatchScalarField&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<scalar>> clone() const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new coupledWallEnergyFvPatchScalarField(*this)
            );
    }

    //- Construct as copy setting internal field reference
    coupledWallEnergyFvPatchScalarField
    (
        const coupledWallEnergyFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new coupledWallEnergyFvPatchScalarField(*this, iF)
            );
    }


    //- Destructor
    virtual ~coupledWallEnergyFvPatchScalarField()
    {}


    // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        ) override;

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        ) override;

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;

    // Evaluation functions

        //- Is this patch on the h/e field rather than temperature?
        virtual bool ishe() const
        {
            setThermos();
            return internalField().name() == thermoPtr_->he().name();
        }

        //- Return patch-normal gradient
        virtual tmp<scalarField> snGrad() const override;

        //- Return patch-normal gradient
        //  Note: the deltaCoeffs supplied are not used
        virtual tmp<scalarField> snGrad
        (
            const scalarField& deltaCoeffs
        ) const override;


        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        ) override;

        virtual void updateCoeffs() override;

        //- Return the matrix coefficients of the neighbouring cell value
        //  corresponding to the evaluation of the value of this patchField
        virtual tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField
        virtual tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientInternalCoeffs() const override;

        //- Return the matrix coefficients of the neighbouring cell value
        //  corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const override;

        //- Convert boundary/internal coeffs to region-coupled boundary/internal
        //  coeffs (from coeffs of boundary face value and internal value to
        //  coeffs of neighbour & internal value)
        void regionCoupledBoundaryCoeffs
        (
            const fvMatrix<scalar>& matrix,
            const scalarField& bouCoeffs,
            const scalarField& intCoeffs,
            scalarField& coupledBouCoeffs,
            scalarField& coupledIntCoeffs
        ) override;

        //- Return explicit correction to new boundary face value
        tmp<scalarField> faceCorr() const override
        {
            return tmp<scalarField>(fCorr_);
        }

        //- Return explicit correction to current boundary face value
        tmp<scalarField> faceCorrEval() const override
        {
            return tmp<scalarField>(fCorrEval_);
        }

        virtual tmp<scalarField> boundarySources
        (
            const scalarField& pf, scalarField& df
        ) const override;

    //- Write
    virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
