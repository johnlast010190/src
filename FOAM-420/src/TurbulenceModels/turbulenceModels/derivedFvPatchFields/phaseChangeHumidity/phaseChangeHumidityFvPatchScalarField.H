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
    Foam::phaseChangeHumidityFvPatchScalarField

Description
    Foam::phaseChangeHumidityFvPatchScalarField

SourceFiles
    phaseChangeHumidityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeHumidityFvPatchScalarField_H
#define phaseChangeHumidityFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class phaseChangeHumidityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeHumidityFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{

protected:

    // Protected data

        //- Temperature field name
        word TName_;

        //- pressure field name
        word PName_;

                //- density field name
        word rhoName_;

        //- Turbulent mass diffusion field name
        word DtName_;

        //- surface liquid mass / massflux (for steady state) (kg/m2 | kg/m2/s)
        scalarField filmMass_;

        //- heat flux of phase change
        scalarField heatFlux_;

        //- wet area ratio, =1 for films
        scalarField wetRatio_;

        //- steady / transient Switch
        Switch transient_;

        //- option to enable evaporation for steady state
        Switch steadyEvaporation_;

        //- option to enable thermal coupling
        Switch thermalCoupling_;

        // physical constants

            //- heat of phase change
            scalar dHv_;

            //- Molecular weight of vapour phase
            scalar Mvap_;

            //- Molecular weight of mixture
            scalar Mmix_;

        //- droplet data

            //- droplet modelling switch
            Switch dropSim_;

            //- last effective condensing droplet diameter
            autoPtr<scalarField> dDropCondensing_;

            //- advancing liquid contact angle
            scalar thetaA_;

            //- retreating liquid contact angle
            scalar thetaR_;

            //- fixed condensing wet/dry area ratio
            scalar dryWetCondensing_;

            //- liquid density
            scalar rhoLiquid_;



private:

    //- true for any transient, false for steady
    bool timeScheme();

    //- calculate mass diffusion coefficient
    tmp<scalarField> massDiffusionCoeff();

    //- lookup patch density
    tmp<scalarField> patchDensity();

    //- calculate effective face droplet diameter given contact angle
    scalar dropDiameter(label i, scalar theta);

    //- calculate effective face contact angle given diameter
    scalar contactAngle(label i, scalar d);

    //- update wet area ratio
    void updateWetDryAreaRatio();

public:

    //- Runtime scalar information
    TypeName("phaseChangeHumidity");


    // Constructors

    //- Construct from patch and internal field
    phaseChangeHumidityFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    phaseChangeHumidityFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given phaseChangeHumidityFvPatchScalarField onto a new patch
    phaseChangeHumidityFvPatchScalarField
    (
        const phaseChangeHumidityFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    phaseChangeHumidityFvPatchScalarField
    (
        const phaseChangeHumidityFvPatchScalarField&
    );

    //- Construct and return a clone

    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return tmp<fvPatchField<scalar>>
        (
            new phaseChangeHumidityFvPatchScalarField(*this)
        );
    }

    //- Construct as copy setting internal field reference
    phaseChangeHumidityFvPatchScalarField
    (
        const phaseChangeHumidityFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference

    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchField<scalar>>
        (
            new phaseChangeHumidityFvPatchScalarField(*this, iF)
        );
    }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


    // Return defining fields

        scalarField& filmMass()
        {
            return filmMass_;
        }

        const scalarField& filmMass() const
        {
            return filmMass_;
        }

        scalarField& heatFlux()
        {
            return heatFlux_;
        }

        const scalarField& heatFlux() const
        {
            return heatFlux_;
        }

        scalarField& wetAreaRatio()
        {
            return wetRatio_;
        }

        const scalarField& wetAreaRatio() const
        {
            return wetRatio_;
        }

        bool dropDvalid() const
        {
            return dDropCondensing_.valid();
        }

    //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
