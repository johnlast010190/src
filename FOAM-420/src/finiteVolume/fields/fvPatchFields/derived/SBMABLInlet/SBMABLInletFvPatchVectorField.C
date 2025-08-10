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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/SBMABLInlet/SBMABLInletFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    void SBMABLInletFvPatchVectorField::scaleFlux
(
    const scalar targetFlowRate,
    scalarField& phiPN
) const
{
    scalar Atotal = gSum(patch().magSf());
    const scalar currentFlowRate = gSum(phiPN);

    scalar addFlux(targetFlowRate - currentFlowRate);
    phiPN += addFlux*patch().magSf()/Atotal;
}

scalar SBMABLInletFvPatchVectorField::calcFlowRate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    //volumetric or mass
    scalar flowRate = gSum(phip);

    return flowRate;
}

tmp<scalarField> SBMABLInletFvPatchVectorField::scaleNormVelocity
(
    const scalar flowRate
) const
{

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    tmp<scalarField> phiPN(new scalarField(this->size(), 0.0));

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        phiPN.ref() = (patch().Sf() & (*this));
        scaleFlux(targetFlowRate_, phiPN.ref());
        phiPN.ref() /= patch().magSf();
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(),
                rhoName_
            );

        phiPN.ref() = rhop*(patch().Sf() & (*this));
        scaleFlux(targetFlowRate_, phiPN.ref());
        phiPN.ref() /= (rhop*patch().magSf());
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    // return velocity
    return phiPN;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SBMABLInletFvPatchVectorField::
SBMABLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    decayingTurbulenceFvPatchVectorField(p, iF),
    profile_(),
    fluctuations_(false),
    writeVortons_(false),
    scaleU_(false),
    updateTFR_(true),
    phiName_("phi"),
    rhoName_("rho"),
    targetFlowRate_()
{}


SBMABLInletFvPatchVectorField::
SBMABLInletFvPatchVectorField
(
    const SBMABLInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    decayingTurbulenceFvPatchVectorField(ptf, p, iF, mapper),
    profile_(ptf.profile_->clone()),
    fluctuations_(ptf.fluctuations_),
    writeVortons_(ptf.writeVortons_),
    scaleU_(ptf.scaleU_),
    updateTFR_(ptf.updateTFR_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    targetFlowRate_(ptf.targetFlowRate_)
{}


SBMABLInletFvPatchVectorField::
SBMABLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    decayingTurbulenceFvPatchVectorField(p, iF),
    profile_(ABLProfile::New(p, dict)),
    fluctuations_(dict.lookupOrDefault<Switch>("fluctuations", true)),
    writeVortons_(dict.lookupOrDefault<Switch>("writeVortons", false)),
    scaleU_(dict.lookupOrDefault<Switch>("scaleVelocity", false)),
    updateTFR_(true),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    targetFlowRate_(-1)
{
    //set U,R,L fields in decay turbulence base class
    profile_->setProfiles(U_, L_, R_);

    // initialise base class data mambers
    decayingTurbulenceFvPatchVectorField::initialiseFromDictionary(dict);

}


SBMABLInletFvPatchVectorField::
SBMABLInletFvPatchVectorField
(
    const SBMABLInletFvPatchVectorField& blpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    decayingTurbulenceFvPatchVectorField(blpvf, iF),
    profile_(blpvf.profile_->clone()),
    fluctuations_(blpvf.fluctuations_),
    writeVortons_(blpvf.writeVortons_),
    scaleU_(blpvf.scaleU_),
    updateTFR_(blpvf.updateTFR_),
    phiName_(blpvf.phiName_),
    rhoName_(blpvf.rhoName_),
    targetFlowRate_(blpvf.targetFlowRate_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void SBMABLInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (fluctuations_)
    {
        decayingTurbulenceFvPatchVectorField::updateCoeffs();
    }

    if (scaleU_)
    {
        // calculate target flow rate at the first loop
        if (updateTFR_)
        {
            targetFlowRate_ = calcFlowRate();
            Info<<"Target flow rate calculation: "<<targetFlowRate_<<endl;

            if (targetFlowRate_ < 0.0)
            {
                updateTFR_  = false;
            }
        }

        //scale only after the targed is known
        if (!updateTFR_)
        {
            // keep the tangential component
            vectorField Ut((*this)- ((*this)&patch().nf())*patch().nf());
            // scale the normal component
            tmp<scalarField> Upn(scaleNormVelocity(targetFlowRate_));

            Field<vector>::operator=(Upn()*patch().nf()+Ut);
        }
    }
}


void SBMABLInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    //base class input settings
    os.writeEntry("minVortonLength", minVortonLength_);

    writeEntryIfDifferent<label>(os, "direction", 1, direction_);
    writeEntryIfDifferent<Switch>(os, "writeVortons", false, writeVortons_);
    writeEntryIfDifferent<Switch>(os, "scaleVelocity", false, scaleU_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<Switch>(os, "fluctuations", true, fluctuations_);

    //profile input
    profile_->write(os);

    //fields
    L_.writeEntry("L", os);
    U_.writeEntry("U", os);
    R_.writeEntry("R", os);

    //Vortons
    if (writeVortons_)
    {
        if (Pstream::master())
        {
            os.writeEntry("vortons", vortons_);
        }
        //instantaneous R field only useful if structures are available
        Rinst_.writeEntry("Rinst", os);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    SBMABLInletFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
