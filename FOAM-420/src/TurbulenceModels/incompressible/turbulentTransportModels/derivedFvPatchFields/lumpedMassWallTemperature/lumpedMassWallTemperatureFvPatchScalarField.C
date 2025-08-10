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
    (c) 2020-2021 Esi Ltd.
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "lumpedMassWallTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Cp_(0.0),
    mass_(0.0),
    curTimeIndex_(-1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Cp_(ptf.Cp_),
    mass_(ptf.mass_),
    curTimeIndex_(-1)
{}


lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    Cp_(readScalar(dict.lookup("Cp"))),
    mass_(readScalar(dict.lookup("mass"))),
    curTimeIndex_(-1)
{
    refGrad() = 0.0;
    valueFraction() = 1.0;
    refValue() = scalarField("value", dict, p.size());

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lumpedMassWallTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }

    const label patchi = patch().index();

    // Calculate heat flux in or out the wall
    scalarField& Tp(*this);
    const scalarField& magSf = patch().magSf();

    const scalar deltaT(db().time().deltaTValue());

    const incompressibleTurbulenceModel& turbModel =
        db().lookupObject<incompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const tmp<volScalarField> talphaEff
        = turbModel.alphaEff()*turbModel.rho()*turbModel.Cp();
    const volScalarField& alphaEff = talphaEff;
    const scalarField& alphaEffp = alphaEff.boundaryField()[patchi];

    // qconv
    scalarField q(alphaEffp*snGrad());

    // Add fvOption sources
    scalarField bSourceDeriv;
    q -= boundarySources(Tp, bSourceDeriv);

    // Total heat in or out of the wall
    const scalar Qnet = gSum(q*magSf);

    Tp += -(Qnet/mass_/Cp_)*deltaT;

    refGrad() = 0.0;
    refValue() = Tp;
    valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qin(0);
        scalar Qout(0);

        forAll(q, facei)
        {
            if (q[facei] > 0.0) // out the wall
            {
                Qout += q[facei]*magSf[facei];
            }
            else if (q[facei] < 0.0) // into the wall
            {
                Qin += q[facei]*magSf[facei];
            }
        }

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Qnet
            << " wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << " Qin [W]:" << Qin
            << " Qout [W]:" << Qout
            << endl;
    }

    curTimeIndex_ = this->db().time().timeIndex();
}


void lumpedMassWallTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);

    os.writeEntry("Cp", Cp_);
    os.writeEntry("mass", mass_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    lumpedMassWallTemperatureFvPatchScalarField
);

// ************************************************************************* //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
