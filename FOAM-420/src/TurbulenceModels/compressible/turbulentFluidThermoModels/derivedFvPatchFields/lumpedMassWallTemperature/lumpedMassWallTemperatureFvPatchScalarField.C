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
    (c) 2016 OpenCFD Ltd.
    (c) 2021 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/lumpedMassWallTemperature/lumpedMassWallTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), "undefined", "undefined", "undefined-K"),
    Cp_(0.0),
    mass_(0.0),
    curTimeIndex_(-1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    boundaryKappa(db(), patch(), ptf),
    Cp_(ptf.Cp_),
    mass_(ptf.mass_),
    curTimeIndex_(-1)
{}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    boundaryKappa(db(), patch(), dict),
    Cp_(readScalar(dict.lookup("Cp"))),
    mass_(readScalar(dict.lookup("mass"))),
    curTimeIndex_(-1)
{
    refGrad() = 0.0;
    valueFraction() = 1.0;
    refValue() = scalarField("value", dict, p.size());

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    boundaryKappa(tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    boundaryKappa(db(), patch(), tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedMassWallTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }

    // Calculate heat flux in or out the wall
    scalarField& Tp(*this);
    const scalarField& magSf = patch().magSf();

    const scalar deltaT(db().time().deltaTValue());

    //this is what kappa() returns
    scalarField CpAlphaEff( kappa() );
    scalarField q(CpAlphaEff*snGrad());

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


void Foam::lumpedMassWallTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    boundaryKappa::write(os);

    os.writeEntry("Cp", Cp_);
    os.writeEntry("mass", mass_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        lumpedMassWallTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
