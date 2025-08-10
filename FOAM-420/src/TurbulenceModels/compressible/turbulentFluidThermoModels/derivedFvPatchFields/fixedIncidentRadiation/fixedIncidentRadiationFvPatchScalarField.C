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
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/fixedIncidentRadiation/fixedIncidentRadiationFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "global/constants/constants.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "submodels/absorptionEmissionModel/absorptionEmissionModel/absorptionEmissionModel.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), "undefined", "undefined", "undefined-K"),
    qrIncident_(p.size(), 0.0)
{}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(psf, p, iF, mapper),
    boundaryKappa(db(), patch(), psf),
    qrIncident_(psf.qrIncident_)
{}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), dict),
    qrIncident_("qrIncident", dict, p.size())
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(psf, iF),
    boundaryKappa(db(), patch(), psf),
    qrIncident_(psf.qrIncident_)
{}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    boundaryKappa(db(), patch(), ptf),
    qrIncident_(ptf.qrIncident_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    m(qrIncident_, qrIncident_);
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(psf, addr);

    const fixedIncidentRadiationFvPatchScalarField& thftpsf =
        refCast<const fixedIncidentRadiationFvPatchScalarField>
        (
            psf
        );

    qrIncident_.rmap(thftpsf.qrIncident_, addr);
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedGradientFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(qrIncident_, scalar(0));
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField intFld(patchInternalField());

    const radiation::radiationModel& radiation =
        db().lookupObject<radiation::radiationModel>("radiationProperties");

    scalarField emissivity
    (
        radiation.absorptionEmission().e()().boundaryField()[patch().index()]
    );

    gradient() =
        emissivity
       *(
            qrIncident_
          - physicoChemical::sigma.value()*pow4(*this)
        )/kappa();

    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar qr = gSum(kappa()*gradient()*patch().magSf());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " -> "
            << " radiativeFlux:" << qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    boundaryKappa::write(os);
    qrIncident_.writeEntry("qrIncident", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedIncidentRadiationFvPatchScalarField
    );
}
}

// ************************************************************************* //
