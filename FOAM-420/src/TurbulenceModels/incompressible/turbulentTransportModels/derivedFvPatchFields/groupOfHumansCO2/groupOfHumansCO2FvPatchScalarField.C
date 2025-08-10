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
    (c) 1991-2009 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "groupOfHumansCO2FvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentTransportModels/derivedFvPatchFields/groupOfHumansHeatFlux/groupOfHumansHeatFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

groupOfHumansCO2FvPatchScalarField::
groupOfHumansCO2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(p, iF),
    RCO2_(188.9)
{}


groupOfHumansCO2FvPatchScalarField::
groupOfHumansCO2FvPatchScalarField
(
    const groupOfHumansCO2FvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(ptf, p, iF, mapper),
    RCO2_(ptf.RCO2_)
{}


groupOfHumansCO2FvPatchScalarField::
groupOfHumansCO2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(p, iF, dict),
    RCO2_(dict.lookupOrDefault<scalar>("RCO2", 188.9))
{}


groupOfHumansCO2FvPatchScalarField::
groupOfHumansCO2FvPatchScalarField
(
    const groupOfHumansCO2FvPatchScalarField& ptf
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(ptf),
    RCO2_(ptf.RCO2_)
{}


groupOfHumansCO2FvPatchScalarField::
groupOfHumansCO2FvPatchScalarField
(
    const groupOfHumansCO2FvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(ptf, iF),
    RCO2_(ptf.RCO2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void groupOfHumansCO2FvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const groupOfHumansHeatFluxFvPatchScalarField& Tp
        = refCast<const groupOfHumansHeatFluxFvPatchScalarField>
        (
            patch().lookupPatchField<volScalarField, scalar>("T")
        );

    scalar rhoCO2 = Tp.Plocal() / Tp.Tlocal() / RCO2_;

    turbulentMassFluxConcentrationFvPatchScalarField::m_
        = ( Tp.Ntotal() * rhoCO2 / 1000.0 / 3600.0 / Tp.Areal() )
        * 0.83 * Tp.Meff() * Tp.Aeff() * Tp.Troom() / 5.617 / 273.15;

    turbulentMassFluxConcentrationFvPatchScalarField::updateCoeffs();
}


void groupOfHumansCO2FvPatchScalarField::write
(
    Ostream& os
) const
{
    writeEntryIfDifferent<scalar>(os, "RCO2", 188.9, RCO2_);
    turbulentMassFluxConcentrationFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    groupOfHumansCO2FvPatchScalarField
);

// Backward compat.
addNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    groupOfHumansCO2FvPatchScalarField,
    groupOfHumansCO2
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam


// ************************************************************************* //
