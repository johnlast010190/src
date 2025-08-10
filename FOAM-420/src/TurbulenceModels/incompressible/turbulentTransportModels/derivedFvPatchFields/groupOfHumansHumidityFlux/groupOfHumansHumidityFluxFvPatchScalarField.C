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

#include "groupOfHumansHumidityFluxFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentTransportModels/derivedFvPatchFields/groupOfHumansHeatFlux/groupOfHumansHeatFluxFvPatchScalarField.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

groupOfHumansHumidityFluxFvPatchScalarField::
groupOfHumansHumidityFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(p, iF)
{}


groupOfHumansHumidityFluxFvPatchScalarField::
groupOfHumansHumidityFluxFvPatchScalarField
(
    const groupOfHumansHumidityFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(ptf, p, iF, mapper)
{}


groupOfHumansHumidityFluxFvPatchScalarField::
groupOfHumansHumidityFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(p, iF, dict)
{}


groupOfHumansHumidityFluxFvPatchScalarField::
groupOfHumansHumidityFluxFvPatchScalarField
(
    const groupOfHumansHumidityFluxFvPatchScalarField& ptf
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(ptf)
{}


groupOfHumansHumidityFluxFvPatchScalarField::
groupOfHumansHumidityFluxFvPatchScalarField
(
    const groupOfHumansHumidityFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentMassFluxConcentrationFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void groupOfHumansHumidityFluxFvPatchScalarField::updateCoeffs()
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

    scalar TroomDegC = Tp.Troom() - 273.15;
    turbulentMassFluxConcentrationFvPatchScalarField::m_
        = ( Tp.Ntotal() / 1000.0 / 3600.0 / Tp.Areal() )
        * (0.4206 * sqr(TroomDegC) -14.27 * TroomDegC + 153.99);

    turbulentMassFluxConcentrationFvPatchScalarField::updateCoeffs();
}


void groupOfHumansHumidityFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    turbulentMassFluxConcentrationFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    groupOfHumansHumidityFluxFvPatchScalarField
);

// Backward compat.
addNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    groupOfHumansHumidityFluxFvPatchScalarField,
    groupOfHumansHumidityFlux
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam


// ************************************************************************* //
