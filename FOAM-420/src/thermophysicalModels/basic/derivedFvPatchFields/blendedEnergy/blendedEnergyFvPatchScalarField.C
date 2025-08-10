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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blendedEnergyFvPatchScalarField::blendedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    blendedFvPatchField<scalar>(p, iF)
{
    if (iF.member() == "h" || iF.member() == "e")
    {
        const basicThermo& thermo =
            iF.mesh().lookupObject<basicThermo>
            (
                iF.groupName(basicThermo::dictName, iF.group())
            );
        const fvPatchScalarField& Tpf = thermo.T().boundaryField()[p.index()];
        const blendedFvPatchField<scalar>& Tbpf =
            refCast<const blendedFvPatchField<scalar>>(Tpf);
        valueFraction_ = Tbpf.valueFraction();
        sensorName_ = Tbpf.sensorName();
        bcValueController_ = Tbpf.bcValueController().clone();
        boundaryOne_.set
        (
            fvPatchField<scalar>::New
            (
                basicThermo::heBoundaryType(Tbpf.boundaryOne()),
                basicThermo::heBoundaryBaseType(Tbpf.boundaryOne()),
                p,
                iF
            ).ptr()
        );
        boundaryTwo_.set
        (
            fvPatchField<scalar>::New
            (
                basicThermo::heBoundaryType(Tbpf.boundaryTwo()),
                basicThermo::heBoundaryBaseType(Tbpf.boundaryTwo()),
                p,
                iF
            ).ptr()
        );
    }
}


Foam::blendedEnergyFvPatchScalarField::blendedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    blendedFvPatchField<scalar>(p, iF, dict)
{}


Foam::blendedEnergyFvPatchScalarField::blendedEnergyFvPatchScalarField
(
    const blendedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    blendedFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::blendedEnergyFvPatchScalarField::blendedEnergyFvPatchScalarField
(
    const blendedEnergyFvPatchScalarField& ptf
)
:
    blendedFvPatchField<scalar>(ptf)
{}


Foam::blendedEnergyFvPatchScalarField::blendedEnergyFvPatchScalarField
(
    const blendedEnergyFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    blendedFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       blendedEnergyFvPatchScalarField
   );
}

// ************************************************************************* //
