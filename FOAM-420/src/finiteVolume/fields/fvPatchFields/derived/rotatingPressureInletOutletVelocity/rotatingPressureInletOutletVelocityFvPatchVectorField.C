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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/rotatingPressureInletOutletVelocity/rotatingPressureInletOutletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
calcTangentialVelocity()
{
    const scalar t = this->db().time().timeOutputValue();
    vector om = omega_->value(t);

    vector axisHat = om/mag(om);
    const vectorField tangentialVelocity
    (
        (-om) ^ (patch().Cf() - axisHat*(axisHat & patch().Cf()))
    );

    const vectorField n(patch().nf());
    refValue() = tangentialVelocity - n*(n & tangentialVelocity);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletOutletVelocityFvPatchVectorField(p, iF),
    omega_()
{}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const rotatingPressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    pressureInletOutletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_, false)
{
    calcTangentialVelocity();
}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    pressureInletOutletVelocityFvPatchVectorField(p, iF, dict),
    omega_(Function1<vector>::New("omega", dict))
{
    patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    calcTangentialVelocity();
}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const rotatingPressureInletOutletVelocityFvPatchVectorField& rppvf
)
:
    pressureInletOutletVelocityFvPatchVectorField(rppvf),
    omega_(rppvf.omega_, false)
{
    calcTangentialVelocity();
}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const rotatingPressureInletOutletVelocityFvPatchVectorField& rppvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletOutletVelocityFvPatchVectorField(rppvf, iF),
    omega_(rppvf.omega_, false)
{
    calcTangentialVelocity();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("phi", phiName());
    omega_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rotatingPressureInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
