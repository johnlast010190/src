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
    (c) 2018-2023 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/freestreamVelocity/freestreamVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//
Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF)
{
}

Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF, dict, false)
{
    referenceFrameFvPatch<vector>::read(dict);
    inletFlux_ = true;
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    freestreamValue() = vectorField("freestreamValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(freestreamValue());
    }

    refGrad() = Zero;
    if (dict.found("valueFraction"))
    {
        valueFraction()  = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        valueFraction() = 1;
    }
}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& wbppsf
)
:
    mixedFvPatchVectorField(wbppsf)
{}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& wbppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    Field<vector> Up(0.5*(patchInternalField() + *this));
    if (coorFramePtr())
    {
        if (coorFramePtr()->isAttachToMRF(this->patch().index()))
        {
            //- this will not work in MRF frame.
            //  Currently this BC will not be included in MRF faces since the
            //  farfield BCs will contain an inletFlux true flag. But you have
            //  to know if U is absolute/relative to the MRF frame, otherwise
            //  the Un-Ut work that follows will be wrong.
            //  Needs extension of the coordinate frames for this.
            const Field<vector> refVel(getFrameVelocity()());
            Up -= refVel;
        }
    }

    const Field<scalar> magUp(mag(Up));

    const Field<vector> nf(patch().nf());

    Field<scalar>& vf = valueFraction();

    forAll(vf, i)
    {
        if (magUp[i] > VSMALL)
        {
            vf[i] = 0.5 - 0.5*(Up[i] & nf[i])/magUp[i];
        }
        else
        {
            vf[i] = 0.5;
        }
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::freestreamVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    freestreamValue().writeEntry("freestreamValue", os);
    valueFraction().writeEntry("valueFraction", os);
    referenceFrameFvPatch<vector>::write(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        freestreamVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
