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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/outflowVelocity/outflowVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
outflowVelocityFvPatchVectorField::
outflowVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    relax_(0.2),
    scale_(1.0),
    massCorr_(0.0)
{}


Foam::
outflowVelocityFvPatchVectorField::
outflowVelocityFvPatchVectorField
(
    const outflowVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    relax_(ptf.relax_),
    scale_(ptf.scale_),
    massCorr_(ptf.massCorr_)
{
}


Foam::
outflowVelocityFvPatchVectorField::
outflowVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.2)),
    scale_(dict.lookupOrDefault<scalar>("scale", 1.0)),
    massCorr_(dict.lookupOrDefault<scalar>("massCorr", 0.0))
{
}


Foam::
outflowVelocityFvPatchVectorField::
outflowVelocityFvPatchVectorField
(
    const outflowVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    relax_(ptf.relax_),
    scale_(ptf.scale_),
    massCorr_(ptf.massCorr_)
{
}


Foam::
outflowVelocityFvPatchVectorField::
outflowVelocityFvPatchVectorField
(
    const outflowVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    relax_(ptf.relax_),
    scale_(ptf.scale_),
    massCorr_(ptf.massCorr_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outflowVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    const vectorField Uin(this->patchInternalField());

    if (scale_ == 1.0 && massCorr_ ==0.0) forceAssign(Uin);

    const vectorField Unn((Uin)&patch().nf()*patch().nf());
    const vectorField Ut(*this - Unn);

    if (massCorr_ == 0)
    {
        forceAssign(Ut + Unn*scale_);
    }
    else
    {
        //- just in first interation if the adjustable flux is zero
        forceAssign(
            Ut + Unn + 0.1*(massCorr_)*patch().nf()
        );
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::outflowVelocityFvPatchVectorField::scaleProfile
(
    const scalar scale,
    const scalar mass
)
{
    scale_ = scale;
    massCorr_ = mass;
    scale_ = max(0.1, scale_);
    scale_ = min(10, scale_);
}


void Foam::outflowVelocityFvPatchVectorField::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>& bEq
) const
{
    bEq.boundaryRelax(relax_, this->patch().index());
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::outflowVelocityFvPatchVectorField
::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<Field<vector>> vict
    (
        new Field<vector>(this->size(), 1.0)
    );

    return vict;
}

Foam::tmp<Foam::Field<Foam::vector>>
Foam::outflowVelocityFvPatchVectorField
::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<Field<vector>> vbct
    (
        new Field<vector>(this->size(), Zero)
    );

    return vbct;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::outflowVelocityFvPatchVectorField
::valueDivInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<Field<vector>> vict
    (
        new Field<vector>(this->size(), Zero)
    );

    return vict;
}

Foam::tmp<Foam::Field<Foam::vector>>
Foam::outflowVelocityFvPatchVectorField
::valueDivBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<Field<vector>> vbct
    (
        new Field<vector>(*this)
    );

    return vbct;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::outflowVelocityFvPatchVectorField::gradientInternalCoeffs() const
{
    return tmp<Field<vector>>
    (
        new Field<vector>(this->size(), Zero)
    );
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::outflowVelocityFvPatchVectorField::gradientBoundaryCoeffs() const
{
    return tmp<Field<vector>>
    (
        new Field<vector>(this->size(), Zero)
    );
}



void Foam::outflowVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    writeEntryIfDifferent<scalar>(os, "relax", 0.2, relax_);

    writeEntryIfDifferent<scalar>(os, "scale", 1, scale_);

    writeEntryIfDifferent<scalar>(os, "massCorr", 0, massCorr_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       outflowVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
