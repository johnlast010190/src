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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/multiphaseGIBvelocityInlet/multiphaseGIBvelocityInletFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "global/unitConversion/unitConversion.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseGIBvelocityInletFvPatchVectorField::
multiphaseGIBvelocityInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fvPatchField<vector>(p, iF),
    injectVel_(vector::zero),
    angle_(1.0)
{}


Foam::multiphaseGIBvelocityInletFvPatchVectorField::
multiphaseGIBvelocityInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const vector& value
)
:
    fvPatchField<vector>(p, iF, value),
    injectVel_(vector::zero),
    angle_(1.0)
{}


Foam::multiphaseGIBvelocityInletFvPatchVectorField::
multiphaseGIBvelocityInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<vector>(p, iF, dict, valueRequired),
    injectVel_(dict.lookup("injectVel")),
    angle_(dict.lookupOrDefault<scalar>("angle", 1.0))
{
    computeBlendingFactor();
}


Foam::multiphaseGIBvelocityInletFvPatchVectorField::
multiphaseGIBvelocityInletFvPatchVectorField
(
    const multiphaseGIBvelocityInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<vector>(p, iF),
    injectVel_(ptf.injectVel_),
    angle_(ptf.angle_)
{
    if (ptf.blendingFactor_.size())
    {
        mapper(blendingFactor_, ptf.blendingFactor_);
    }
}


Foam::multiphaseGIBvelocityInletFvPatchVectorField::
multiphaseGIBvelocityInletFvPatchVectorField
(
    const multiphaseGIBvelocityInletFvPatchVectorField& pivpvf
)
:
    fvPatchField<vector>(pivpvf),
    blendingFactor_(pivpvf.blendingFactor_),
    injectVel_(pivpvf.injectVel_),
    angle_(pivpvf.angle_)
{}


Foam::multiphaseGIBvelocityInletFvPatchVectorField::
multiphaseGIBvelocityInletFvPatchVectorField
(
    const multiphaseGIBvelocityInletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fvPatchField<vector>(pivpvf, iF),
    blendingFactor_(pivpvf.blendingFactor_),
    injectVel_(pivpvf.injectVel_),
    angle_(pivpvf.angle_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseGIBvelocityInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<vector>::autoMap(m);

    if (blendingFactor_.size())
    {
        m(blendingFactor_, blendingFactor_);
    }
}


void Foam::multiphaseGIBvelocityInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fvPatchField<vector>::rmap(ptf, addr);

    if (blendingFactor_.size())
    {
        const multiphaseGIBvelocityInletFvPatchVectorField& tiptf =
            refCast<const multiphaseGIBvelocityInletFvPatchVectorField>(ptf);

        blendingFactor_.rmap(tiptf.blendingFactor_, addr);
    }
}


void Foam::multiphaseGIBvelocityInletFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fvPatchField<vector>::autoMapGIB(mapper);
    mapper.map(blendingFactor_, scalar(0));
}


void Foam::multiphaseGIBvelocityInletFvPatchVectorField::computeBlendingFactor()
{
    tmp<vectorField> nHat = this->patch().nf();

    blendingFactor_ = pos0(((-injectVel_/mag(injectVel_)) & nHat) - Foam::cos(degToRad(angle_)));
}


void Foam::multiphaseGIBvelocityInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    computeBlendingFactor();

    const fvMesh& mesh = internalField().mesh();

    vectorField::operator=
    (
        blendingFactor_*injectVel_
    );

    if (mesh.changing())
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        const scalar deltaT = mesh.time().deltaTValue();

        const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());

        scalarField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );

        const vectorField n(p.nf());
        const scalarField& magSf = p.magSf();
        tmp<scalarField> Un = phip/(magSf + VSMALL);


        vectorField::operator+=(Up + n*(Un - (n & Up)));
    }

    fvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::multiphaseGIBvelocityInletFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<vector>>
    (
        new Field<vector>(this->size(), Zero)
    );
}

Foam::tmp<Foam::vectorField>
Foam::multiphaseGIBvelocityInletFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}


Foam::tmp<Foam::vectorField>
Foam::multiphaseGIBvelocityInletFvPatchVectorField::gradientInternalCoeffs() const
{
    return -pTraits<vector>::one*this->patch().deltaCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::multiphaseGIBvelocityInletFvPatchVectorField::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


void Foam::multiphaseGIBvelocityInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("injectVel", injectVel_);
    os.writeEntry("angle", angle_);
    blendingFactor_.writeEntry("blendingFactor", os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::multiphaseGIBvelocityInletFvPatchVectorField::operator=
(
    const fvPatchField<vector>& ptf
)
{
    this->evaluate();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        multiphaseGIBvelocityInletFvPatchVectorField
    );
}

// ************************************************************************* //
