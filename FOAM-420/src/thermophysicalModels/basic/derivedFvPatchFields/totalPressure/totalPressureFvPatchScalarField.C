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
    (c) 2020-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/totalPressure/totalPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    psiName_("none"),
    gamma_(0.0),
    p0_(p.size(), 0.0),
    staticHead_(*this),
    freeStreamVelocity_(nullptr)
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? readScalar(dict.lookup("gamma")) : 1),
    p0_("p0", dict, p.size()),
    staticHead_(*this, p, dict),
    freeStreamVelocity_
    (
        dict.found("freeStreamVelocity")
        ? new point(dict.lookup("freeStreamVelocity"))
        : nullptr
    )
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(mapper(ptf.p0_)),
    staticHead_(*this, ptf.staticHead_, mapper),

    freeStreamVelocity_
    (
         ptf.freeStreamVelocity_.valid()
       ? new vector(ptf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const totalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(p0_, p0_);
    staticHead_.autoMap(m);
}


void Foam::totalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const totalPressureFvPatchScalarField& tiptf =
        refCast<const totalPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);

    staticHead_.rmap(tiptf.staticHead_, addr);
}


void Foam::totalPressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(p0_, gAverage(p0_));

    staticHead_.autoMapGIB(mapper);
}


void Foam::totalPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    tmp<scalarField> pp;

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            const fvPatchField<scalar>& rho =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(this->db(), staticHead_.getRhoName());

            pp = p0p - 0.5*rho*(1.0 - pos0(phip))*magSqr(Up);
        }
        else
        {
            // High-speed compressible flow

            const fvPatchField<scalar>& psip =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(this->db(), psiName_);

            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                pp =
                    (
                        p0p
                       /pow
                        (
                            (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                            1.0/gM1ByG
                        )
                    );
            }
            else
            {
                pp = p0p/(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up));
            }
        }

    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow
        pp = p0p - 0.5*neg(phip)*magSqr(Up);
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }


    staticHead_.addStaticHead(pp.ref());


    forceAssign(pp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::totalPressureFvPatchScalarField::updateCoeffs()
{
    if (freeStreamVelocity_.valid())
    {
        vectorField Up
        (
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName())
            - freeStreamVelocity_()
        );

        updateCoeffs(p0(), Up);
    }
    else
    {
        updateCoeffs
        (
            p0(),
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName())
        );
    }
}


void Foam::totalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
    p0_.writeEntry("p0", os);

    staticHead_.write(os);

    if (freeStreamVelocity_.valid())
    {
        os.writeEntry("freeStreamVelocity", freeStreamVelocity_());
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
