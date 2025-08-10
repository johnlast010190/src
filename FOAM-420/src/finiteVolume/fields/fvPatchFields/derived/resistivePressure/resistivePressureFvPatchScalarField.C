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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/resistivePressure/resistivePressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::resistivePressureFvPatchScalarField::resistivePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    rAUName_("rAU"),
    p0_(p.size(), 0.0),
    C1_(0.0),
    C2_(0.0),
    Ct_(0.0)
{}


Foam::resistivePressureFvPatchScalarField::resistivePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rAUName_(dict.lookupOrDefault<word>("rAU", "(1|A(U))")),
    p0_("p0", dict, p.size()),
    C1_(dict.lookupOrDefault<scalar>("C1", 0.0)),
    C2_(dict.lookupOrDefault<scalar>("C2", 0.0)),
    Ct_(dict.lookupOrDefault<scalar>("Ct", 0.0))
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


Foam::resistivePressureFvPatchScalarField::resistivePressureFvPatchScalarField
(
    const resistivePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rAUName_(ptf.rAUName_),
    p0_(mapper(ptf.p0_)),
    C1_(ptf.C1_),
    C2_(ptf.C2_),
    Ct_(ptf.Ct_)
{}


Foam::resistivePressureFvPatchScalarField::resistivePressureFvPatchScalarField
(
    const resistivePressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rAUName_(ptf.rAUName_),
    p0_(ptf.p0_),
    C1_(ptf.C1_),
    C2_(ptf.C2_),
    Ct_(ptf.Ct_)
{}


Foam::resistivePressureFvPatchScalarField::resistivePressureFvPatchScalarField
(
    const resistivePressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rAUName_(ptf.rAUName_),
    p0_(ptf.p0_),
    C1_(ptf.C1_),
    C2_(ptf.C2_),
    Ct_(ptf.Ct_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::resistivePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(p0_, p0_);
}


void Foam::resistivePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const resistivePressureFvPatchScalarField& tiptf =
        refCast<const resistivePressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::resistivePressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(p0_, gAverage(p0_));
}


void Foam::resistivePressureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{

    if (updated())
    {
        return;
    }

    word rAUlookup = rAUName_;
    if (!db().foundObject<volScalarField>(rAUlookup))
    {
        rAUlookup = "(1|A(U))";

        WarningInFunction
            << "Could not find specified rUA field: " << rAUName_
            << ", defaulting to (1|A(U)) for this update."
            << endl;
    }

    const fvPatchField<scalar>& rAp =
        patch().lookupPatchFieldInDb<volScalarField, scalar>(this->db(), rAUlookup);

    const surfaceScalarField& phi =
        this->db().lookupObject<surfaceScalarField>(phiName_);

    scalarField Un(Up & patch().nf()); //(phip / patch().magSf());
    scalarField sqrUt( magSqr(Up - Un *patch().nf()) );

    scalarField& p = *this;

    scalarField pi(this->patchInternalField());

    // calculate base gradient before modifying pressure
    scalarField dp0(this->snGrad());

    const scalarField& rd = patch().deltaCoeffs();

    //scalarField dp1 = rd * (p - pi);

    scalarField dUn(Un * 0.0);

    //iterative pressure solver based on
    // Un_new = Un_old - 1/A * (grad(p) - grad(p_old))
    // and
    // p = p0 + (C1 + 0.5*C2*|Un|)Un + 0.5*Ct*sqr(Ut) :: incompressible
    // and
    // p = p0 + (rho*C1 + 0.5*rho*C2*|Un|)Un + 0.5*Ct*rho*sqr(Ut)
    // :: compressible

    if (phi.dimensions() == dimVelocity*dimArea)
    {

        forAll(p, i)
        {

            label n = 0;
            scalar err = GREAT;

            do
            {
                scalar pold = p[i];

                dUn[i] = - (p[i] - pi[i])*rAp[i]*rd[i] + dp0[i]*rAp[i];

                scalar UN = Un[i] + dUn[i];

                scalar f = p[i] - p0_[i] - (C1_ + 0.5*C2_*mag(UN))*UN
                    - 0.5*Ct_*sqrUt[i];

                scalar df
                    = 1 + 0.25 * C2_*UN * (2 * UN * rd[i] * rAp[i])
                    /max(mag(UN), VSMALL)
                    + C1_ * rAp[i] * rd[i]
                    + 0.5 * C2_ * rAp[i] * rd[i] * mag(UN);

                p[i] = pold - f/(sign(df)*max(mag(df), VSMALL));

                err = mag(pold - p[i])/max(SMALL, max(mag(p[i]), mag(pold)));

                n++;

                if (n == 100)
                {
                    WarningInFunction
                        << "Boundary pressure solver did not converge. error ="
                        << err << endl;
                }

            } while (err > 0.000001 && n < 100);

        }

        //convergence acceleration
        //correct outlet pressure and patch adjacent pressure using means
        //to accelerate convergence of boundary
        scalar pMean = gSum(p*patch().magSf()) / gSum(patch().magSf());
        scalar tpMean = gSum((p0_ + (C1_ + 0.5*C2_*mag(Un))*Un + 0.5*Ct_*sqrUt)
                *patch().magSf()) / gSum(patch().magSf());

        p += - 0.1 * (pMean - tpMean);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        // density for compressible
           const fvPatchField<scalar>& rhop =
              patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), rhoName_);

        forAll(p, i)
        {

            label n = 0;
            scalar err = GREAT;

            do
            {
                scalar pold = p[i];

                dUn[i] = - (p[i] - pi[i])*rAp[i]*rd[i] + dp0[i]*rAp[i];

                scalar UN = Un[i] + dUn[i];

                scalar f
                    = p[i] - p0_[i]
                    - (rhop[i] * C1_ + 0.5 * rhop[i] * C2_*mag(UN))*UN
                    - 0.5*rhop[i]*Ct_*sqrUt[i];

                scalar df
                    = 1 + 0.25 * rhop[i] * C2_ * UN * (2 * UN * rd[i] * rAp[i])
                    /max(mag(UN), VSMALL)
                    + rhop[i] * C1_ * rAp[i] * rd[i]
                    + 0.5 * rhop[i] * C2_ * rAp[i] * rd[i] * mag(UN);

                p[i] = pold - f/(sign(df)*max(mag(df), VSMALL));

                err = mag(pold - p[i])/max(SMALL, max(mag(p[i]), mag(pold)));

                n++;

                if (n == 100)
                {
                    WarningInFunction
                        << "Boundary pressure solver did not converge. error ="
                        << err << endl;
                }

            } while (err > 0.000001 && n < 100);
        }

        //convergence acceleration
        //correct outlet pressure using mean target value
        //to accelerate convergence of boundary
        scalar pMean = gSum(p*patch().magSf()) / gSum(patch().magSf());
        scalar tpMean = gSum((p0_ + (rhop*C1_ + 0.5*rhop*C2_*mag(Un))*Un
                + 0.5*rhop*Ct_*sqrUt)
                *patch().magSf()) / gSum(patch().magSf());

        p += - 0.1 * (pMean - tpMean);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are not correct"
            << abort(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();

}


void Foam::resistivePressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs(patch().lookupPatchFieldInDb<volVectorField, vector>(this->db(), UName_));
}


void Foam::resistivePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "rAU", "(1|A(U))", rAUName_);
    writeEntryIfDifferent<scalar>(os, "C1", 0.0, C1_);
    writeEntryIfDifferent<scalar>(os, "C2", 0.0, C2_);
    writeEntryIfDifferent<scalar>(os, "Ct", 0.0, Ct_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        resistivePressureFvPatchScalarField
    );
}

// ************************************************************************* //
