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
    (c) 2020 Rudolf Hellmuth - Vascular Flow Technologies
    (c) 2020 Esi Ltd.
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/windkesselPressure/windkesselPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WindkesselPressureFvPatchScalarField::
WindkesselPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(0.0),
    Rp_(0.0),
    Rd_(0.0),
    C_(0.0),
    Qold_(0.0),
    pOld_(0.0),
    Qdummy_(Qold_),
    curTimeIndex_(-1),
    alpha_(1.0)
{
}


Foam::WindkesselPressureFvPatchScalarField::
WindkesselPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(dict.lookupOrDefault<scalar>("p0", 0.0)),
    Rp_(dict.lookupOrDefault<scalar>("Rp", 0.0)),
    Rd_(readScalar(dict.lookup("Rd"))),
    C_(readScalar(dict.lookup("C"))),
    Qold_(dict.lookupOrDefault<scalar>("Qold", 0.0)),
    pOld_(dict.lookupOrDefault<scalar>("pOld", 0.0)),
    Qdummy_(Qold_),
    curTimeIndex_(-1),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0))
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
        fvPatchScalarField::operator=
        (
            scalarField(patchInternalField())
        );
    }

    //if pOld is not specified in the p dictionary, use the current patch
    //average value as pOld.
    pOld_ = dict.lookupOrDefault<scalar>
    (
        "pOld",
        patchAvg(*this)
    );
}


Foam::WindkesselPressureFvPatchScalarField::
WindkesselPressureFvPatchScalarField
(
    const WindkesselPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_),
    Rp_(ptf.Rp_),
    Rd_(ptf.Rd_),
    C_(ptf.C_),
    Qold_(ptf.Qold_),
    pOld_(ptf.pOld_),
    Qdummy_(Qold_),
    curTimeIndex_(-1),
    alpha_(ptf.alpha_)
{}


Foam::WindkesselPressureFvPatchScalarField::
WindkesselPressureFvPatchScalarField
(
    const WindkesselPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    p0_(ptf.p0_),
    Rp_(ptf.Rp_),
    Rd_(ptf.Rd_),
    C_(ptf.C_),
    Qold_(ptf.Qold_),
    pOld_(ptf.pOld_),
    Qdummy_(Qold_),
    curTimeIndex_(-1),
    alpha_(ptf.alpha_)
{}


Foam::WindkesselPressureFvPatchScalarField::
WindkesselPressureFvPatchScalarField
(
    const WindkesselPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    p0_(ptf.p0_),
    Rp_(ptf.Rp_),
    Rd_(ptf.Rd_),
    C_(ptf.C_),
    Qold_(ptf.Qold_),
    pOld_(ptf.pOld_),
    Qdummy_(Qold_),
    curTimeIndex_(-1),
    alpha_(ptf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::WindkesselPressureFvPatchScalarField::patchAvg(scalarField& var)
{
    return gSum(patch().magSf()*var)/gSum(patch().magSf());
}


Foam::scalar Foam::WindkesselPressureFvPatchScalarField::pLumpedModel(scalar Q)
{
    scalar dt = this->db().time().deltaTValue();

    // If this is the first inner step inside a time step:
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // Patch current pressure
        scalarField& p = *this;

        // Patch average pressure
        pOld_ = patchAvg(p);

        // Flow rate of the previous time step stored in the dummy variable.
        Qold_ = Qdummy_;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    Qdummy_ = Q;

    if (debug)
    {
        Info<< "* " << patch().name() << ":  "
             << "Qold = " << Qold_ << ", "
             << "pOld - p0 = " << pOld_ - p0_ << endl;
    }

    // Calculate equivalent pressure to current flow rate
    //return p0_ + (Rd_*C_*(meanpOld - p0_) + Rd_*dt*Q)/(dt + Rd_*C_);
    return
    (
          p0_
        + (
              Rd_*C_*(pOld_ - p0_)
            + ((Rd_ + Rp_)*dt + Rp_*Rd_*C_)*Q
            - Rp_*Rd_*C_*Qold_
          )
        / (dt + Rd_*C_)
    );
}

void Foam::WindkesselPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), "U");

    const vectorField& Sf = patch().Sf();

    // Patch normal velocity
    // Obs.: the patch field is used to take advantage of relaxation factors of
    // the SIMPLE algorithm.
    scalarField Upn( (Up & Sf)/patch().magSf() );

    // Calculate current flow rate
    scalar Q = gSum(Upn*patch().magSf());

    // Calculate dynamic head for inflow
    scalarField Kin( alpha_*neg(Upn)*0.5*mag(Upn)*Upn );
    scalar Kmean = patchAvg(Kin);

    // Calculate equivalent pressure to current flow rate
    scalar p = pLumpedModel(Q);

    if (debug)
    {
        Info<< "* " << patch().name() << ":  "
             << "Q = " << Q << ", "
             << "p - p0 = " << p - p0_ << endl;
    }

    fixedValueFvPatchScalarField::forceAssign(p - Kmean + Kin);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::WindkesselPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<scalar>(os, "p0", 0.0, p0_);
    writeEntryIfDifferent<scalar>(os, "Rp", 0.0, Rp_);
    os.writeEntry("Rd", Rd_);
    os.writeEntry("C", C_);
    os.writeEntry("Qold", Qold_);
    os.writeEntry("pOld", pOld_);
    writeEntryIfDifferent<scalar>(os, "alpha", 1.0, alpha_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        WindkesselPressureFvPatchScalarField
    );
}

// ************************************************************************* //
