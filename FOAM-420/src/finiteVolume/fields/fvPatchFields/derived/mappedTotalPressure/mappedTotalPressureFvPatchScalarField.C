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
    (c) 1991-2007 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedTotalPressure/mappedTotalPressureFvPatchScalarField.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mappedTotalPressureFvPatchScalarField::
mappedTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("undefined"),
    phiName_("undefined"),
    rhoName_("undefined"),
    psiName_("undefined"),
    gamma_(0.0),
    p0_(p.size(), 0.0),
    t_(-1)
 {}


mappedTotalPressureFvPatchScalarField::
mappedTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookup("U")),
    phiName_(dict.lookup("phi")),
    rhoName_(dict.lookup("rho")),
    psiName_(dict.lookup("psi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    p0_("p0", dict, p.size()),
    t_(-1)
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

    if (!isType<mappedFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


mappedTotalPressureFvPatchScalarField::
mappedTotalPressureFvPatchScalarField
(
    const mappedTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(mapper(ptf.p0_)),
    t_(-1)
{
    if (!isType<mappedFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}



mappedTotalPressureFvPatchScalarField::
mappedTotalPressureFvPatchScalarField
(
    const mappedTotalPressureFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    UName_(pivpvf.UName_),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    psiName_(pivpvf.psiName_),
    gamma_(pivpvf.gamma_),
    p0_(pivpvf.p0_),
    t_(-1)
{}


mappedTotalPressureFvPatchScalarField::
mappedTotalPressureFvPatchScalarField
(
    const mappedTotalPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    UName_(pivpvf.UName_),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    psiName_(pivpvf.psiName_),
    gamma_(pivpvf.gamma_),
    p0_(pivpvf.p0_),
    t_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mappedTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    m(*this, *this);
    m(p0_, p0_);
}


void mappedTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const mappedTotalPressureFvPatchScalarField& tiptf =
        refCast
        <const mappedTotalPressureFvPatchScalarField>
        (ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void mappedTotalPressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(p0_, gAverage(p0_));
}


// Update the coefficients associated with the patch field
void mappedTotalPressureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{
    if (this->updated())
    {
        return;
    }

    //only update pressure once per timestep
    if (t_ != this->db().time().value())
    {
        t_ = this->db().time().value();

        // Get the scheduling information from the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
           mappedTotalPressureFvPatchScalarField::patch().patch()
        );
        const mapDistribute& distMap = mpp.map();

        // Force recalculation of schedule
        (void)distMap.schedule();

        const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
        //map pressure
        //====================================================================//
        scalarField mappedPressure(this->size());
        switch (mpp.mode())
        {
           case mappedPatchBase::NEARESICELL:
           {
              if (mpp.sameRegion())
              {
                 mappedPressure = this->internalField();
              }
              else
              {
                 mappedPressure = nbrMesh.lookupObject<volScalarField>("p").internalField();
              }
              mapDistributeBase::distribute
              (
                 Pstream::defaultCommsType,
                 distMap.schedule(),
                 distMap.constructSize(),
                 distMap.subMap(),
                 false,
                 distMap.constructMap(),
                 false,
                 mappedPressure,
                 flipOp()
              );

              break;
           }
           case mappedPatchBase::NEARESIPATCHFACE:
           {
              const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
              (
                 mpp.samplePatch()
              );
              if (nbrPatchID < 0)
              {
                  FatalErrorInFunction
                      << "Unable to find sample patch " << mpp.samplePatch()
                      << " in region " << mpp.sampleRegion()
                      << " for patch " << this->patch().name() << nl
                      << abort(FatalError);
              }

              const volScalarField& nbrField = nbrMesh.lookupObject<volScalarField>("p");

              mappedPressure = nbrField.boundaryField()[nbrPatchID];
              mapDistributeBase::distribute
              (
                 Pstream::defaultCommsType,
                 distMap.schedule(),
                 distMap.constructSize(),
                 distMap.subMap(),
                 false,
                 distMap.constructMap(),
                 false,
                 mappedPressure,
                 flipOp()
              );

              break;
           }
           case mappedPatchBase::NEARESIFACE:
           {
              scalarField allValues(nbrMesh.nFaces(), pTraits<scalar>::zero);

              const volScalarField& nbrField = nbrMesh.lookupObject<volScalarField>("p");
              forAll(nbrField.boundaryField(), patchI)
              {
                 const fvPatchField<scalar>& pf =
                    nbrField.boundaryField()[patchI];
                 label faceStart = pf.patch().patch().start();

                 forAll(pf, faceI)
                 {
                    allValues[faceStart++] = pf[faceI];
                 }
              }

              mapDistributeBase::distribute
              (
                 Pstream::defaultCommsType,
                 distMap.schedule(),
                 distMap.constructSize(),
                 distMap.subMap(),
                 false,
                 distMap.constructMap(),
                 false,
                 allValues,
                 flipOp()
              );

              mappedPressure = this->patch().patchSlice(allValues);

              break;
           }
           default:
           {
               FatalErrorInFunction
                   << "Unknown sampling mode: " << mpp.mode()
                   << nl << abort(FatalError);
           }
        }

        //map velocity
        //====================================================================//
        vectorField mappedVelocity(this->size());

        const volVectorField& U
            = patch().boundaryMesh().mesh().lookupObject<volVectorField>("U");

        switch (mpp.mode())
        {
           case mappedPatchBase::NEARESICELL:
           {
              if (mpp.sameRegion())
              {
                 mappedVelocity = U.internalField();
              }
              else
              {
                 mappedVelocity =  nbrMesh.lookupObject<volVectorField>("U").internalField();
              }
              mapDistributeBase::distribute
              (
                 Pstream::defaultCommsType,
                 distMap.schedule(),
                 distMap.constructSize(),
                 distMap.subMap(),
                 false,
                 distMap.constructMap(),
                 false,
                 mappedVelocity,
                 flipOp()
              );

              break;
           }
           case mappedPatchBase::NEARESIPATCHFACE:
           {
              const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
              (
                 mpp.samplePatch()
              );
              if (nbrPatchID < 0)
              {
                 FatalErrorInFunction
                     << "Unable to find sample patch " << mpp.samplePatch()
                     << " in region " << mpp.sampleRegion()
                     << " for patch " << this->patch().name() << nl
                     << abort(FatalError);
              }


              const volVectorField& nbrField = nbrMesh.lookupObject<volVectorField>("U");

              mappedVelocity = nbrField.boundaryField()[nbrPatchID];
              mapDistributeBase::distribute
              (
                 Pstream::defaultCommsType,
                 distMap.schedule(),
                 distMap.constructSize(),
                 distMap.subMap(),
                 false,
                 distMap.constructMap(),
                 false,
                 mappedVelocity,
                 flipOp()
              );

              break;
           }
           case mappedPatchBase::NEARESIFACE:
           {
              vectorField allValues(nbrMesh.nFaces(), pTraits<vector>::zero);

              const volVectorField& nbrField = nbrMesh.lookupObject<volVectorField>("U");
              forAll(nbrField.boundaryField(), patchI)
              {
                 const fvPatchField<vector>& pf =
                    nbrField.boundaryField()[patchI];
                 label faceStart = pf.patch().patch().start();

                 forAll(pf, faceI)
                 {
                    allValues[faceStart++] = pf[faceI];
                 }
              }

              mapDistributeBase::distribute
              (
                 Pstream::defaultCommsType,
                 distMap.schedule(),
                 distMap.constructSize(),
                 distMap.subMap(),
                 false,
                 distMap.constructMap(),
                 false,
                 allValues,
                 flipOp()
              );

              mappedVelocity = this->patch().patchSlice(allValues);

              break;
           }
           default:
           {
               FatalErrorInFunction
               << "Unknown sampling mode: " << mpp.mode()
               << nl << abort(FatalError);
           }
        }

        //calc pressure from mapped pressure totalPressure and flux
        scalar magPatch = gSum(patch().magSf());

        scalar pMapMean = gSum(mappedPressure* patch().magSf())/magPatch;

        scalar U2MapMean
            = gSum(magSqr(mappedVelocity)*patch().magSf())/magPatch;

        scalarField p0Map
        (
            mappedPressure + 0.5*magSqr(mappedVelocity)
            + (p0_ - pMapMean - 0.5*U2MapMean)
        );

        const fvsPatchField<scalar>& phip =
            patch().lookupPatchFieldInDb<surfaceScalarField, scalar>(db(), phiName_);

        if (psiName_ == "none" && rhoName_ == "none")
        {
            forceAssign(p0Map - 0.5*(1.0 - pos0(phip))*magSqr(Up));
        }
        else if (rhoName_ == "none")
        {
            const fvPatchField<scalar>& psip =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), psiName_);

            if (gamma_ > 1.0)
            {
                scalar gM1ByG = (gamma_ - 1.0)/gamma_;

                forceAssign
                (
                    p0Map
                   /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    )
                );
            }
            else
            {
                forceAssign(p0Map/(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up)));
            }
        }
        else if (psiName_ == "none")
        {
            const fvPatchField<scalar>& rho =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), rhoName_);

            forceAssign(p0Map - 0.5*rho*(1.0 - pos0(phip))*magSqr(Up));
        }
        else
        {
            FatalErrorInFunction
                << " rho or psi set inconsitently, rho = " << rhoName_
                << ", psi = " << psiName_ << '.' << nl
                << "    Set either rho or psi or neither depending on the "
                   "definition of total pressure." << nl
                << "    Set the unused variables to 'none'."
                << "\n    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }

        //some reporting
        scalar pMean = gSum(*this * patch().magSf())/magPatch;

        scalar p0Mean = pMean + 0.5 * gSum(magSqr(Up)*patch().magSf())/magPatch;

        scalar UMean = -gSum(phip)/magPatch;

        Info<< "Inlet velocity : " << UMean << ", pressure: " << pMean
             << ", total pressure: " << p0Mean << endl;

    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void mappedTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
        patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName_)
    );
}


// Write
void mappedTotalPressureFvPatchScalarField::write
(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("U", UName_);
    os.writeEntry("phi", phiName_);
    os.writeEntry("rho", rhoName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mappedTotalPressureFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
