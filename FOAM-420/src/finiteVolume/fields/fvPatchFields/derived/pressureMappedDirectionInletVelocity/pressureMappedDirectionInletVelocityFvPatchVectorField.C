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

#include "fields/fvPatchFields/derived/pressureMappedDirectionInletVelocity/pressureMappedDirectionInletVelocityFvPatchVectorField.H"
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

pressureMappedDirectionInletVelocityFvPatchVectorField::
pressureMappedDirectionInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletDir_(p.size())
{}


pressureMappedDirectionInletVelocityFvPatchVectorField::
pressureMappedDirectionInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

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


pressureMappedDirectionInletVelocityFvPatchVectorField::
pressureMappedDirectionInletVelocityFvPatchVectorField
(
    const pressureMappedDirectionInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    inletDir_(mapper(ptf.inletDir_))
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



pressureMappedDirectionInletVelocityFvPatchVectorField::
pressureMappedDirectionInletVelocityFvPatchVectorField
(
    const pressureMappedDirectionInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    inletDir_(pivpvf.inletDir_)
{}


pressureMappedDirectionInletVelocityFvPatchVectorField::
pressureMappedDirectionInletVelocityFvPatchVectorField
(
    const pressureMappedDirectionInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureMappedDirectionInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    m(*this, *this);
    m(inletDir_, inletDir_);
}


void pressureMappedDirectionInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const pressureMappedDirectionInletVelocityFvPatchVectorField& tiptf =
        refCast
        <const pressureMappedDirectionInletVelocityFvPatchVectorField>
        (ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void Foam::pressureMappedDirectionInletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMapGIB(mapper);
    mapper.map(inletDir_, vector::zero);
}


// Update the coefficients associated with the patch field
void pressureMappedDirectionInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    // Get the scheduling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
       pressureMappedDirectionInletVelocityFvPatchVectorField::patch().patch()
    );
    const mapDistribute& distMap = mpp.map();
    // Force recalculation of schedule
    (void)distMap.schedule();

    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());

    vectorField mappedVelocity(this->size());

    switch (mpp.mode())
    {
       case mappedPatchBase::NEARESICELL:
       {
          if (mpp.sameRegion())
          {
             mappedVelocity = this->internalField();
          }
          else
          {
             mappedVelocity = nbrMesh.lookupObject<volVectorField>("U").internalField();
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

    inletDir_ = mappedVelocity/mag(mappedVelocity);

    //apply pressure inlet
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>("phi");

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    vectorField n( patch().nf() );
    scalarField ndmagS( (n & inletDir_)*patch().magSf() );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        forceAssign(inletDir_*phip/ndmagS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), "rho");

        forceAssign(inletDir_*phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void pressureMappedDirectionInletVelocityFvPatchVectorField::write
(Ostream& os) const
{
    fvPatchVectorField::write(os);
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureMappedDirectionInletVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
