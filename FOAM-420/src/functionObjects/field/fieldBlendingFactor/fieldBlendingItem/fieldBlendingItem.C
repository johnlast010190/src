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
    (c) 2011 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldBlendingFactor/fieldBlendingItem/fieldBlendingItem.H"
#include "finiteVolume/fvc/fvcAverage.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fieldBlendingItem, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fieldBlendingItem::blendingFactorTransfer(scalar& currentBF, scalar newBF)
{
    if (newBF < currentBF) currentBF = newBF;
    else currentBF = min(currentBF + stab_, newBF);

    currentBF = min(1.0, max(0.0, currentBF));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fieldBlendingItem::fieldBlendingItem
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict,
    const word blendingFieldName
)
:
    obr_(obr),
    mesh_(mesh),
    fieldName_(blendingFieldName),
    stab_(1.0/max(1.0, dict.lookupOrDefault<scalar>("stabilise", 1))),
    blendSources_()
{
    readBlendingItem(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldBlendingItem::~fieldBlendingItem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fieldBlendingItem::readBlendingItem
(
    const dictionary& dict
)
{
    // create blending sources
    PtrList<entry> sourceDicts(dict.lookup("blendingSources"));
    blendSources_.setSize(sourceDicts.size());

    forAll(sourceDicts, sI)
    {
        const entry& sourceDict = sourceDicts[sI];
        blendSources_.set
        (
            sI,
            blendingSource::New(obr_, mesh_, sourceDict.dict())
        );
    }

    if (!obr_.foundObject<surfaceScalarField>(fieldName_))
    {
        //create and check in blending factor field
        autoPtr<surfaceScalarField> fieldBlendingFactor;

        if (debug)
        {
            fieldBlendingFactor.reset
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        fieldName_,
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(fieldName_, dimless, 1.0)
                )
            );
        }
        else
        {
            fieldBlendingFactor.reset
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        fieldName_,
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(fieldName_, dimless, 1.0)
                )
            );
        }

        obr_.store(fieldBlendingFactor);
    }

}

void fieldBlendingItem::update(const surfaceScalarField* bsbfPtr)
{
    surfaceScalarField& blendingFactorField =
        const_cast<surfaceScalarField&>
        (
            obr_.lookupObject<surfaceScalarField>(fieldName_)
        );
    // Avoid recalculation if field has been loaded from disk
    if
    (
        mesh_.time().timeIndex() == mesh_.time().startTimeIndex()
     && exists(blendingFactorField.objectPath())
    )
    {
        return;
    }

    Info<< tab << fieldName_ << " blended faces:" << flush;

    // create intermediate field from base and blending sources
    autoPtr<surfaceScalarField> newBFF
    (
        new surfaceScalarField
        (
            IOobject
            (
                "newBF",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(fieldName_, dimless, 1.0)
        )
    );

    if (bsbfPtr != nullptr)
    {
        newBFF().forceAssign(*bsbfPtr);
    }

    scalarField& internalNBFF = newBFF();
    surfaceScalarField::Boundary& boundaryNBFF = newBFF->boundaryFieldRef();

    const labelUList& owner = mesh_.owner();

    forAll(blendSources_, bsI)
    {
        tmp<surfaceScalarField> bsiField(blendSources_[bsI].blendingFactor());

        forAll(owner, facei)
        {
            internalNBFF[facei] = min(bsiField()[facei], internalNBFF[facei]);
        }

        forAll(newBFF->boundaryField(), patchi)
        {
            //if (newBFF->boundaryField()[patchi].coupled())
            {
                forAll(newBFF->boundaryField()[patchi], facei)
                {
                    boundaryNBFF[patchi][facei] =
                        min(newBFF->boundaryField()[patchi][facei],
                        bsiField().boundaryField()[patchi][facei]);
                }
            }
        }
    }

    // apply new blending factor to bf field, limit blending factors
    // and stabilise
    scalarField& internalBFF = blendingFactorField;
    surfaceScalarField::Boundary& blendingFactorbf =
        blendingFactorField.boundaryFieldRef();

    forAll(owner, facei)
    {
        blendingFactorTransfer(internalBFF[facei], internalNBFF[facei]);
    }

    forAll(blendingFactorField.boundaryField(), patchi)
    {
        //if (blendingFactorField.boundaryField()[patchi].coupled())
        {
            forAll(blendingFactorField.boundaryField()[patchi], facei)
            {
                blendingFactorTransfer
                (
                    blendingFactorbf[patchi][facei],
                    newBFF->boundaryField()[patchi][facei]
                );
            }
        }
    }

    Info<< endl;
}

void fieldBlendingItem::write()
{
    forAll(blendSources_, bsI)
    {
        blendSources_[bsI].write();
    }

    this->field().write();
}

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
