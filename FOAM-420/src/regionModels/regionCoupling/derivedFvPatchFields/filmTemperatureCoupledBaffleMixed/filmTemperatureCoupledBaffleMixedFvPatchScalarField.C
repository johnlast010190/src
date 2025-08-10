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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/filmTemperatureCoupledBaffleMixed/filmTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "surfaceFilmModel/surfaceFilmModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::filmModelType&
Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }

    DynamicList<word> modelNames;
    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        modelNames.append(iter()->regionMesh().name());
    }

    FatalErrorInFunction
        << "Unable to locate film region " << filmRegionName_
        << ".  Available regions include: " << modelNames
        << abort(FatalError);

    return **models.begin();
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::
filmTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    mappedPatchBase(p.patch()),
    boundaryKappa(db(), patch(), "undefined", "undefined", "undefined-K"),
    filmRegionName_("surfaceFilmProperties"),
    TnbrName_("undefined-Tnbr")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::
filmTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const filmTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    boundaryKappa(db(), patch(), ptf),
    filmRegionName_(ptf.filmRegionName_),
    TnbrName_(ptf.TnbrName_)
{}

Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::
filmTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    mappedPatchBase(p.patch(), dict),
    boundaryKappa(db(), patch(), dict),
    filmRegionName_
    (
           dict.lookup("filmRegion")
    ),
    TnbrName_(dict.lookup("Tnbr"))
{

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}

Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::
filmTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const filmTemperatureCoupledBaffleMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    mappedPatchBase(wtcsf.patch().patch(), wtcsf),
    boundaryKappa(db(), patch(), wtcsf),
    filmRegionName_(wtcsf.filmRegionName_),
    TnbrName_(wtcsf.TnbrName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Upcasting for cleaner code
    const mappedPatchBase& mpp = *this;
    const label nbrPatchI = mpp.samplePolyPatch().index();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    label patchI = patch().index();

    if (mpp.sampleRegion() == "region0")
    {
        patchI = nbrPatch.index();
    }

    // Retrieve film model
    const filmModelType& film = filmModel();
    const label& filmPatch = film.regionPatchID(patchI);

    scalarField alphaFilm = film.alpha().boundaryField()[filmPatch];
    scalarField TFilm = film.T().boundaryField()[filmPatch];
    scalarField deltaFilm = film.delta().boundaryField()[filmPatch];
    scalarField filmKappa = film.kappa().boundaryField()[filmPatch];

    filmTemperatureCoupledBaffleMixedFvPatchScalarField
        nbrField =
    refCast
    <
        const filmTemperatureCoupledBaffleMixedFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            TnbrName_
        )
    );


    if (mpp.sampleRegion() != "region0")
    {
        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mppb =
            refCast<const mappedPatchBase>(patch().patch());

        mppb.distribute(alphaFilm);
        mppb.distribute(TFilm);
        mppb.distribute(deltaFilm);
        mppb.distribute(filmKappa);

        tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
        tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

        nbrIntFld.ref() =  nbrField.patchInternalField();
        nbrKDelta.ref() =  nbrField.kappa()*nbrPatch.deltaCoeffs();

        mpp.distribute(nbrIntFld.ref());
        mpp.distribute(nbrKDelta.ref());

        tmp<scalarField> myKDelta = kappa()*patch().deltaCoeffs();

        nbrKDelta.ref() = 2.0*alphaFilm*filmKappa/(mag(deltaFilm)+ROOTVSMALL)
        + nbrKDelta.ref()*(1.0-alphaFilm);

        this->refValue() = nbrIntFld()*(1-alphaFilm)+alphaFilm*TFilm;
        this->refGrad() = 0.0;
        this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());
    }
    else
    {
        const mappedPatchBase& mppb =
            refCast<const mappedPatchBase>(nbrField.patch().patch());

        mppb.distribute(alphaFilm);
        mppb.distribute(TFilm);
        mppb.distribute(deltaFilm);
        mppb.distribute(filmKappa);

        tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
        tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

        nbrIntFld.ref() = alphaFilm*TFilm
        + nbrField.patchInternalField().ref()*(1.0-alphaFilm);

        nbrKDelta.ref() = 2.0*alphaFilm*filmKappa/(mag(deltaFilm)+ROOTVSMALL)
        + nbrField.kappa()*nbrPatch.deltaCoeffs()*(1.0-alphaFilm);

        mpp.distribute(nbrIntFld.ref());
        mpp.distribute(nbrKDelta.ref());

        tmp<scalarField> myKDelta = kappa()*patch().deltaCoeffs();

        this->refValue() = nbrIntFld();
        this->refGrad() = 0.0;
        this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());
    }

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::filmTemperatureCoupledBaffleMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    mappedPatchBase::write(os);
    os.writeEntry("filmRegion", filmRegionName_);
    os.writeEntry("Tnbr", TnbrName_);

    boundaryKappa::write(os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        filmTemperatureCoupledBaffleMixedFvPatchScalarField
    );
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
