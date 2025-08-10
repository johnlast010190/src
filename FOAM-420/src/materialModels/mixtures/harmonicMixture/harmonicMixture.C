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
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "harmonicMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(harmonicMixture, 0);
    addToRunTimeSelectionTable(materialModel, harmonicMixture, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::volumeMassFractions& Foam::harmonicMixture::lookupOrConstructBase()
{
    const word name = groupName("massVolumeFractions", phaseName_);
    const dictionary& dict = materialsDict();
    if (!obr_.foundObject<volumeMassFractions>(name))
    {
        if (dict.found("phases"))
        {
            obr_.store(new phaseVolumeFractions(dict, obr_));
        }
        else
        {
            obr_.store(new speciesMassFractions(dict, obr_, phaseName_));
        }
    }
    return obr_.lookupObject<volumeMassFractions>(name);
}


Foam::harmonicMixture::harmonicMixture
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    frac_(lookupOrConstructBase())
{}


Foam::autoPtr<Foam::harmonicMixture> Foam::harmonicMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<harmonicMixture>
    (
        new harmonicMixture
        (
            obr,
            dict,
            phaseName,
            specieName,
            name
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::harmonicMixture::updateTable(const word& modelName)
{
    const wordList names =
        dict_->found("species")
      ? dict_->lookup<wordList>("species")
      : dict_->lookup<wordList>("phases");
    const label nModels = names.size();
    const word ttName = materialTables_.tableName(phaseName_, word::null);

    const HashTable<matScalarTable>& sModels = materialTables_.sTable();
    word funcName = IOobject::group();

    if (IOobject::name().find("Model") != string::npos)
    {
        funcName = string(this->name()).replace("Model", "");
    }
    const auto i = funcName.find('.');
    if (i != 0)
    {
        funcName = funcName.substr(0, i);
    }
    if (materialTables_.foundModel<scalar>(ttName, funcName))
    {
        sMod_.setSize(nModels);
        dep_.setSize(1);
        dep_[0].model = sModels[ttName][funcName];
        dep_[0].dependencies.setSize(nModels);
        forAll(names, I)
        {
            const word tName =
                dict_->found("species")
              ? materialTables_.tableName(phaseName_, names[I])
              : materialTables_.tableName(names[I], word::null);
            sMod_.set(I, sModels[tName][funcName]);
            dep_[0].dependencies.set(I, sModels[tName][funcName]);
        }
    }
}


Foam::tmp<Foam::volScalarField>
Foam::harmonicMixture::harmonicMixtureGeometric() const
{
    const dimensionedScalar s("s", sMod_[0].dimensions(), VSMALL);
    tmp<volScalarField> tMixture
    (
        new volScalarField
        (
            IOobject
            (
                sMod_[0].funcType(),
                frac_.fractions()[0].mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            frac_.fractions()[0].mesh(),
            dimensionedScalar
            (
                sMod_[0].funcType(),
                dimless/sMod_[0].dimensions(),
                0
            )
        )
    );
    volScalarField& mixture = tMixture.ref();
    volScalarField fractionsSum
    (
        IOobject
        (
            "fractionsSum",
            frac_.fractions()[0].mesh().time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        frac_.fractions()[0].mesh(),
        dimensionedScalar("fractionsSum", dimless, 0)
    );
    forAll(sMod_, i)
    {
        if (frac_.active()[i])
        {
            fractionsSum += frac_.fractions()[i];
            mixture += frac_.fractions()[i]/(sMod_[i]() + s);
        }
    }
    return fractionsSum/tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::harmonicMixture::harmonicMixtureInternal() const
{
    const label fieldSize(frac_.fractions()[0].primitiveField().size());
    tmp<scalarField> tMixture(new scalarField(fieldSize, 0));
    scalarField& mixture = tMixture.ref();
    scalarField fractionsSum(mixture.size(), 0);
    forAll(sMod_, i)
    {
        if (frac_.active()[i])
        {
            fractionsSum += frac_.fractions()[i].primitiveField();
            mixture +=
                frac_.fractions()[i].primitiveField()
               /(sMod_[i].primitiveField() + VSMALL);
        }
    }
    mixture = fractionsSum/mixture;
    return tMixture;
}


Foam::tmp<Foam::scalarField>
Foam::harmonicMixture::harmonicMixturePatch
(
    const label patchi
) const
{
    const label fieldSize(frac_.fractions()[0].boundaryField()[patchi].size());
    tmp<scalarField> tMixture(new scalarField(fieldSize, 0));
    scalarField& mixture = tMixture.ref();
    scalarField fractionsSum(mixture.size(), 0);
    forAll(sMod_, i)
    {
        if (frac_.active()[i])
        {
            fractionsSum += frac_.fractions()[i].boundaryField()[patchi];
            mixture +=
                frac_.fractions()[i].boundaryField()[patchi]
               /(sMod_[i].boundaryField()[patchi] + VSMALL);
        }
    }
    mixture = fractionsSum/mixture;
    return tMixture;
}


Foam::scalar Foam::harmonicMixture::harmonicMixtureCell
(
    const label celli
) const
{
    scalar mixture = 0;
    scalar fractionsSum = 0;
    forAll(sMod_, i)
    {
        if (frac_.active()[i])
        {
            fractionsSum += frac_.fractions()[i][celli];
            mixture += frac_.fractions()[i][celli]/(sMod_[i][celli] + VSMALL);
        }
    }
    return fractionsSum/mixture;
}


// ************************************************************************* //
