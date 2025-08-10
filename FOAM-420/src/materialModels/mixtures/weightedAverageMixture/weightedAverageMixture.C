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

#include "weightedAverageMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(weightedAverageMixture, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        weightedAverageMixture,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::volumeMassFractions&
Foam::weightedAverageMixture::lookupOrConstructBase()
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


Foam::weightedAverageMixture::weightedAverageMixture
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


Foam::autoPtr<Foam::weightedAverageMixture> Foam::weightedAverageMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<weightedAverageMixture>
    (
        new weightedAverageMixture
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

void Foam::weightedAverageMixture::updateTable(const word& modelName)
{
    const wordList names =
        dict_->found("species")
      ? dict_->lookup<wordList>("species")
      : dict_->lookup<wordList>("phases");
    const label nModels = names.size();
    const word ttName = materialTables_.tableName(phaseName_, word::null);

    const HashTable<matScalarTable>& sModels = materialTables_.sTable();
    const HashTable<matVectorTable>& vModels = materialTables_.vTable();
    const HashTable<matTensorTable>& tModels = materialTables_.tTable();
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
    else if (materialTables_.foundModel<vector>(ttName, funcName))
    {
        vMod_.setSize(nModels);
        dep_.setSize(1);
        dep_[0].model = vModels[ttName][funcName];
        dep_[0].dependencies.setSize(nModels);
        forAll(names, I)
        {
            const word tName =
                dict_->found("species")
              ? materialTables_.tableName(phaseName_, names[I])
              : materialTables_.tableName(names[I], word::null);
            vMod_.set(I, vModels[tName][funcName]);
            dep_[0].dependencies.set(I, vModels[tName][funcName]);
        }
    }
    else if (materialTables_.foundModel<tensor>(ttName, funcName))
    {
        tMod_.setSize(nModels);
        dep_.setSize(1);
        dep_[0].model = tModels[ttName][funcName];
        dep_[0].dependencies.setSize(nModels);
        forAll(names, I)
        {
            const word tName =
                dict_->found("species")
              ? materialTables_.tableName(phaseName_, names[I])
              : materialTables_.tableName(names[I], word::null);
            tMod_.set(I, tModels[tName][funcName]);
            dep_[0].dependencies.set(I, vModels[tName][funcName]);
        }
    }
}


Foam::baseModels<Foam::scalar>* Foam::weightedAverageMixture::castScalarModel
(
    const word& modelName
)
{
    if (weightedAverageMixtureModel::typeName == modelName)
    {
        return dynamic_cast<weightedAverageMixtureModel*>(this);
    }
    return nullptr;
}


Foam::baseModels<Foam::vector>* Foam::weightedAverageMixture::castVectorModel
(
    const word& modelName
)
{
    if (weightedAverageMixtureModel::typeName == modelName)
    {
        return dynamic_cast<vectorWeightedAverageMixtureModel*>(this);
    }
    return nullptr;
}


Foam::baseModels<Foam::tensor>* Foam::weightedAverageMixture::castTensorModel
(
    const word& modelName
)
{
    if (weightedAverageMixtureModel::typeName == modelName)
    {
        return dynamic_cast<tensorWeightedAverageMixtureModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::weightedAverageMixture::weightedAverageMixtureGeometric() const
{
    return weightedAverageMixtureTypeGeometric<scalar>(sMod_);
}


Foam::tmp<Foam::scalarField>
Foam::weightedAverageMixture::weightedAverageMixtureInternal() const
{
    return weightedAverageMixtureTypeInternal<scalar>(sMod_);
}


Foam::tmp<Foam::scalarField>
Foam::weightedAverageMixture::weightedAverageMixturePatch
(
    label const patchi
) const
{
    return weightedAverageMixtureTypePatch<scalar>(sMod_, patchi);
}


Foam::scalar Foam::weightedAverageMixture::weightedAverageMixtureCell
(
    label const celli
) const
{
    return weightedAverageMixtureTypeCell<scalar>(sMod_, celli);
}


Foam::tmp<Foam::volVectorField>
Foam::weightedAverageMixture::vectorWeightedAverageMixtureGeometric() const
{
    return weightedAverageMixtureTypeGeometric<vector>(vMod_);
}


Foam::tmp<Foam::vectorField>
Foam::weightedAverageMixture::vectorWeightedAverageMixtureInternal() const
{
    return weightedAverageMixtureTypeInternal<vector>(vMod_);
}


Foam::tmp<Foam::vectorField>
Foam::weightedAverageMixture::vectorWeightedAverageMixturePatch
(
    label const patchi
) const
{
    return weightedAverageMixtureTypePatch<vector>(vMod_, patchi);
}


Foam::vector Foam::weightedAverageMixture::vectorWeightedAverageMixtureCell
(
    label const celli
) const
{
    return weightedAverageMixtureTypeCell<vector>(vMod_, celli);
}


Foam::tmp<Foam::volTensorField>
Foam::weightedAverageMixture::tensorWeightedAverageMixtureGeometric() const
{
    return weightedAverageMixtureTypeGeometric<tensor>(tMod_);
}


Foam::tmp<Foam::tensorField>
Foam::weightedAverageMixture::tensorWeightedAverageMixtureInternal() const
{
    return weightedAverageMixtureTypeInternal<tensor>(tMod_);
}


Foam::tmp<Foam::tensorField>
Foam::weightedAverageMixture::tensorWeightedAverageMixturePatch
(
    label const patchi
) const
{
    return weightedAverageMixtureTypePatch<tensor>(tMod_, patchi);
}


Foam::tensor Foam::weightedAverageMixture::tensorWeightedAverageMixtureCell
(
    label const celli
) const
{
    return weightedAverageMixtureTypeCell<tensor>(tMod_, celli);
}


// ************************************************************************* //
