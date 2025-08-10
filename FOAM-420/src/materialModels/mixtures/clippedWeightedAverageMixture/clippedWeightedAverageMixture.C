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

#include "clippedWeightedAverageMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(clippedWeightedAverageMixture, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        clippedWeightedAverageMixture,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::volumeMassFractions&
Foam::clippedWeightedAverageMixture::lookupOrConstructBase()
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


Foam::clippedWeightedAverageMixture::clippedWeightedAverageMixture
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


Foam::autoPtr<Foam::clippedWeightedAverageMixture>
Foam::clippedWeightedAverageMixture::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<clippedWeightedAverageMixture>
    (
        new clippedWeightedAverageMixture
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

void Foam::clippedWeightedAverageMixture::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>*
Foam::clippedWeightedAverageMixture::castScalarModel
(
    const word& modelName
)
{
    if (clippedWeightedAverageMixtureModel::typeName == modelName)
    {
        return dynamic_cast<clippedWeightedAverageMixtureModel*>(this);
    }
    return nullptr;
}


Foam::baseModels<Foam::vector>*
Foam::clippedWeightedAverageMixture::castVectorModel
(
    const word& modelName
)
{
    if (clippedWeightedAverageMixtureModel::typeName == modelName)
    {
        return dynamic_cast<vectorWeightedAverageMixtureModel*>(this);
    }
    return nullptr;
}


Foam::baseModels<Foam::tensor>*
Foam::clippedWeightedAverageMixture::castTensorModel
(
    const word& modelName
)
{
    if (clippedWeightedAverageMixtureModel::typeName == modelName)
    {
        return dynamic_cast<tensorWeightedAverageMixtureModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::clippedWeightedAverageMixture::clippedWeightedAverageMixtureGeometric()
const
{
    return clippedWeightedAverageMixtureTypeGeometric<scalar>(sMod_);
}


Foam::tmp<Foam::scalarField>
Foam::clippedWeightedAverageMixture::clippedWeightedAverageMixtureInternal()
const
{
    return clippedWeightedAverageMixtureTypeInternal<scalar>(sMod_);
}


Foam::tmp<Foam::scalarField>
Foam::clippedWeightedAverageMixture::clippedWeightedAverageMixturePatch
(
    label const patchi
) const
{
    return clippedWeightedAverageMixtureTypePatch<scalar>(sMod_, patchi);
}


Foam::scalar
Foam::clippedWeightedAverageMixture::clippedWeightedAverageMixtureCell
(
    label const celli
) const
{
    return clippedWeightedAverageMixtureTypeCell<scalar>(sMod_, celli);
}


Foam::tmp<Foam::volVectorField>
Foam::clippedWeightedAverageMixture::
vectorClippedWeightedAverageMixtureGeometric()
const
{
    return clippedWeightedAverageMixtureTypeGeometric<vector>(vMod_);
}


Foam::tmp<Foam::vectorField>
Foam::clippedWeightedAverageMixture::
vectorClippedWeightedAverageMixtureInternal()
const
{
    return clippedWeightedAverageMixtureTypeInternal<vector>(vMod_);
}


Foam::tmp<Foam::vectorField>
Foam::clippedWeightedAverageMixture::vectorClippedWeightedAverageMixturePatch
(
    label const patchi
) const
{
    return clippedWeightedAverageMixtureTypePatch<vector>(vMod_, patchi);
}


Foam::vector
Foam::clippedWeightedAverageMixture::vectorClippedWeightedAverageMixtureCell
(
    label const celli
) const
{
    return clippedWeightedAverageMixtureTypeCell<vector>(vMod_, celli);
}


Foam::tmp<Foam::volTensorField>
Foam::clippedWeightedAverageMixture::
tensorClippedWeightedAverageMixtureGeometric()
const
{
    return clippedWeightedAverageMixtureTypeGeometric<tensor>(tMod_);
}


Foam::tmp<Foam::tensorField>
Foam::clippedWeightedAverageMixture::
tensorClippedWeightedAverageMixtureInternal()
const
{
    return clippedWeightedAverageMixtureTypeInternal<tensor>(tMod_);
}


Foam::tmp<Foam::tensorField>
Foam::clippedWeightedAverageMixture::tensorClippedWeightedAverageMixturePatch
(
    label const patchi
) const
{
    return clippedWeightedAverageMixtureTypePatch<tensor>(tMod_, patchi);
}


Foam::tensor
Foam::clippedWeightedAverageMixture::tensorClippedWeightedAverageMixtureCell
(
    label const celli
) const
{
    return clippedWeightedAverageMixtureTypeCell<tensor>(tMod_, celli);
}


// ************************************************************************* //
