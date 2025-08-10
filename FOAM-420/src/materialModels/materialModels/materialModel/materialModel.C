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

#include "materialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materialModel, 0);
    defineRunTimeSelectionTable(materialModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModel::materialModel
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    regIOobject
    (
        IOobject
        (
            name,
            obr.time().timeName(),
            obr.subRegistry("materialModels", true, false)
        )
    ),
    matObr_(obr.subRegistry("materialModels")),
    materialTables_(matObr_.lookupObjectRef<materialTables>("materialTables")),
    phaseName_(phaseName),
    specieName_(specieName),
    name_(name),
    mesh_(meshFromRegistry(obr)),
    obr_(obr),
    dict_(&dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModel::~materialModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::materialModel::meshFromRegistry
(
    const objectRegistry& obr
) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}


Foam::word Foam::materialModel::phasePropertyName
(
    const word fieldName,
    const word& phaseName
)
{
    // Function will need to account for different field names in
    // different solvers

    // TODO Test sequential inversion for each phase T
    if (fieldName == "p" || fieldName == "U" || fieldName == "T")
    {
        return fieldName;
    }
    return IOobject::groupName(fieldName, phaseName);
}


void Foam::materialModel::updateDictPtr()
{
    const bool isPhaseMixture = (phaseName_ == word::null) ? false : true;
    const bool isSpecieMixture = (specieName_ == word::null) ? false : true;

    const dictionary& matDict = materialsDict();

    const auto i = name_.rfind('.');

    word coeffDictName;

    if (i == std::string::npos || i == 0)
    {
        coeffDictName = name_;
    }
    else
    {
        coeffDictName = name_.substr(0, i);
    }
    coeffDictName += "Coeffs";

    // Model dictionary
    dict_ =
        isPhaseMixture ?
        (
            isSpecieMixture
          ? &matDict.subDict(phaseName_).subDict(specieName_).optionalSubDict(coeffDictName)
          : &matDict.subDict(phaseName_).optionalSubDict(coeffDictName)
        )
      : (
            isSpecieMixture
          ? &matDict.subDict(specieName_).optionalSubDict(coeffDictName)
          : &matDict.optionalSubDict(coeffDictName)
        );

    if (dict_ == nullptr)
    {
        FatalErrorInFunction
            << "Sub-dictionary not found."
            << exit(FatalError);
    }
}


void Foam::materialModel::fill
(
    const word& checkModel,
    const word& modelName,
    const wordList& modelNames,
    const labelList& inds,
    const boolList& modelsRequired,
    const label depNumber
)
{
    if (checkModel == modelName)
    {
        boolList modelsAdeed(modelNames.size(), false);
        const word tName(materialTables_.tableName(phaseName_, specieName_));

        depSubList oneDep;
        label s = 0;
        if (materialTables_.sTable().found(tName))
        {
            const matScalarTable& sModels =
                materialTables_.sTable(phaseName_, specieName_);
            if (sModels.found(modelName))
            {
                forAll(modelNames, i)
                {
                    const word& name = modelNames[i];
                    if (sModels.found(name))
                    {
                        oneDep.model = sModels[modelName];
                        oneDep.dependencies.setSize(++s);
                        modelsAdeed[i] = true;
                        sMod_.set(inds[i], sModels[name]);
                        oneDep.dependencies.set(s-1, sModels[name]);
                    }
                }
            }
        }

        if (materialTables_.vTable().found(tName))
        {
            const matVectorTable& vModels =
                materialTables_.vTable(phaseName_, specieName_);
            if (vModels.found(modelName))
            {
                forAll(modelNames, i)
                {
                    const word& name = modelNames[i];
                    if (vModels.found(name))
                    {
                        oneDep.model = vModels[modelName];
                        oneDep.dependencies.setSize(++s);
                        modelsAdeed[i] = true;
                        vMod_.set(inds[i], vModels[name]);
                        oneDep.dependencies.set(i, vModels[name]);
                    }
                }
            }
        }

        if (materialTables_.tTable().found(tName))
        {
            const matTensorTable& tModels =
                materialTables_.tTable(phaseName_, specieName_);
            if (tModels.found(modelName))
            {
                forAll(modelNames, i)
                {
                    const word& name = modelNames[i];
                    if (tModels.found(name))
                    {
                        oneDep.model = tModels[modelName];
                        oneDep.dependencies.setSize(++s);
                        modelsAdeed[i] = true;
                        tMod_.set(inds[i], tModels[name]);
                        oneDep.dependencies.set(i, tModels[name]);
                    }
                }
            }
        }
        if (oneDep.model != nullptr)
        {
            dep_[depNumber] = oneDep;
        }
        bool allUpdated = true;
        forAll(modelsAdeed, i)
        {
            allUpdated =
                (modelsRequired[i] == true && modelsAdeed[i] == false)
              ? false
              : allUpdated;
        }
        if (!allUpdated)
        {
            FatalErrorInFunction
                << "At least one of the required models for \""
                << modelName << "\""
                << " from the list: " << nl << modelNames << " not found."
                << exit(FatalError);
        }
    }
}


Foam::scalar Foam::materialModel::value
(
    const word& modelName,
    const label modelIndex,
    const label celli,
    const bool update,
    const bool readDict
)
{
    if (readDict)
    {
        sMod_[modelIndex].read();
    }
    if (update)
    {
        sMod_[modelIndex].updateTable(modelName);
    }
    return sMod_[modelIndex][celli];
}


const Foam::dictionary& Foam::materialModel::materialsDict() const
{
    return materialTables_.dict();
}


// ************************************************************************* //
