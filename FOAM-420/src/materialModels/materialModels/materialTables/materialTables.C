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

#include "materialTables.H"
#include "materialModels/materialModel/materialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(materialTables, 0);
}

const Foam::word Foam::materialTables::speciesListName("species");
const Foam::word Foam::materialTables::phasesListName("phases");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::materialTables::collectSolverFlagsInfo
(
    const word& phaseName,
    const word& specieName,
    const word& funcName,
    const materialModel* modelPtr
)
{
    if
    (
        funcName == kappaModel::typeName
     || funcName == CpvModel::typeName
     || funcName == rhoModel::typeName
    )
    {
        // Name not so important at the moment
        // (it just have to be different for each item)
        const word modName
        (
            IOobject::groupName(funcName, tableName(phaseName, specieName))
        );
        solverFlagsInfo_.set
        (
            modName,
            solverFlags
            ({
                funcName,
                phaseName,
                specieName,
                modelPtr->incompressible(),
                modelPtr->isochoric(),
                modelPtr->isCpvConst(),
                modelPtr->isotropic(),
                modelPtr->mixture()
            })
        );
    }
}


void Foam::materialTables::setSolverFlags()
{
    const wordList models(solverFlagsInfo_.toc());
    forAll(models, modelI)
    {
        const word& name = models[modelI];
        const solverFlags& mInfo = solverFlagsInfo_[name];
        if (isMultiSpecies_.lookup(mInfo.phaseName, false))
        {
            if (mInfo.specieName == word::null && mInfo.mixture)
            {
                if (mInfo.funcName == rhoModel::typeName)
                {
                    bool incompressibleMixture = true;
                    forAll(models, i)
                    {
                        const solverFlags& mixtureInfo =
                            solverFlagsInfo_[models[i]];
                        if
                        (
                            mixtureInfo.specieName != word::null
                         && mInfo.phaseName == mixtureInfo.phaseName
                        )
                        {
                            if (mixtureInfo.funcName == rhoModel::typeName)
                            {
                                isochoric_.set(mixtureInfo.phaseName, false);
                                if (!mixtureInfo.incompressible)
                                {
                                    incompressibleMixture = false;
                                }
                            }
                        }
                    }
                    incompressible_.set
                    (
                        mInfo.phaseName, incompressibleMixture
                    );
                }
                else if (mInfo.funcName == CpvModel::typeName)
                {
                    forAll(models, i)
                    {
                        const solverFlags& mixtureInfo =
                            solverFlagsInfo_[models[i]];
                        if
                        (
                            mixtureInfo.specieName != word::null
                         && mInfo.phaseName == mixtureInfo.phaseName
                        )
                        {
                            if (mixtureInfo.funcName == CpvModel::typeName)
                            {
                                // As soon as we have a mixure, we can assume
                                // Cpv is not uniform
                                constCpv_.set(mixtureInfo.phaseName, false);
                            }
                        }
                    }
                }
                else if (mInfo.funcName == kappaModel::typeName)
                {
                    bool isotropicMixture = true;
                    forAll(models, i)
                    {
                        const solverFlags& mixtureInfo =
                            solverFlagsInfo_[models[i]];
                        if
                        (
                            mixtureInfo.specieName != word::null
                         && mInfo.phaseName == mixtureInfo.phaseName
                        )
                        {
                            if
                            (
                                mixtureInfo.funcName == kappaModel::typeName
                            )
                            {
                                if (!mixtureInfo.isotropic)
                                {
                                    isotropicMixture = false;
                                }
                            }
                        }
                    }
                    isotropic_.set(mInfo.phaseName, isotropicMixture);
                }
            }
        }
        else
        {
            if (mInfo.funcName == rhoModel::typeName)
            {
                incompressible_.set(mInfo.phaseName, mInfo.incompressible);
                isochoric_.set(mInfo.phaseName, mInfo.isochoric);
            }
            else if (mInfo.funcName == CpvModel::typeName)
            {
                constCpv_.set(mInfo.phaseName, mInfo.constCpv);
            }
            else if (mInfo.funcName == kappaModel::typeName)
            {
                isotropic_.set(mInfo.phaseName, mInfo.isotropic);
            }
        }
    }
}


Foam::wordList Foam::materialTables::listModels
(
    const word& phaseName,
    const wordList& groupModelNames
)
{
    wordList speciesPhasesList;
    if (phaseName == word::null)
    {
        speciesPhasesList =
            dict_.found(speciesListName)
          ? dict_.lookup<wordList>(speciesListName)
          : dict_.lookup<wordList>(phasesListName);
    }
    else
    {
        speciesPhasesList =
            dict_.subDict(phaseName).lookup<wordList>(speciesListName);
    }

    const word tName0 = tableName(phaseName, speciesPhasesList.first());
    wordList initModelList = groupModelNames;

    if (groupModelNames == wordList::null())
    {
        if (scalarModels_.found(tName0))
        {
            initModelList.append(scalarModels_[tName0].toc());
        }
        if (vectorModels_.found(tName0))
        {
            initModelList.append(vectorModels_[tName0].toc());
        }
        if (tensorModels_.found(tName0))
        {
            initModelList.append(tensorModels_[tName0].toc());
        }
    }

    wordList modelList;
    bool allSpeciesPhases = true;
    forAll(initModelList, modelI)
    {
        const word modelName = initModelList[modelI];
        forAll(speciesPhasesList, speciesPhasei)
        {
            const word tName =
                tableName(phaseName, speciesPhasesList[speciesPhasei]);

            if (!foundModel(tName, modelName))
            {
                allSpeciesPhases = false;
            }
        }
        if (allSpeciesPhases && !modelList.found(modelName))
        {
            modelList.append(modelName);
        }
        else
        {
            allSpeciesPhases = true;
        }
    }

    return modelList;
}


void Foam::materialTables::groupNotFoundError
(
    const word& modelName,
    const word& noDictDepGroup
)
{
    if (noDictDepGroup != word::null)
    {
        FatalErrorInFunction
            << "If \"" << noDictDepGroup
            << "\" models aren't specified "
            << "in the \"materialProperties\" dictionary, the \""
            << modelName << "\" models have to be defined."
            << nl << exit(FatalError);
    }
    FatalErrorInFunction
        << "Settings for model type: \"" << modelName
        << "\" wasn't found in the \"materialProperties\" dictionary."
        << nl << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialTables::materialTables
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            materialTables::type(),
            obr.time().timeName(),
            obr.subRegistry("materialModels", true, false)
        )
    ),
    enableLinkingModels_
    (
        dict.lookupOrDefault<Switch>("enableLinkingModels", true)
    ),
    obr_(obr),
    matObr_(obr.subRegistry("materialModels")),
    dict_(dict),
    isMultiphase_(false),
    isMultiSpecies_()
{
    if
    (
        dict_.found(phasesListName)
     && dict_.lookup<wordList>(phasesListName).size() > 1
    )
    {
        isMultiphase_ = true;
        const wordList phases(dict_.lookup<wordList>(phasesListName));
        forAll(phases, phaseI)
        {
            const dictionary& phaseDict = dict_.subDict(phases[phaseI]);
            if
            (
                phaseDict.found(speciesListName)
             && phaseDict.lookup<wordList>(speciesListName).size() > 1
            )
            {
                isMultiSpecies_.set(phases[phaseI], true);
            }
            else
            {
                isMultiSpecies_.set(phases[phaseI], false);
            }
        }
    }
    else if
    (
        dict_.found(speciesListName)
     && dict_.lookup<wordList>(speciesListName).size() > 1
    )
    {
        isMultiSpecies_.set(word::null, true);
    }
    else
    {
        isMultiSpecies_.set(word::null, false);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::materialTables::tableName
(
    const word& phaseName,
    const word& specieName
) const
{
    const word phaseSpecieName =
        (specieName != word::null && phaseName == word::null)
      ? specieName
      : groupName(phaseName, specieName);

    if (phaseSpecieName == word::null)
    {
        return "mixture";
    }
    return phaseSpecieName;
}


const Foam::HashTable<Foam::matScalarTable>&
Foam::materialTables::sTable() const
{
    return scalarModels_;
}


const Foam::HashTable<Foam::matVectorTable>&
Foam::materialTables::vTable() const
{
    return vectorModels_;
}


const Foam::HashTable<Foam::matTensorTable>&
Foam::materialTables::tTable() const
{
    return tensorModels_;
}


Foam::HashTable<Foam::matScalarTable>&
Foam::materialTables::sTable()
{
    return scalarModels_;
}


Foam::HashTable<Foam::matVectorTable>&
Foam::materialTables::vTable()
{
    return vectorModels_;
}


Foam::HashTable<Foam::matTensorTable>&
Foam::materialTables::tTable()
{
    return tensorModels_;
}


Foam::matScalarTable& Foam::materialTables::sTable
(
    const word& phaseName,
    const word& specieName
)
{
    const word tName(tableName(phaseName, specieName));
    return scalarModels_[tName];
}


Foam::matVectorTable& Foam::materialTables::vTable
(
    const word& phaseName,
    const word& specieName
)
{
    const word tName(tableName(phaseName, specieName));
    return vectorModels_[tName];
}


Foam::matTensorTable& Foam::materialTables::tTable
(
    const word& phaseName,
    const word& specieName
)
{
    const word tName(tableName(phaseName, specieName));
    return tensorModels_[tName];
}


bool Foam::materialTables::addScalarModel
(
    const word& funcName,
    const word& phaseName,
    const word& specieName,
    materialModel* modelPtr,
    const word& castType
)
{
    // Try to add mixture types only for correct type
    bool addModel =
        modelPtr->mixture()
      ? foundModel<scalar>(firstSpeciesPhase(phaseName), funcName)
      : true;

    baseModels<scalar>* scalarPtr = modelPtr->castScalarModel(castType);
    if (scalarPtr != nullptr && addModel)
    {
        const word tName = tableName(phaseName, specieName);
        if (!scalarModels_.found(tName))
        {
            scalarModels_.set(tName, matScalarTable());
        }
        scalarModels_[tName].set(funcName, scalarPtr);
        return true;
    }
    return false;
}


bool Foam::materialTables::addVectorModel
(
    const word& funcName,
    const word& phaseName,
    const word& specieName,
    materialModel* modelPtr,
    const word& castType
)
{
    // Try to add mixture types only for correct type
    bool addModel =
        modelPtr->mixture()
      ? foundModel<vector>(firstSpeciesPhase(phaseName), funcName)
      : true;

    baseModels<vector>* vectorPtr = modelPtr->castVectorModel(castType);
    if (vectorPtr != nullptr && addModel)
    {
        const word tName = tableName(phaseName, specieName);
        if (!vectorModels_.found(tName))
        {
            vectorModels_.set(tName, matVectorTable());
        }
        vectorModels_[tName].set(funcName, vectorPtr);
        return true;
    }
    return false;
}


bool Foam::materialTables::addTensorModel
(
    const word& funcName,
    const word& phaseName,
    const word& specieName,
    materialModel* modelPtr,
    const word& castType
)
{
    // Try to add mixture types only for correct type
    bool addModel =
        modelPtr->mixture()
      ? foundModel<tensor>(firstSpeciesPhase(phaseName), funcName)
      : true;

    baseModels<tensor>* tensorPtr = modelPtr->castTensorModel(castType);
    if (tensorPtr != nullptr && addModel)
    {
        const word tName = tableName(phaseName, specieName);
        if (!tensorModels_.found(tName))
        {
            tensorModels_.set(tName, matTensorTable());
        }
        tensorModels_[tName].set(funcName, tensorPtr);
        return true;
    }
    return false;
}


Foam::word Foam::materialTables::firstSpeciesPhase
(
    const word& phaseName
)
{
    if (phaseName == word::null)
    {
        return
            dict_.found(speciesListName)
          ? dict_.lookup<wordList>(speciesListName).first()
          : dict_.lookup<wordList>(phasesListName).first();
    }

    const word speciesName =
        dict_.subDict
        (
            phaseName
        ).lookup<wordList>(speciesListName).first();

    return tableName(phaseName, speciesName);
}


const Foam::dictionary& Foam::materialTables::processDict
(
    const word& phaseName,
    const word& specieName,
    const word& modelName,
    const word& funcName,
    word& modelType,
    word& obrModelName,
    const word& noDictModelType
)
{
    const bool isPhaseMixture = (phaseName == word::null) ? false : true;
    const bool isSpecieMixture = (specieName == word::null) ? false : true;

    // Model dictionary
    const dictionary& modelDict =
            isPhaseMixture ?
            (
                isSpecieMixture
              ? dict_.subDict(phaseName).subDict(specieName)
              : dict_.subDict(phaseName)
            )
           : (isSpecieMixture ? dict_.subDict(specieName) : dict_);

    const word specModelName(funcName + "Model");

    if (modelDict.found(specModelName)) // Single model
    {
        modelType = modelDict.lookup<word>(specModelName);
        obrModelName =
            groupName(groupName(specModelName, phaseName), specieName);
        return modelDict.optionalSubDict(specModelName + "Coeffs");
    }
    else if (isMixture(modelName)) // Mixture
    {
        modelType = modelName;
        obrModelName =
            groupName(groupName(specModelName, phaseName), specieName);
        return modelDict;
    }
    else if (modelDict.found(modelName)) //  Group NO mixture
    {
        modelType = modelDict.lookup<word>(modelName);
        obrModelName = groupName(groupName(modelName, phaseName), specieName);
        return modelDict.optionalSubDict(modelName + "Coeffs");
    }

    else if (noDictModelType != word::null) // models/mixtures
    {
        modelType = noDictModelType;
        obrModelName = groupName(groupName(modelName, phaseName), specieName);
        return modelDict.optionalSubDict(noDictModelType + "Coeffs");
    }
    else
    {
        FatalErrorInFunction
            << "If type \"" << modelName << "\" isn't \"mixture\" "
            << "it needs to be defined in " << modelDict.name() << " dict."
            << exit(FatalError);
    }
    return modelDict;
}


word Foam::materialTables::standardMixture
(
    const word& funcName,// Can we check for cast type?
    const word& phaseName,
    const word& defaultType
)
{
    const dictionary& modelDict =
        (phaseName == word::null)
      ? dict_
      : dict_.subDict(phaseName);

    if (defaultType != word::null)
    {
        return defaultType;
    }

    const bool isSupportedType = modelDict.found("mixture");
    if (isSupportedType)
    {
        word mixtureType = modelDict.lookup<word>("mixture");

        // Check if the mixture is one of the 3 supported types
        if
        (
            mixtureType != "standardMixture"
         && mixtureType != "volumeMixture"
         && mixtureType != "massMixture"
        )
        {
            FatalErrorInFunction
                << "The \"mixture\" \"" << mixtureType
                << "\" is not supported. " << nl << modelDict.name()
                << nl << exit(FatalError);
        }

        // Translate standardMixture to volume/mass mixing
        if (mixtureType == "standardMixture")
        {
            if (phaseName == word::null && dict_.found(phasesListName))
            {
                mixtureType = "volumeMixture";
            }
            else
            {
                mixtureType = "massMixture";
            }
        }

        if (mixtureType == "massMixture")
        {
            if
            (
                   funcName == rhoModel::typeName
                || funcName == psiModel::typeName
                || funcName == WModel::typeName
                || funcName == "buoyancy"
            )
            {
                return harmonicMixtureModel::typeName;
            }
            else
            {
                return weightedAverageMixtureModel::typeName;
            }
        }
        else if (mixtureType == "volumeMixture")
        {
            return clippedWeightedAverageMixtureModel::typeName;
        }
    }
    else
    {
        FatalErrorInFunction
            << "The \"mixture\" is not defined in: " << modelDict.name() << nl
            << nl << exit(FatalError);
    }
    return word::null;
}


void Foam::materialTables::add
(
    const word& modelName,
    const word& funcName,
    const word& phaseName,
    const word& specieName,
    const word& noDictModelType,
    const word& modelCastType
)
{
    // Model type loaded from dictionary. Looked up from construction table.
    // Example: "perfectGas"
    word modelType;

    // The model will be stored in object registry
    // under this name
    word obrModelName;
    const dictionary& modelDict =
        processDict
        (
            phaseName,
            specieName,
            modelName, // rhoModel, equationOfState and !!harmonicMixture etc!!
            funcName, // rho/Cp/Cv etc.
            modelType, // construct. table models rhoConst/harmonicMixture etc.
            obrModelName, // Name in object registry
            noDictModelType //rho/Cp/Cv
        );

    if ((modelType.find("Model") != std::string::npos) && enableLinkingModels_)
    {
        const word modName
        (
            IOobject::groupName(funcName, tableName(phaseName, specieName))
        );
        const word linkToModel(modelType.replace("Model", ""));
        linkingModelsInfo_.set
        (
            modName,
            linkingInfo
            ({
                funcName,
                phaseName,
                specieName,
                linkToModel,
                modelCastType
            })
        );
        return;
    }

    if (!matObr_.foundObject<materialModel>(obrModelName))
    {
        materialModel::New
        (
            obr_,
            modelDict,
            modelType,
            phaseName,
            specieName,
            obrModelName
        ).ptr()->store();
    }
    materialModel* modelPtr =
        matObr_.lookupObjectRefPtr<materialModel>(obrModelName);

    // Change cast type in case it it explicitly defined as different type
    // or if it is mixture type model
    const word castType
    (
        (modelCastType != word::null)
      ? modelCastType // Explicitly defined castType
      : (
            modelPtr->mixture()
          ? modelPtr->type() // Cast type for mixtures
          : funcName         // Default cast type
        )
    );

    const bool sAdded =
        addScalarModel(funcName, phaseName, specieName, modelPtr, castType);

    const bool vAdded =
        addVectorModel(funcName, phaseName, specieName, modelPtr, castType);

    const bool tAdded =
        addTensorModel(funcName, phaseName, specieName, modelPtr, castType);

    collectSolverFlagsInfo(phaseName, specieName, funcName, modelPtr);

    if (debug && !sAdded && !vAdded && !tAdded)
    {
        Warning<< "";
        if (phaseName != word::null)
        {
            Info<< "phase \"" << phaseName << "\" ";
        }
        if (specieName != word::null)
        {
            Info<< "specie \"" << specieName << "\" ";
        }
        Info<< "type \"" << modelType << "\" doesn't contain \""
            << funcName  << "\"." << endl;
    }
}


void Foam::materialTables::addGroup
(
    const wordList& modelNames,
    const wordList& funcNames,
    const word& phaseName,
    const word& specieName,
    const word& noDictModelType,
    bool isModelMixture,
    const word& mixModelName,
    const wordList& modelCastTypes
)
{
    // List models already requires submodels to be available
    wordList corrFuncNames =
        isModelMixture ? listModels(phaseName, funcNames) : funcNames;

    forAll(corrFuncNames, funcI)
    {
        // Support for two ways of defining the groups:
        // 1) By group name (one entry). Example: equationOfState
        // 2) By listing all the group models. Example: rhoModel, psiModel
        word modelName =
            (modelNames.size() == 1) ? modelNames.first() : modelNames[funcI];
        if (isModelMixture)
        {
            modelName =
                standardMixture(corrFuncNames[funcI], phaseName, mixModelName);
        }
        word modelCastType = word::null;
        if
        (
            modelCastTypes.size()
         && (modelCastTypes != wordList::null())
         && !isModelMixture
        )
        {
            modelCastType = modelCastTypes[funcI];
        }
        add
        (
            modelName,
            corrFuncNames[funcI],
            phaseName,
            specieName,
            isModelMixture ? word::null : noDictModelType,
            modelCastType
        );
    }
}


bool Foam::materialTables::addSpecies
(
    const word& modelName,
    const wordList& funcNames,
    const word& phaseName,
    const word& noDictModelType,
    const word& noDictDepGroup,
    const wordList& modelCastTypes
)
{
    bool groupExhist = false;
    const dictionary& phaseDict =
        phaseName != word::null ? dict_.subDict(phaseName) : dict_;
    if (phaseDict.found(speciesListName))
    {
        const wordList species =
            phaseDict.lookup<wordList>(speciesListName);
        forAll(species, specieI)
        {
            const word& specieName = species[specieI];
            const dictionary& speciesDictI = phaseDict.subDict(specieName);
            const bool isDefaultDefined
            (
                noDictModelType != word::null
             && speciesDictI.found(noDictDepGroup)
            );

            if (speciesDictI.found(modelName) || isDefaultDefined)
            {
                addGroup
                (
                    {modelName},
                    funcNames,
                    phaseName,
                    specieName,
                    noDictModelType,
                    false,
                    word::null,
                    modelCastTypes
                );
                groupExhist = true;
            }
        }
    }
    return groupExhist;
}


bool Foam::materialTables::addPhases
(
    const word& modelName,
    const wordList& funcNames,
    const word& noDictModelType,
    const word& noDictDepGroup,
    const wordList& modelCastTypes
)
{
    bool groupExhist = false;
    if (dict_.found(phasesListName))
    {
        const wordList phases =
            dict_.lookup<wordList>(phasesListName);
        forAll(phases, phasei)
        {
            const word& phaseName = phases[phasei];
            const dictionary& phaseiDict = dict_.subDict(phaseName);
            const bool isDefaultDefined
            (
                noDictModelType != word::null
             && phaseiDict.found(noDictDepGroup)
            );

            if (phaseiDict.found(modelName) || isDefaultDefined)
            {
                addGroup
                (
                    {modelName},
                    funcNames,
                    phaseName,
                    word::null,
                    noDictModelType,
                    false,
                    word::null,
                    modelCastTypes
                );
                groupExhist = true;
            }
        }
    }
    return groupExhist;
}


bool Foam::materialTables::addSpeciesMixtures
(
    const word& modelName,
    const wordList& funcNames,
    const word& phaseName,
    const word& noDictModelType,
    const word& noDictDepGroup,
    const wordList& modelCastTypes,
    bool groupExhist
)
{
    const dictionary& phaseDict =
        (phaseName != word::null) ? dict_.subDict(phaseName) : dict_;
    const bool isDefaultDefined
    (
        noDictModelType != word::null
     && phaseDict.found(noDictDepGroup)
    );
    const bool isGroupInDict(phaseDict.found(modelName));
    if (isGroupInDict || isDefaultDefined || groupExhist)
    {
        const word groupType =
            isGroupInDict ? phaseDict.lookup<word>(modelName) : word::null;

        // Mixture from the dict
        const word mixFromDict(isMixture(groupType) ? groupType : word::null);

        // Is the group mixture (even undefined one)
        const bool isModelMixture =
            isMixture(groupType)
         || (groupExhist && !isGroupInDict && !isDefaultDefined)
         || (isDefaultDefined && isMixture(noDictModelType));
        addGroup
        (
            {modelName},
            funcNames,
            phaseName,
            word::null,
            noDictModelType,
            isModelMixture,
            mixFromDict,
            modelCastTypes
        );

        return true;
    }

    return groupExhist;
}


bool Foam::materialTables::addPhaseMixtures
(
    const word& modelName,
    const wordList& funcNames,
    const word& noDictModelType,
    const word& noDictDepGroup,
    const wordList& modelCastTypes,
    bool groupExhist
)
{
    const bool isDefaultDefined
    (
        noDictModelType != word::null
     && dict_.found(noDictDepGroup)
    );
    const bool isGroupInDict(dict_.found(modelName));
    if (isGroupInDict || isDefaultDefined || groupExhist)
    {
        const word groupType =
            isGroupInDict ? dict_.lookup<word>(modelName) : word::null;

        // Mixture from the dict
        const word mixFromDict(isMixture(groupType) ? groupType : word::null);

        // Is the group mixture (even undefined one)
        const bool isModelMixture =
            isMixture(groupType)
         || (groupExhist && !isGroupInDict && !isDefaultDefined)
         || (isDefaultDefined && isMixture(noDictModelType));

        addGroup
        (
            {modelName},
            funcNames,
            word::null,
            word::null,
            noDictModelType,
            isModelMixture,
            mixFromDict,
            modelCastTypes
        );

        return true;
    }

    return groupExhist;
}


void Foam::materialTables::addSpeciesAndSpeciesMixtures
(
    const word& modelName,
    const wordList& modelNames,
    const word& phaseName,
    const word& noDictModelType,
    const word& noDictDepGroup,
    const wordList& modelCastTypes
)
{
    bool groupExhist =
        addSpecies
        (
            modelName,
            modelNames,
            phaseName,
            noDictModelType,
            noDictDepGroup,
            modelCastTypes
        );

    groupExhist =
        addSpeciesMixtures
        (
            modelName,
            modelNames,
            phaseName,
            noDictModelType,
            noDictDepGroup,
            modelCastTypes,
            groupExhist
        );

    // Convenience location to add phase mixtures
    // (After models for last phase)
    if
    (
        dict_.found(phasesListName)
     && phaseName == dict_.lookup<wordList>(phasesListName).last()
    )
    {
        groupExhist =
            addPhaseMixtures
            (
                modelName,
                modelNames,
                noDictModelType,
                noDictDepGroup,
                modelCastTypes,
                groupExhist
            );

        // Give error if group doesn't exhist and defaults are not provided
        // TODO: This has to do with adding phase mixtures.
        // Should be handled better.
        if (!groupExhist && modelName != "energy")
        {
            groupNotFoundError(modelName, noDictDepGroup);
        }

    }

}


bool Foam::materialTables::isModelInDict
(
    const word& modelName, const word& phaseName
) const
{
    const dictionary& phaseDict =
        phaseName == word::null ? dict_ : dict_.subDict(phaseName);
    if (phaseDict.found(modelName))
    {
        return true;
    }
    else if (phaseDict.found(speciesListName))
    {
        const wordList species =
            phaseDict.lookup<wordList>(speciesListName);
        forAll(species, specieI)
        {
            const word& specieName = species[specieI];
            if (phaseDict.subDict(specieName).found(modelName))
            {
                return true;
            }
        }
    }
    if (phaseDict.found(phasesListName))
    {
        const wordList phases = phaseDict.lookup<wordList>(phasesListName);
        forAll(phases, phaseI)
        {
            if (isModelInDict(modelName, phases[phaseI]))
            {
                return true;
            }
        }
    }

    return false;
}


void Foam::materialTables::addModelReferences()
{
    const wordList infoNames(linkingModelsInfo_.toc());
    if (infoNames.size())
    {
        Info<< endl;
        forAll(infoNames, nameI)
        {
            const linkingInfo& linkInfo = linkingModelsInfo_[infoNames[nameI]];
            const word& funcName = linkInfo.funcName;
            const word& phaseName = linkInfo.phaseName;
            const word& specieName = linkInfo.specieName;
            const word& linkToModel = linkInfo.linkToModel;
            const word tName = tableName(phaseName, specieName);

            Info<< "Linking \"" << infoNames[nameI]
                << "\" to: \"" << IOobject::groupName(linkToModel, tName)
                << "\" model." << endl;

            if (foundModel<scalar>(tName, linkToModel))
            {
                if (!foundModel<scalar>(tName, funcName))
                {
                    scalarModels_[tName].set
                    (
                        funcName,
                        scalarModels_[tName][linkToModel]
                    );
                }
            }
            else if (foundModel<vector>(tName, linkToModel))
            {
                if (!foundModel<vector>(tName, funcName))
                {
                    vectorModels_[tName].set
                    (
                        funcName,
                        vectorModels_[tName][linkToModel]
                    );
                }
            }
            else if (foundModel<tensor>(tName, linkToModel))
            {
                if (!foundModel<tensor>(tName, funcName))
                {
                    tensorModels_[tName].set
                    (
                        funcName,
                        tensorModels_[tName][linkToModel]
                    );
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Linking model for: \"" << funcName << "\" type: \""
                    << linkToModel << "\" wasn't found."
                    << nl << exit(FatalError);
            }
        }

        // Add phase mixtures
        forAll(infoNames, nameI)
        {
            const linkingInfo& linkInfo = linkingModelsInfo_[infoNames[nameI]];
            const word& phaseName = linkInfo.phaseName;
            if (phaseName != word::null)
            {
                const word& funcName = linkInfo.funcName;
                const word tName = tableName(phaseName, word::null);

                // Check if phase mixture exhists
                if (!foundModel(tName, funcName))
                {
                    // Check if the model is defined in each species
                    if (listModels(phaseName, wordList({funcName})).size())
                    {
                        // Third we have to check which mixing
                        // for the model should be used
                        const word mixtureType
                        (
                            standardMixture(funcName, phaseName, word::null)
                        );
                        add
                        (
                            mixtureType,
                            funcName,
                            phaseName,
                            word::null,
                            word::null,
                            word::null
                        );
                    }
                }
            }
        }

        // Add global mixture
        forAll(infoNames, nameI)
        {
            const linkingInfo& linkInfo = linkingModelsInfo_[infoNames[nameI]];
            const word& funcName = linkInfo.funcName;
            const word tName = tableName(word::null, word::null);

            // Check if global mixture exhists
            if (!foundModel(tName, funcName))
            {
                // Check if the model is defined in each species
                if (listModels(word::null, wordList({funcName})).size())
                {
                    // Third we have to check which mixing
                    // for the model should be used
                    const word mixtureType
                    (
                        standardMixture(funcName, word::null, word::null)
                    );
                    add
                    (
                        mixtureType,
                        funcName,
                        word::null,
                        word::null,
                        word::null,
                        word::null
                    );
                }
            }
        }
    }
}


void Foam::materialTables::linkModelsAndUpdateTables()
{
    setSolverFlags();
    if (enableLinkingModels_)
    {
        addModelReferences();
    }
    updateTable<matScalarTable>(scalarModels_);
    updateTable<matVectorTable>(vectorModels_);
    updateTable<matTensorTable>(tensorModels_);
}


void Foam::materialTables::checkDepenencies()
{
    const wordList mToCheck = matObr_.names<materialModel>();
    forAll(mToCheck, modelI)
    {
        const materialModel& model =
            matObr_.lookupObject<materialModel>(mToCheck[modelI]);
        const depList& d = model.dep();
        forAll(d, depI)
        {
            baseBaseModels* model1 = d[depI].model;
            const UPtrList<baseBaseModels>& deps1 = d[depI].dependencies;
            forAll(deps1, I)
            {
                const depList& dd = deps1[I].dep();
                forAll(dd, depJ)
                {
                    baseBaseModels* model2 = dd[depJ].model;
                    forAll(deps1, J)
                    {
                        if (model2 == deps1(J) && model2 != nullptr)
                        {
                            const UPtrList<baseBaseModels>& deps2 =
                                dd[depJ].dependencies;
                            forAll(deps2, K)
                            {
                                if (deps2(K) == model1 && model1 != nullptr)
                                {
                                    FatalErrorInFunction
                                        << "Model \"" << mToCheck[modelI]
                                        << "\""
                                        << " type: \"" << model1->type()
                                        << "\" function "
                                        << "\"" << model1->funcType() << "\" "
                                        << nl
                                        << "is circularly dependent on: "
                                        << "\"" << model2->type() << "\""
                                        << " function: \""
                                        << model2->funcType() << "\""
                                        << nl << exit(FatalError);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void Foam::materialTables::updateScalarField
(
    const word& fieldName,
    const volScalarField& field
)
{
    // Material models update still needed for temperature inverse
    wordList matNames = matObr_.names<materialModel>();
    forAll(matNames, matI)
    {
        materialModel& model =
            matObr_.lookupObjectRef<materialModel>(matNames[matI]);
        model.updateScalarField(fieldName, field);
    }

    wordList refNames = obr_.names<refScalarField>();
    forAll(refNames, refI)
    {
        refScalarField& model =
            obr_.lookupObjectRef<refScalarField>(refNames[refI]);
        model.updateScalarField(fieldName, field);
    }
}


void Foam::materialTables::updateVectorField
(
    const word& fieldName,
    const volVectorField& field
)
{
    wordList matNames = matObr_.names<materialModel>();
    forAll(matNames, matI)
    {
        materialModel& model =
            matObr_.lookupObjectRef<materialModel>(matNames[matI]);
        model.updateVectorField(fieldName, field);
    }

    wordList refNames = obr_.names<refVectorField>();
    forAll(refNames, refI)
    {
        refVectorField& model =
            obr_.lookupObjectRef<refVectorField>(refNames[refI]);
        model.updateVectorField(fieldName, field);
    }
}


void Foam::materialTables::updateTensorField
(
    const word& fieldName,
    const volTensorField& field
)
{
    wordList matNames = matObr_.names<materialModel>();
    forAll(matNames, matI)
    {
        materialModel& model =
            matObr_.lookupObjectRef<materialModel>(matNames[matI]);
        model.updateTensorField(fieldName, field);
    }

    wordList refNames = obr_.names<refTensorField>();
    forAll(refNames, refI)
    {
        refTensorField& model =
            obr_.lookupObjectRef<refTensorField>(refNames[refI]);
        model.updateTensorField(fieldName, field);
    }
}


void Foam::materialTables::updateSubDictsPtrs()
{
    wordList matNames = matObr_.names<materialModel>();
    forAll(matNames, matI)
    {
        materialModel& model =
            matObr_.lookupObjectRef<materialModel>(matNames[matI]);
        model.updateDictPtr();
    }
}


bool Foam::materialTables::read() const
{
    // Would it be better to do update through object registry?
    wordList topLevel = scalarModels_.toc();
    forAll(topLevel, tableI)
    {
        const word tName(topLevel[tableI]);
        const wordList modelsToUpdate = scalarModels_[tName].toc();
        forAll(modelsToUpdate, modelI)
        {
            scalarModels_[tName][modelsToUpdate[modelI]]->read();
        }
    }
    topLevel = vectorModels_.toc();
    forAll(topLevel, tableI)
    {
        const word tName(topLevel[tableI]);
        const wordList modelsToUpdate = vectorModels_[tName].toc();
        forAll(modelsToUpdate, modelI)
        {
            vectorModels_[tName][modelsToUpdate[modelI]]->read();
        }
    }
    topLevel = tensorModels_.toc();
    forAll(topLevel, tableI)
    {
        const word tName(topLevel[tableI]);
        const wordList modelsToUpdate = tensorModels_[tName].toc();
        forAll(modelsToUpdate, modelI)
        {
            tensorModels_[tName][modelsToUpdate[modelI]]->read();
        }
    }
    return true;
}


void Foam::materialTables::reportModelsInfo() const
{
    if (debug)
    {
        Info<< nl << nl << "Scalar models" << nl << "-------------" << endl;
        reportModelsInfo<matScalarTable>(scalarModels_);

        Info<< nl << nl << "Vector models" << nl << "-------------" << endl;
        reportModelsInfo<matVectorTable>(vectorModels_);

        Info<< nl << nl << "Tensor models" << nl << "-------------" << endl;
        reportModelsInfo<matTensorTable>(tensorModels_);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

const Foam::baseModels<Foam::scalar>& Foam::materialTables::operator()
(
    const word& modelName,
    const word& phaseName,
    const word& specieName
) const
{
    const word tName = tableName(phaseName, specieName);
    return *scalarModels_[tName][modelName];
}


const Foam::baseModels<Foam::scalar>& Foam::materialTables::operator()
(
    const word& modelName,
    const word& phaseName
) const
{
    const word tName = tableName(phaseName, word::null);
    return *scalarModels_[tName][modelName];
}


const Foam::baseModels<Foam::scalar>& Foam::materialTables::operator()
(
    const word& modelName
) const
{
    return *scalarModels_["mixture"][modelName];
}


// ************************************************************************* //
