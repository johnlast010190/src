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
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "foamSolve.H"
#include "cfdTools/general/include/fvCFD.H"
#include "primitives/strings/fileName/fileName.H"
#include "regionProperties/regionProperties.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "monolithicSolve.H"
#include "algorithms/serialThreads/serialThreads.H"
#include "algorithms/serialThreads/serialThreadsTemplates.C"
#include "global/etcFiles/etcFiles.H"
#include "db/Time/timeSelector.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(foamSolve, 0);

} // End namespace Foam


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

bool scheduleItem::compareFunc::operator()
(
    const scheduleItem* ps1, const scheduleItem* ps2
)
{
    // If one of them is a corrector, base the comparisons on the
    // first item in the corrector
    const solveInfo& s1 = ps1->firstSolveItem();
    const solveInfo& s2 = ps2->firstSolveItem();

    // Differentiate by coupled group membership first, to keep together.
    // Use (region, field) key of first member as a unique
    // identifier
    if (s1.coupledGroup.size() > 1 || s2.coupledGroup.size() > 1)
    {
        if (s1.coupledGroup[0] != s2.coupledGroup[0])
        {
            return (s1.coupledGroup[0] < s2.coupledGroup[0]);
        }
    }

    // Next, differentiate by membership of region group
    if (s1.groupName != s2.groupName)
    {
        return (s1.groupName < s2.groupName);
    }

    // Tie-breakers - original region
    // ordering in regionProperties ...
    if (s1.regionOrder != s2.regionOrder)
    {
        return (s1.regionOrder < s2.regionOrder);
    }

    // ... then original ordering in fvOptions
    if (s1.solverOrder != s2.solverOrder)
    {
        return (s1.solverOrder < s2.solverOrder);
    }

    // ... then ordering of solves as reported by solver object
    if (s1.fieldOrder != s2.fieldOrder)
    {
        return (s1.fieldOrder < s2.fieldOrder);
    }

    return false;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void foamSolve::getRegions()
{
    // Read region groups from regionProperties file. It defaults to single
    // global mesh if not found
    regionProperties rp(time_());

    groupNames_ = rp.groupNames();
    groupConsolidated_.resize(groupNames_.size());


    const solutionInstance& instance = scheduler_->getCurrentInstance();
    regionNames_ = instance.regionNames();
    label k = 0;
    regionGroup_.resize(regionNames_.size());

    forAll(groupNames_, i)
    {
        const wordList& regionNamesGroup = rp[groupNames_[i]];
        forAll(regionNamesGroup, j)
        {
            if (regionNames_.found(regionNamesGroup[j]))
            {
                regionGroup_[k] = i;
                k++;
            }
        }
        groupConsolidated_[i] =
            rp.consolidatedGroupNames().found(groupNames_[i]);
    }

}

void Foam::foamSolve::createOptions(const dictionary& configOptionsDict)
{
    // Create fvOptions for each mesh
    fvOptions_.clear();
    fvOptions_.resize(solutionRegistries_.size());

    forAll(solutionRegistries_, i)
    {
        Info<< "Creating fvOption for mesh "
             << solutionRegistries_[i].mesh().name()
             << nl << "for region  "
             << solutionRegistries_[i].regionName()
             << endl;

        dictionary optionsDict =
            configOptionsDict.subOrEmptyDict("regionGroups");

        optionsDict =
            optionsDict.subOrEmptyDict(groupNames_[regionGroup_[i]]);
        if (groupNames_[regionGroup_[i]]!="multipointAdjoint")
        {
            optionsDict.merge(configOptionsDict.subOrEmptyDict("allRegions"));
        }

        fvOptions_.set
        (
            i,
            &fv::options::New
            (
                solutionRegistries_[i].mesh(),
                solutionRegistries_[i].registry(),
                optionsDict,
                false
            )
        );
    }
    Info<< endl;

    // Initialise
    forAll(fvOptions_, i)
    {
        fv::optionList& opts = fvOptions_[i];
        orderFvOptions(opts);
        fvOptions_[i].initialise();
    }

    // Pick out the solverObjects
    solverObjects_.clear();
    solverObjects_.resize(fvOptions_.size());
    forAll(fvOptions_, i)
    {
        solverObjects_.set
        (
            i,
            new solverObjectList(fvOptions_[i], false)
        );
    }
}

void Foam::foamSolve::orderFvOptions(fv::optionList& options)
{
    // First, reorder so that non-solver objects are initialised first
    labelList oldToNew(options.size());
    label newIdx = 0;
    forAll(options, j)
    {
        if (!isA<fv::solverOption>(options[j]))
        {
            oldToNew[j] = newIdx++;
        }
        else
        {
            refCast<fv::solverOption>(options[j]).setSolveController
            (
                *this
            );
        }
    }
    // Second, reorder so that primal objects are initialised first
    forAll(options, j)
    {
        if (isA<fv::solverOption>(options[j]))
        {
            bool isAdjoint = refCast<fv::solverOption>(options[j]).solverObject().adjoint();
            if (!isAdjoint)
            {
                oldToNew[j] = newIdx++;
            }
        }
    }
    // Last, reorder so that adjoint objects are initialised last
    forAll(options, j)
    {
        if (isA<fv::solverOption>(options[j]))
        {
            bool isAdjoint = refCast<fv::solverOption>(options[j]).solverObject().adjoint();
            if (isAdjoint)
            {
                oldToNew[j] = newIdx++;
            }
        }
    }
    options.reorder(oldToNew);
}

void Foam::foamSolve::gatherDependencyInfo()
{
    correctorNames_.clear();
    correctorMembers_.clear();
    allSolves_.clear();
    allCorrectors_.clear();

    // Store all regions' dependencies for later processing below after
    // initial processing on all regions
    List<SolveTable<List<Tuple2<solveID, bool>>>> allDependencies
    (
        solverObjects_.size()
    );

    forAll(solverObjects_, i)
    {
        solverObjectList& objects = solverObjects_[i];
        word regionName = solutionRegistries_[i].regionName();

        DynamicList<word> solveNames;
        HashTable<wordList> derivedFields;
        DynamicList<label> solverOrder;
        HashTable<wordList> corrMembers;

        objects.getSolveGraph
        (
            solveNames,
            derivedFields,
            allDependencies[i],
            solverOrder,
            corrMembers
        );

        // Consolidate any correctors with the same name across different solver
        // objects
        forAllConstIters(corrMembers, corr)
        {
            if (!correctorNames_.found(corr.key()))
            {
                correctorNames_.append(corr.key());
                correctorMembers_.append(solveList());
            }
            forAll(correctorNames_, k)
            {
                if (correctorNames_[k] == corr.key())
                {
                    // Merge lists
                    forAll(corr(), si)
                    {
                        solveID key(corr()[si], regionName);
                        if (!correctorMembers_[k].found(key))
                        {
                            correctorMembers_[k].append(key);
                        }
                    }
                    break;
                }
            }
        }

        fvMesh& mesh = solutionRegistries_[i].mesh();
        const objectRegistry& solReg = solutionRegistries_[i].registry();
        forAll(solveNames, j)
        {
            autoPtr<solveInfo> infoPtr(new solveInfo);
            label regioni = regionNames_[regionName];
            infoPtr->solveName = solveNames[j];
            infoPtr->solveNameAndFields.insert(solveNames[j]);
            infoPtr->solveNameAndFields.insert
            (
                derivedFields.lookup(solveNames[j], wordList())
            );
            infoPtr->regionName = regionName;
            infoPtr->groupName = groupNames_[regionGroup_[regioni]];
            infoPtr->regionOrder = regioni;
            infoPtr->solverOrder = solverOrder[j];
            infoPtr->fieldOrder = j;

            // Get cross-region matrix coupling from this field's BCs
            getCoupledFields<scalar>
            (
                mesh,
                solReg,
                solveNames[j],
                infoPtr->coupledFields
            );
            getCoupledFields<vector>
            (
                mesh,
                solReg,
                solveNames[j],
                infoPtr->coupledFields
            );
            getCoupledFields<sphericalTensor>
            (
                mesh,
                solReg,
                solveNames[j],
                infoPtr->coupledFields
            );
            getCoupledFields<symmTensor>
            (
                mesh,
                solReg,
                solveNames[j],
                infoPtr->coupledFields
            );
            getCoupledFields<tensor>
            (
                mesh,
                solReg,
                solveNames[j],
                infoPtr->coupledFields
            );

            // Record field to solve name mapping
            forAllConstIters(infoPtr->solveNameAndFields, f)
            {
                if
                (
                    !fieldToSolveName_.insert
                    (
                        solveID(f.key(), regionName),
                        solveID(infoPtr->solveName, regionName)
                    )
                )
                {
                    FatalErrorInFunction
                        << "Derived field " << f.key() << " is duplicated "
                        << "in region " << regionName << ". "
                        << "Clashing solve names: " << infoPtr->solveName
                        << " and "
                        << fieldToSolveName_
                           [
                               solveID(f.key(), regionName)
                           ].first()
                        << exit(FatalError);
                }
            }

            // Place in schedule
            solveID key(infoPtr->solveName, infoPtr->regionName);
            if (allSolves_.found(key))
            {
                FatalErrorInFunction
                    << "Solve " << solveNames[j] << " is duplicated "
                    << "in region " << regionName
                    << exit(FatalError);
            }
            allSolves_.set
            (
                key,
                infoPtr.ptr()
            );
        }
    }

    forAll(solverObjects_, i)
    {
        // Allocate dependencies
        forAllConstIters(allDependencies[i], depSolves)
        {
            // Convert to primary solve names
            solveID key = depSolves.key();
            key = fieldToSolveName_.lookup(key, key);
            forAll(depSolves(), si)
            {
                solveID depKey(depSolves()[si].first());
                bool compulsory = depSolves()[si].second();
                depKey = fieldToSolveName_.lookup(depKey, depKey);
                if (allSolves_.found(depKey))
                {
                    if (allSolves_.found(key))
                    {
                        allSolves_[key]->dependencies.append(depKey);
                    }
                    else if (compulsory)
                    {
                        // We ignore optional dependencies of solves that are
                        // not found.
                        FatalErrorInFunction
                            << "Solve " << key
                            << ", which was not found, has a compulsory "
                            << "dependency on solve " << depKey
                            << nl << exit(FatalError);
                    }
                }
                else if (compulsory)
                {
                    // We ignore optional dependencies on solves that are
                    // not found.
                    FatalErrorInFunction
                        << "Solve " << key << " has a compulsory "
                        << "dependency on solve " << depKey
                        << ", which was not found."
                        << nl << exit(FatalError);
                }
            }
        }
    }

    if (debug)
    {
        Info<< "\nGathered solve information:" << nl;
        forAllConstIters(allSolves_, iter)
        {
            Info<< "Solve: " << iter.key();
            Info<< " Dependencies: " << flatOutput(iter()->dependencies);
            Info<< " Coupled group: " << flatOutput(iter()->coupledFields);
            Info<< nl;
        }
        Info<< endl;

        Info<< "\nGathered corrector information:" << nl;
        forAll(correctorNames_, i)
        {
            Info<< "Corrector: " << correctorNames_[i];
            Info<< " Members: " << flatOutput(correctorMembers_[i]);
            Info<< nl;
        }
        Info<< endl;
    }
}


void Foam::foamSolve::checkSolveInfo()
{
    // At this point, allSolves_ only contains solve names.
    // Check that none of them collide with corrector names
    forAllConstIters(allSolves_, solve)
    {
        if (correctorNames_.found(solve.key().first()))
        {
            FatalErrorInFunction
                << "Corrector name '" << solve.key().first()
                << "' is the same as the name of a solve in region '"
                << solve.key().second()
                << "'. Solve names and corrector names may not correspond."
                << nl << exit(FatalError);
        }
    }
}


void Foam::foamSolve::substituteCorrectorNames
(
    const label correctori,
    DynamicList<label>& correctorsBeingProcessed
)
{
    correctorsBeingProcessed.append(correctori);
    DynamicList<solveID> newCorrectorMembers;
    forAll(correctorMembers_[correctori], memberj)
    {
        const label correctorj =
            correctorNames_.find
            (
                correctorMembers_[correctori][memberj].first()
            );
        if (correctorj >= 0)
        {
            if (correctorsBeingProcessed.found(correctorj))
            {
                wordList cycleNames(correctorsBeingProcessed.size());
                forAll(cycleNames, i)
                {
                    cycleNames[i] =
                        correctorNames_[correctorsBeingProcessed[i]];
                }
                FatalErrorInFunction
                    << "Detected circular inclusion of the corrector named "
                    << correctorNames_[correctorj] << ". "
                    << "The correctors involved in the cycle are: "
                    << cycleNames << nl << exit(FatalError);
            }
            substituteCorrectorNames(correctorj, correctorsBeingProcessed);
            // Now correctorj contains no corrector names. Substitute the solve
            // names into the new corrector members list
            newCorrectorMembers.append(correctorMembers_[correctorj]);
        }
        else
        {
            newCorrectorMembers.append(correctorMembers_[correctori][memberj]);
        }
    }
    correctorMembers_[correctori] = newCorrectorMembers;
    correctorsBeingProcessed.resize(correctorsBeingProcessed.size()-1);
}


template<class Type>
void Foam::foamSolve::getCoupledFields
(
    const fvMesh& mesh,
    const objectRegistry& solReg,
    const word& fieldName,
    List<solveID>& coupledFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;
    if (solReg.foundObject<GeoField>(fieldName))
    {
        GeoField& f = solReg.lookupObjectRef<GeoField>(fieldName);
        forAll(f.boundaryField(), patchi)
        {
            fvPatchField<Type>& pf = f.boundaryFieldRef()[patchi];
            if (pf.regionCoupled())
            {
                regionCoupledFvPatchField<Type>& rcpf =
                    refCast<regionCoupledFvPatchField<Type>>(pf);
                const word otherField = rcpf.coupledField();
                word otherRegion = rcpf.neighbourDb().name();
                if (otherRegion == word::null)
                {
                    otherRegion = rcpf.sampleRegion();
                }
                coupledFields.append
                (
                    solveID
                    (
                        otherField,
                        otherRegion
                    )
                );
            }
        }
    }
}


void Foam::foamSolve::gatherCoupledGroups
(
    solveInfoTable& tree,
    const solveID& nodeKey
)
{
    // Recursively gather the combined list of all contiguous coupled solves
    // so that they can be treated as a group

    solveInfo& node = *(tree[nodeKey]);
    if (node.flag == 0)
    {
        // Node not visited yet

        node.flag = 2;

        // Search all other coupled fields and keep a record of this
        // contiguously connected group of fields
        DynamicList<solveID> coupledGroup;
        coupledGroup.append(nodeKey);
        forAll(node.coupledFields, i)
        {
            // Coupling is bi-directional, so don't come back to one we already
            // started
            solveInfo& otherInfo =
                dynamic_cast<solveInfo&>(*(tree[node.coupledFields[i]]));
            const label otherVisited = otherInfo.flag;
            if (otherVisited != 2)
            {
                gatherCoupledGroups(tree, node.coupledFields[i]);
                // Only need to add to the group if the other node was a new
                // discovery
                if (otherVisited == 0)
                {
                    coupledGroup.append(otherInfo.coupledGroup);
                }
            }
        }
        node.coupledGroup = coupledGroup;

        // Spread the complete information to the others in the coupled group
        // we encountered so far
        forAll(node.coupledGroup, i)
        {
            const solveID& otherKey = node.coupledGroup[i];
            if (otherKey != nodeKey)
            {
                solveInfo& otherNode =
                    dynamic_cast<solveInfo&>(*(tree[otherKey]));
                otherNode.coupledGroup = node.coupledGroup;
            }
        }

        // Mark as completed
        node.flag = 1;
        return;
    }
    else if (node.flag == 1)
    {
        // Already done
        return;
    }
    else
    {
        // Started from this node and haven't got back yet
        FatalErrorInFunction
            << "Unexpected circular dependencies detected involving "
            << "field " << nodeKey.first()
            << " in region " << nodeKey.second() << exit(FatalError);
    }
}



void Foam::foamSolve::nestCorrectors()
{
    // Try to nest the correctors
    List<HashSet<solveID, solveID::hash>> correctorSolves
    (
        correctorNames_.size()
    );
    List<HashSet<word>> correctorCorrectors(correctorNames_.size());
    forAll(correctorNames_, corri)
    {
        correctorSolves[corri].insert(correctorMembers_[corri]);
    }
    // For each corrector, see if it can entirely replace the solves within any
    // of the other correctors (doing these in reverse just to give the more
    // natural result that the correctors specified first end up on the outside
    // in the edge case where they correspond exactly)
    forAllReverse(correctorNames_, corri)
    {
        forAllReverse(correctorNames_, corrj)
        {
            if (corrj == corri)
            {
                continue;
            }
            bool allItemsFound = true;
            forAllConstIters(correctorSolves[corri], iter)
            {
                if (!correctorSolves[corrj].found(iter.key()))
                {
                    allItemsFound = false;
                    break;
                }
            }
            forAllConstIter(HashSet<word>, correctorCorrectors[corri], iter)
            {
                if (!correctorCorrectors[corrj].found(iter.key()))
                {
                    allItemsFound = false;
                    break;
                }
            }
            if (allItemsFound)
            {
                // Substitute all members for the matched corrector
                // in our intermediate variables
                correctorSolves[corrj].erase(correctorSolves[corri]);
                correctorCorrectors[corrj].erase(correctorCorrectors[corri]);
                correctorCorrectors[corrj].insert(correctorNames_[corri]);
            }
        }
    }

    // Nest the solve schedule

    // Create schedule items for all the correctors
    // These are initially created in the top-level schedule and later moved if
    // they are nested inside other correctors
    scheduleUTable& topLevelItems = topLevelCorrector_.items;
    forAll(correctorNames_, corri)
    {
        correctorInfo *newItem = new correctorInfo(correctorNames_[corri]);
        if (!allCorrectors_.insert(correctorNames_[corri], newItem))
        {
            delete newItem;
            FatalErrorInFunction
                << "Duplication of corrector " << correctorNames_[corri]
                << nl << exit(FatalError);
        }
        else
        {
            topLevelItems.insert
            (
                solveID(correctorNames_[corri], word::null),
                newItem
            );
        }
    }

    forAll(correctorNames_, corri)
    {
        // Fill corrector schedule items with sub-correctors
        forAllConstIter(HashSet<word>, correctorCorrectors[corri], iter)
        {
            solveID key(iter.key(), word::null);
            if (!topLevelItems.found(key))
            {
                FatalErrorInFunction
                    << "Circular reference found in correctors "
                    << "while processing corrector " << correctorNames_[corri]
                    << " with members: " << correctorMembers_[corri]
                    << nl << exit(FatalError);
            }
            scheduleUTable::iterator removeIter = topLevelItems.find(key);
            if
            (
                !allCorrectors_[correctorNames_[corri]]->items.insert
                (
                    key, removeIter()
                )
            )
            {
                FatalErrorInFunction
                    << "Unexpected duplication of corrector."
                    << nl << exit(FatalError);
            }
            topLevelItems.erase(removeIter);
        }
    }

    // Go through and link the solves into their corresponding correctors
    // First create pointers in the top-level corrector, then move
    forAllConstIter(solveInfoTable, allSolves_, iter)
    {
        topLevelItems.insert(iter.key(), iter());
    }
    forAll(correctorNames_, corri)
    {
        forAllConstIters(correctorSolves[corri], iter)
        {
            if (!topLevelItems.found(iter.key()))
            {
                FatalErrorInFunction
                    << "Correctors could not be nested. "
                    << "The problem is related to solve " << iter.key()
                    << " and corrector " << correctorNames_[corri] << "." << nl;
                FatalErrorInFunction
                    << "    The correctors were specified as follows:" << nl;
                forAll(correctorNames_, corrj)
                {
                    FatalErrorInFunction
                    << "    " << correctorNames_[corrj] << ": "
                    << flatOutput(correctorMembers_[corrj]) << nl;
                }
                FatalErrorInFunction
                    << nl << exit(FatalError);
            }
            scheduleUTable::iterator removeIter =
                topLevelItems.find(iter.key());
            if
            (
                !allCorrectors_[correctorNames_[corri]]->items.insert
                (
                    iter.key(), removeIter()
                )
            )
            {
                FatalErrorInFunction
                    << "Unexpected duplication of corrector."
                    << nl << exit(FatalError);
            }
            allSolves_[iter.key()]->parentCorrectorName =
                correctorNames_[corri];
            topLevelItems.erase(removeIter);
        }
    }
}


void Foam::foamSolve::addCoupledGroupDependencies
(
    correctorInfo& corr
)
{
    // Make every solve depend on the entire coupled group of each
    // of its dependencies, in order to treat coupled groups as one thing from a
    // dependency point of view

    forAllConstIter(scheduleUTable, corr.items, node)
    {
        if (node()->isCorrector())
        {
            addCoupledGroupDependencies(dynamic_cast<correctorInfo&>(*node()));
        }
        else
        {
            solveInfo& solveItem = dynamic_cast<solveInfo&>(*node());
            solveHashSet dependencies(solveItem.dependencies);
            forAll(solveItem.dependencies, k)
            {
                solveInfo& depInfo = *allSolves_[solveItem.dependencies[k]];
                dependencies.insert(depInfo.coupledGroup);
            }
            // To be safe, make sure we haven't made it depend on itself
            dependencies.erase
            (
                solveID(solveItem.solveName, solveItem.regionName)
            );
            solveItem.dependencies = dependencies.toc();
        }
    }
}


void Foam::foamSolve::addCoupledGroupNestingDependencies
(
    correctorInfo& corr
)
{
    // The solves within coupled groups need to go from least-nested to most-
    // nested since only the final one actually represents the solve; the
    // intermediate ones only consist of the assembly step. We insert
    // dependencies to ensure this.
    // For the same reason, each subsequent solve in the group must be nested
    // in the same or a deeper level of corrector than the previous ones,
    // otherwise a corrector will terminate without the solve actually happening

    forAllConstIter(scheduleUTable, corr.items, node)
    {
        if (node()->isCorrector())
        {
            addCoupledGroupNestingDependencies(dynamic_cast<correctorInfo&>(*node()));
        }
        else
        {
            // Go through coupled group members. If any one is not in this
            // corrector or one nested under it, then check if the current
            // corrector is nested under its parent corrector; if not: fail.
            // If so: add a dependency to ensure it is scheduled first

            solveInfo* solveItem = dynamic_cast<solveInfo*>(node());
            forAll(solveItem->coupledGroup, k)
            {
                solveInfo& otherItem = *allSolves_[solveItem->coupledGroup[k]];
                if (&otherItem != solveItem)
                {
                    solveID foundItem;
                    if (!corr.find(solveItem->coupledGroup[k], foundItem))
                    {
                        if
                        (
                            allCorrectors_[otherItem.parentCorrectorName]->find
                            (
                                solveID
                                (
                                    corr.correctorName, word::null
                                ),
                                foundItem
                            )
                        )
                        {
                            solveItem->dependencies.append
                            (
                                solveItem->coupledGroup[k]
                            );
                        }
                        else
                        {
                            FatalErrorInFunction
                                << "Solves in region-coupled groups must be in "
                                << "nested correctors, but solve "
                                << solveItem->solveName
                                << " in region " << solveItem->regionName << nl
                                << "is in corrector " << corr.correctorName
                                << " which is not nested under corrector "
                                << otherItem.parentCorrectorName << nl
                                << "which contains the coupled solve "
                                << otherItem.solveName << " in region "
                                << otherItem.regionName << nl
                                << exit(FatalError);
                        }
                    }
                }
            }

        }
    }
}


void Foam::foamSolve::constructSolveSchedule
(
    correctorInfo& corr
)
{
    // Construct variable and iteration schedule for a particular corrector

    // First, pre-sort corrector items by predefined ordering according to input
    // files - dependencies will later override this where necessary.
    // Point from list into table so that we can sort;
    // can't use UPtrList for sort - use List of pointers first
    List<scheduleItem*> sortSolveSchedule(corr.items.size());
    label i = 0;
    forAllConstIter(scheduleUTable, corr.items, node)
    {
        if (node()->isCorrector())
        {
            constructSolveSchedule(dynamic_cast<correctorInfo&>(*node()));
        }
        sortSolveSchedule[i++] = node();
    }

    scheduleItem::compareFunc compare;
    sort(sortSolveSchedule, compare);

    // Now work out dependencies local to this corrector.
    // This transforms global dependencies in 'dependencies' member to
    // dependencies local to this corrector in 'localDependencies' member,
    // assigning them to correctors where applicable
    forAll(sortSolveSchedule, j)
    {
        scheduleItem& item = *sortSolveSchedule[j];
        item.flag = false;
        solveInfo* solveItem = dynamic_cast<solveInfo*>(&item);
        solveHashSet dependencies;
        if (solveItem)
        {
            dependencies.insert(solveItem->dependencies);
        }
        else
        {
            // Also make it depend on any unresolved dependencies already moved
            // to correctors
            correctorInfo& corrItem = dynamic_cast<correctorInfo&>(item);
            dependencies.insert(corrItem.localDependencies);
            corrItem.localDependencies.clear();
        }

        forAllConstIters(dependencies, dep)
        {
            if (corr.items.found(dep()))
            {
                item.localDependencies.append(dep());
            }
            else
            {
                // Any dependency that isn't found in the current corrector
                // either points to a solve nested inside a sub-corrector, or
                // inside a parent corrector that we are nested in.
                // In the first case, find out which corrector and replace with
                // a link to it:
                solveID foundIn;
                if (corr.find(dep(), foundIn))
                {
                    item.localDependencies.append(foundIn);
                }
                else
                {
                    // In the second case, move the dependencies to the parent
                    // corrector
                    corr.localDependencies.append(dep());
                }
            }
        }
    }

    // Change ordering to respect dependencies:
    // Instead of walking the tree, we repeatedly loop the list looking for the
    // first allowable solve so as to preserve the sorting as far as possible

    corr.schedule.clear();
    corr.schedule.resize(sortSolveSchedule.size());
    forAll(corr.schedule, i)
    {
        bool foundOne = false;
        forAll(sortSolveSchedule, j)
        {
            scheduleItem& item = *sortSolveSchedule[j];
            if (!item.flag)
            {
                // Check if all dependencies have been done already
                bool dependenciesDone = true;
                forAll(item.localDependencies, k)
                {
                    if (!corr.items[item.localDependencies[k]]->flag)
                    {
                        dependenciesDone = false;
                        break;
                    }
                }
                if (dependenciesDone)
                {
                    foundOne = true;
                    corr.schedule.set(i, &item);
                    item.flag = true;
                    break;
                }
            }
        }
        if (!foundOne)
        {
            FatalErrorInFunction
                << "Unable to resolve solution schedule for corrector "
                << corr.correctorName << "." << nl
                << "This probably indicates circular "
                << "dependencies among the following items:" << nl;
            forAll(sortSolveSchedule, j)
            {
                scheduleItem& item = *sortSolveSchedule[j];
                if (!item.flag)
                {
                    if (item.isCorrector())
                    {
                        FatalErrorInFunction
                            << "Corrector "
                            << dynamic_cast<correctorInfo&>(item).correctorName;
                    }
                    else
                    {
                        solveInfo& solveItem = dynamic_cast<solveInfo&>(item);
                        FatalErrorInFunction
                            << "Solve "
                            <<  solveID
                                (
                                    solveItem.solveName,
                                    solveItem.regionName
                                );
                    }
                    FatalErrorInFunction
                        << " [depends on: "
                        << flatOutput(item.localDependencies)
                        << "]" << nl;
                }
            }
            FatalErrorInFunction << exit(FatalError);
        }
    }
}


void Foam::foamSolve::setScheduleOrder(correctorInfo& corr)
{
    forAll(corr.schedule, stepi)
    {
        if (stepi == 0)
        {
            corr.schedule[stepi].scheduleOrder = corr.scheduleOrder;
        }
        else
        {
            corr.schedule[stepi].scheduleOrder =
                corr.schedule[stepi-1].lastSolveItem().scheduleOrder+1;
        }

        if (corr.schedule[stepi].isCorrector())
        {
            setScheduleOrder
            (
                dynamic_cast<correctorInfo&>(corr.schedule[stepi])
            );
        }
    }
}


void Foam::foamSolve::processSolveSchedule(correctorInfo& corr)
{
    // Record if all required matrices have been assembled
    // at this point for region-coupled solve
    forAll(corr.schedule, stepi)
    {
        if (corr.schedule[stepi].isCorrector())
        {
            processSolveSchedule
            (
                dynamic_cast<correctorInfo&>(corr.schedule[stepi])
            );
        }
        else
        {
            solveInfo& infoi = dynamic_cast<solveInfo&>(corr.schedule[stepi]);

            infoi.readyToSolve = true;
            forAll(infoi.coupledGroup, j)
            {
                solveInfo& infoj = *allSolves_[infoi.coupledGroup[j]];
                if (infoj.scheduleOrder > infoi.scheduleOrder)
                {
                    infoi.readyToSolve = false;
                    break;
                }
            }
        }
    }

    // See if there are any consecutive steps at this corrector level
    // which do not depend on each other, use the same solver, are in the same
    // consolidated region group, and not split by a coupled solve, which can be
    // assembled using serial threads
    forAll(corr.schedule, stepi)
    {
        if (!corr.schedule[stepi].isCorrector())
        {
            solveInfo& infoi = dynamic_cast<solveInfo&>(corr.schedule[stepi]);
            const label regioni = regionNames_[infoi.regionName];
            if (groupConsolidated_[regionGroup_[regioni]])
            {
                solverObjectList& regionObjectsi = solverObjects_[regioni];
                solverObject& solveri = regionObjectsi[infoi.solverOrder];

                label stepj;
                for (stepj = stepi+1; stepj < corr.schedule.size(); stepj++)
                {
                    if (corr.schedule[stepj].isCorrector())
                    {
                        break;
                    }
                    solveInfo& infoj =
                        dynamic_cast<solveInfo&>(corr.schedule[stepj]);
                    solveInfo& infojm1 =
                        dynamic_cast<solveInfo&>(corr.schedule[stepj-1]);
                    const label regionj = regionNames_[infoj.regionName];
                    solverObjectList& regionObjectsj = solverObjects_[regionj];
                    solverObject& solverj = regionObjectsj[infoj.solverOrder];
                    solveID keyi(infoi.solveName, infoi.regionName);
                    if
                    (
                        infoj.dependencies.found(keyi)
                     || solverj.type() != solveri.type()
                     || regionGroup_[regionj] != regionGroup_[regioni]
                     || regioni == regionj
                     || infojm1.readyToSolve
                    )
                    {
                        break;
                    }
                    // Set to zero for all except the first one in the team
                    // (set below)
                    infoj.nTeamThreads = 0;
                }
                infoi.nTeamThreads = stepj-stepi;
                stepi = stepj-1;
            }
            else
            {
                infoi.nTeamThreads = 1;
            }
        }
    }
}


void Foam::foamSolve::runSchedule(correctorInfo& corr)
{
    // Run the items in schedule that are within the specified corrector
    const word& correctorName = corr.correctorName;
    label corrector = 0;
    bool finalIter = false;
    label stepi = 0;

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(correctorProf, "foamSolve::runCorrector "+correctorName);
    #endif

    while
    (
        correctorName == solverObject::timeLoopCorrectorName
      ? time_().run()
      : !finalIter
    )
    {
        updateSolControls();
        // Update the solution control settings from file
        scalar prevCpuTime(0);
        if (correctorName == solverObject::timeLoopCorrectorName)
        {
            forAll(solnControls_, i)
            {
                solnControls_[i].readControls();
            }

            prevCpuTime = time_().elapsedCpuTime();

            // Update global time controls
            bool adjustTimeStep =
                time_().controlDict().lookupOrDefault("adjustTimeStep", false);
            scalar maxDeltaT =
                time_().controlDict().lookupOrDefault<scalar>
                (
                    "maxDeltaT", GREAT
                );
            scalar minDeltaT =
                time_().controlDict().lookupOrDefault<scalar>
                (
                    "minDeltaT", SMALL
                );

            // Set deltaT
            if (adjustTimeStep)
            {
                // Maximum accepted time step
                scalar dt = GREAT;
                forAll(solverObjects_, i)
                {
                    dt = min(dt, solverObjects_[i].getMaxTimeStep());
                }

                // Limit the increase in deltaT
                scalar deltaTFactor = dt/time_().deltaTValue();
                deltaTFactor = min(min(deltaTFactor, 1.0 + 0.1*deltaTFactor), 1.2);

                time_().setDeltaT
                (
                    min
                    (
                        deltaTFactor*time_().deltaTValue(),
                        maxDeltaT
                    )
                );

                Info<< "deltaT = " <<  time_().deltaTValue() << endl;

                //exit Foam if deltaT is below user secified minDeltaT
                if (time_().deltaTValue() < minDeltaT)
                {
                    FatalErrorInFunction
                        << "minDeltaT = " << minDeltaT << nl
                        << "The computed deltaT is lower than the minDeltaT"
                        << exit(FatalError);
                }

            }

            time_()++;

            Info<< nl << nl
                << "########################################"
                << "########################################"
                << nl << endl;
            Info<< "Time = " << time_().timeName() << nl << endl;

            Info<< "Outer Loop Iteration = "
                 << time_().timeIndex() << endl;

            forAll(solutionRegistries_, regi)
            {
                coordinateFrame::updateStates
                (
                    solutionRegistries_[regi].mesh().thisDb()
                );
            }
        }

        finalIter = true;
        bool converged = false;
        forAll(solverObjects_, i)
        {
            bool satisfied =
                solverObjects_[i].isFinalCorrector(corrector, correctorName);

            // If any convergence controls are specified, terminate outer
            // correctors if convergence is detected
            if (correctorName == solverObject::outerCorrectorName && corrector > 0)
            {
                bool conv = false;
                if (isA<pimpleControl>(solnControls_[i]))
                {
                    conv = allConverged(i);
                }

                converged |= (conv && !satisfied);
                satisfied |= conv;
                converged &= satisfied;
            }
            else
            {
                converged = false;
            }

            finalIter &= satisfied;
        }
        if (converged && finalIter)
        {
            Info<< nl << "Residual convergence tolerance reached" << endl;
        }
        if
        (
            correctorName == solverObject::timeLoopCorrectorName
         && time_() > time_().endTime() - 0.5*time_().deltaT()
        )
        {
            finalIter = true;
        }

        if
        (
            correctorName == solverObject::outerCorrectorName
         && !(corrector == 0 && finalIter)
        )
        {
            Info<< nl << "Outer corrector: " << corrector << nl << endl;
        }

        if (debug)
        {
            Info<< "\nIteration: " << corrector << nl
                << "    of corrector: " << correctorName << nl
                << "    final iter: " << Switch(finalIter) << nl << endl;
        }
        forAll(solverObjects_, regioni)
        {
            fvMesh& mesh =  solutionRegistries_[regioni].mesh();
            if (correctorName == solverObject::outerCorrectorName)
            {
                if (isA<pimpleControl>(solnControls_[regioni]))
                {
                    if (finalIter)
                    {
                        mesh.data::add("finalIteration", true);
                    }
                    else
                    {
                        mesh.data::remove("finalIteration");
                    }
                }
                else if (isA<foamCoupledControl>(solnControls_[regioni]))
                {
                    const foamCoupledControl& hcc =
                        refCast<foamCoupledControl>(solnControls_[regioni]);

                    if (hcc.nCorr() > 1)
                    {
                        if (finalIter)
                        {
                            mesh.data::add("finalIteration", true);
                        }
                        else
                        {
                            mesh.data::remove("finalIteration");
                        }
                    }
                }
                solnControls_[regioni].storePrevIterFields();
            }

            solverObjects_[regioni].beginIteration
            (
                corrector,
                correctorName,
                finalIter
            );
            if (corrector == 0)
            {
                mesh.data::set("firstIteration", true);
            }
            else
            {
                mesh.data::remove("firstIteration");
            }
        }

        for (stepi = 0; stepi < corr.schedule.size(); stepi++)
        {
            if (corr.schedule[stepi].isCorrector())
            {
                runSchedule
                (
                    dynamic_cast<correctorInfo&>(corr.schedule[stepi])
                );
            }
            else
            {
                solveInfo& info =
                    dynamic_cast<solveInfo&>(corr.schedule[stepi]);

                // nTeamThreads is positive for first solver in team, and zero
                // for the others, so this block only gets executed for the
                // first one
                const label threadTeamLen = info.nTeamThreads;
                if (threadTeamLen > 0)
                {
                    #if !defined( WIN32 ) && !defined( WIN64 )
                    addProfiling
                    (
                        assembleProf, "foamSolve::assemble " + info.solveName
                    );
                    #endif

                    if (debug && threadTeamLen > 1)
                    {
                        Info<< "Activating team of "
                            << threadTeamLen << " threads" << endl;
                    }

                    if (threadTeamLen > 1)
                    {
                        // assemble gets called in each thread
                        // sequentially
                        serialThreads::run
                        (
                            threadTeamLen,
                            *this,
                            &foamSolve::assemble,
                            corr,
                            static_cast<label>(stepi),
                            static_cast<bool>(finalIter)
                        );
                    }
                    else
                    {
                        assemble(corr, stepi, finalIter);
                    }
                }

                const label regioni = regionNames_[info.regionName];
                const fvMesh& mesh =  solutionRegistries_[regioni].mesh();
                solverObjectList& regionObjects =
                    solverObjects_[regioni];
                solverObject& solver = regionObjects[info.solverOrder];

                // Check if need to wait for further matrices before coupled
                // solve

                // No coupled dependencies: solve segregated
                if (info.coupledFields.empty())
                {
                    {
                        #if !defined( WIN32 ) && !defined( WIN64 )
                        addProfiling
                        (
                            solveProf, "foamSolve::solve " + info.solveName
                        );
                        #endif

                        const word solverDictName =
                            info.solverDictName
                            + (info.solverDictFinal ? "Final" : "");

                        if (info.smx.valid())
                        {
                            printRegionName(regioni);
                            info.smx->solve(mesh.solution().solverDict(solverDictName));
                        }
                        else if (info.vmx.valid())
                        {
                            printRegionName(regioni);
                            info.vmx->solve(mesh.solution().solverDict(solverDictName));
                        }
                        else if (info.stmx.valid())
                        {
                            printRegionName(regioni);
                            info.stmx->solve(mesh.solution().solverDict(solverDictName));
                        }
                        else if (info.tmx.valid())
                        {
                            printRegionName(regioni);
                            info.tmx->solve(mesh.solution().solverDict(solverDictName));
                        }
                        else if (info.v4mx.valid())
                        {
                            printRegionName(regioni);
                            info.v4mx->solve(mesh.solution().solverDict(solverDictName));
                        }
                        info.smx.clear();
                        info.vmx.clear();
                        info.stmx.clear();
                        info.tmx.clear();
                        info.v4mx.clear();
                    }
                    {
                        #if !defined( WIN32 ) && !defined( WIN64 )
                        addProfiling
                        (
                            corrProf, "foamSolve::correct " + info.solveName
                        );
                        #endif
                        solver.correctBase
                        (
                            info.solveName,
                            info.regionName
                            //groupNames_[regionGroup_[regioni]]
                        );
                    }
                }
                else
                {
                    // Try to combine all the dependent matrices and solve
                    // coupled if all available
                    DynamicList<tmp<fvScalarMatrix>> smxs;
                    DynamicList<tmp<fvVectorMatrix>> vmxs;
                    DynamicList<tmp<fvSymmTensorMatrix>> stmxs;
                    DynamicList<tmp<fvTensorMatrix>> tmxs;
                    DynamicList<tmp<fvBlockMatrix<vector4>>> v4mxs;

                    if (info.readyToSolve)
                    {
                        word solverDictName;
                        if
                        (
                            !collectCoupledMatrices
                            (
                                corr,
                                info,
                                smxs, vmxs, stmxs, tmxs, v4mxs,
                                solverDictName
                            )
                        )
                        {
                            FatalErrorInFunction
                                << "All required coupled matrices should exist,"
                                << " but don't."
                                << endl << exit(FatalError);
                        }

                        #if !defined( WIN32 ) && !defined( WIN64 )
                        addProfiling
                        (
                            solveProf,
                            "foamSolve::solveRegionCoupled " + solverDictName
                        );
                        #endif

                        if
                        (
                            (!smxs.empty()) + (!vmxs.empty())
                          + (!stmxs.empty()) + (!tmxs.empty()) > 1
                        )
                        {
                            FatalErrorInFunction
                                << "Currently, region-coupling is supported "
                                << "only between the same type of fields"
                                << nl << endl << exit(FatalError);
                        }

                        if (!smxs.empty())
                        {
                            callCoupledSolver(smxs, solverDictName);
                        }
                        else if (!vmxs.empty())
                        {
                            callCoupledSolver(vmxs, solverDictName);
                        }
                        else if (!stmxs.empty())
                        {
                            callCoupledSolver(stmxs, solverDictName);
                        }
                        else if (!tmxs.empty())
                        {
                            callCoupledSolver(tmxs, solverDictName);
                        }

                        // Call correct on each involved field
                        // First do those in unconsolidated group, then activate
                        // threads for the consolidated ones
                        label nConsolidated = 0;
                        forAll(info.coupledGroup, i)
                        {
                            const solveID& keyi = info.coupledGroup[i];
                            solveInfo& infoi =
                                dynamic_cast<solveInfo&>
                                (
                                    *(allSolves_[keyi])
                                );
                            const label regioni = regionNames_[infoi.regionName];
                            if (groupConsolidated_[regionGroup_[regioni]])
                            {
                                nConsolidated++;
                            }
                            else
                            {
                                #if !defined( WIN32 ) && !defined( WIN64 )
                                addProfiling
                                (
                                    corrProf,
                                    "foamSolve::correct " + info.solveName
                                );
                                #endif

                                solverObjectList& regionObjects =
                                    solverObjects_[regioni];
                                solverObject& solver =
                                    regionObjects[infoi.solverOrder];
                                solver.correct
                                (
                                    infoi.solveName,
                                    //groupNames_[regionGroup_[regioni]]
                                    infoi.regionName
                                );
                            }
                        }
                        if (nConsolidated == 1)
                        {
                            #if !defined( WIN32 ) && !defined( WIN64 )
                            addProfiling
                            (
                                corrProf,
                                "foamSolve::correct " + info.solveName
                            );
                            #endif

                            label groupIdx = 0;
                            correctNext(corr, groupIdx, info);
                        }
                        else if (nConsolidated > 1)
                        {
                            #if !defined( WIN32 ) && !defined( WIN64 )
                            const label regioni =
                                regionNames_[info.regionName];
                            addProfiling
                            (
                                corrProf,
                                "foamSolve::correctConsolidated "
                              + groupNames_[regionGroup_[regioni]]
                            );
                            #endif

                            label groupIdx = 0;
                            serialThreads::run
                            (
                                nConsolidated,
                                *this,
                                &foamSolve::correctNext,
                                corr,
                                static_cast<label&>(groupIdx),
                                const_cast<const solveInfo&>(info)
                            );
                        }
                    }
                }
            }
        }

        forAll(solverObjects_, regioni)
        {
            solverObjects_[regioni].endIteration
            (
                corrector, correctorName, finalIter
            );
        }

        if (correctorName == solverObject::outerCorrectorName)
        {
            // add correct() call for fvOptions
            forAll(fvOptions_, regioni)
            {
                fvOptions_[regioni].correct();
            }
        }

        if (correctorName == solverObject::timeLoopCorrectorName)
        {
            time_().write();

            scalar cpuTime = time_().elapsedCpuTime();
            Info<< nl << "ExecutionTime = " << cpuTime << " s"
                << "  ExecutionStepTime = " << cpuTime-prevCpuTime << " s"
                << "  ClockTime = " << time_().elapsedClockTime() << " s"
                << nl << endl;

            // Check residualControls
            bool converged = true;
            forAll(solnControls_, regioni)
            {
                if (isA<simpleControl>(solnControls_[regioni]))
                {
                    converged &= allConverged(regioni);
                    if (converged && solnControls_[regioni].adjointOptimization())
                    {
                        forAll(solverObjects_, i)
                        {
                            solverObjects_[i].updateTimings();
                        }
                        converged = false;
                    }
                }
                else
                {
                    converged = false;
                }
            }
            if (converged)
            {
                Info<< "Residual convergence tolerance reached" << nl << endl;
                time_().writeAndEnd();
            }
        }

        corrector++;
    }
}


void Foam::foamSolve::assemble
(
    correctorInfo& corr,
    const label stepi,
    const bool finalIter
)
{
    const label stepj =
        serialThreads::active() ? stepi+serialThreads::thisThreadNum() : stepi;

    solveInfo& infoj = dynamic_cast<solveInfo&>(corr.schedule[stepj]);
    const label regionj = regionNames_[infoj.regionName];
    fvMesh& meshj = solutionRegistries_[regionj].mesh();
    solverObjectList& regionObjects = solverObjects_[regionj];
    solverObject& solverj = regionObjects[infoj.solverOrder];

    bool finalSolve =
        meshj.data::template lookupOrDefault<bool>
        (
            "finalIteration",
            false
        );
    word dictName = word::null;

    if (isA<foamCoupledControl>(solnControls_[regionj]))
    {
        const foamCoupledControl& hcc =
        refCast<foamCoupledControl>(solnControls_[regionj]);
        if (hcc.nCorr() == 1)
        {
            meshj.data::remove("finalIteration");
        }
    }

    tmp<fvScalarMatrix> smx =
        solverj.assembleScalarMatrixBase(infoj.solveName, finalSolve, dictName);
    tmp<fvVectorMatrix> vmx =
        solverj.assembleVectorMatrixBase(infoj.solveName, finalSolve, dictName);
    tmp<fvSymmTensorMatrix> stmx =
        solverj.assembleSymmTensorMatrixBase
        (
            infoj.solveName, finalSolve, dictName
        );
    tmp<fvTensorMatrix> tmx =
        solverj.assembleTensorMatrixBase(infoj.solveName, finalSolve, dictName);
    tmp<fvBlockMatrix<vector4>> v4mx = solverj.assembleVector4MatrixBase
        (
            infoj.solveName, finalSolve, dictName
        );

    if (smx.valid())
    {
        infoj.smx = smx;
        if (dictName == word::null)
        {
            dictName = infoj.smx().psi().name();
        }
    }
    else
    {
        infoj.smx.clear();
    }
    if (vmx.valid())
    {
        infoj.vmx = vmx;
        if (dictName == word::null)
        {
            dictName = infoj.vmx().psi().name();
        }
    }
    else
    {
        infoj.vmx.clear();
    }
    if (stmx.valid())
    {
        infoj.stmx = stmx;
        if (dictName == word::null)
        {
            dictName = infoj.stmx().psi().name();
        }
    }
    else
    {
        infoj.stmx.clear();
    }
    if (tmx.valid())
    {
        infoj.tmx = tmx;
        if (dictName == word::null)
        {
            dictName = infoj.tmx().psi().name();
        }
    }
    else
    {
        infoj.tmx.clear();
    }
    if (v4mx.valid())
    {
        infoj.v4mx = v4mx;
        if (dictName == word::null)
        {
            dictName = infoj.v4mx().psi().name();
        }
    }
    else
    {
        infoj.v4mx.clear();
    }
    infoj.finalIter = finalIter;
    infoj.solverDictName = dictName;
    infoj.solverDictFinal = finalSolve;

    if
    (
        label(infoj.smx.valid()) + label(infoj.vmx.valid())
      + label(infoj.stmx.valid()) + label(infoj.tmx.valid())
      + label(infoj.v4mx.valid()) > 1
    )
    {
        FatalErrorInFunction
            << "Multiple matrices of different types returned for "
            << "field " << infoj.solveName
            << nl << endl << exit(FatalError);
    }
}


void Foam::foamSolve::correctNext
(
    correctorInfo& corr, label& groupIdx, const solveInfo& info
)
{
    while (groupIdx < info.coupledGroup.size())
    {
        const solveID& keyi = info.coupledGroup[groupIdx];
        const solveInfo& infoi =
            dynamic_cast<solveInfo&>(*(allSolves_[keyi]));
        const label regioni = regionNames_[infoi.regionName];
        if (groupConsolidated_[regionGroup_[regioni]])
        {
            solverObjectList& regionObjects = solverObjects_[regioni];
            solverObject& solver =
                regionObjects[infoi.solverOrder];
            groupIdx++;
            solver.correctBase
            (
                infoi.solveName,
                infoi.regionName
            );
            break;
        }
        else
        {
            groupIdx++;
        }
    }
}


bool Foam::foamSolve::collectCoupledMatrices
(
    correctorInfo& corr,
    solveInfo& info,
    DynamicList<tmp<fvScalarMatrix>>& smxs,
    DynamicList<tmp<fvVectorMatrix>>& vmxs,
    DynamicList<tmp<fvSymmTensorMatrix>>& stmxs,
    DynamicList<tmp<fvTensorMatrix>>& tmxs,
    DynamicList<tmp<fvBlockMatrix<vector4>>>& v4mxs,
    word& solverDictName
)
{
    // Collect all coupled matrices in the group if they exist. They will no
    // longer be needed and will be moved if this is the final iter for all of
    // them. Return success and empty lists if all are empty, but fail if only
    // some are found.

    bool allFound = true;
    bool anyFound = false;
    bool finalIters = true;
    bool solverDictFinal = false;

    forAll(info.coupledGroup, i)
    {
        const solveInfo& infoi =
            dynamic_cast<solveInfo&>(*(allSolves_[info.coupledGroup[i]]));

        bool found =
            infoi.smx.valid() || infoi.vmx.valid()
         || infoi.stmx.valid() || infoi.tmx.valid()
         || infoi.v4mx.valid();

        allFound &= found;
        anyFound |= found;

        finalIters &= infoi.finalIter;

        solverDictFinal |= infoi.solverDictFinal;
    }

    if (!allFound)
    {
        return !anyFound;
    }
    else
    {
        stringList solverDictNames(info.coupledGroup.size());
        forAll(info.coupledGroup, i)
        {
            const solveInfo& infoi =
                dynamic_cast<solveInfo&>
                (
                    *(allSolves_[info.coupledGroup[i]])
                );

            if (infoi.smx.valid())
            {
                smxs.append(tmp<fvScalarMatrix>(infoi.smx, finalIters));
            }
            if (infoi.vmx.valid())
            {
                vmxs.append(tmp<fvVectorMatrix>(infoi.vmx, finalIters));
            }
            if (infoi.stmx.valid())
            {
                stmxs.append(tmp<fvSymmTensorMatrix>(infoi.stmx, finalIters));
            }
            if (infoi.tmx.valid())
            {
                tmxs.append(tmp<fvTensorMatrix>(infoi.tmx, finalIters));
            }
            if (infoi.v4mx.valid())
            {
                v4mxs.append
                (
                    tmp<fvBlockMatrix<vector4>>(infoi.v4mx, finalIters)
                );
            }

            solverDictNames[i] = infoi.solverDictName;
        }

        solverDictName =
            monolithicSolve<scalar>::makeSuperName(solverDictNames);
        if (solverDictFinal)
        {
            solverDictName += "Final";
        }

        return true;
    }
}


template<class Type>
void Foam::foamSolve::callCoupledSolver
(
    UList<tmp<fvMatrix<Type>>> mxs,
    const word& solverDictName
)
{
    UPtrList<fvMatrix<Type>> eqns(mxs.size());
    bool moving = false;
    bool topoChanging = false;
    DynamicList<word> fieldNames;
    forAll(mxs, i)
    {
        eqns.set(i, &mxs[i].ref());
        moving |= mxs[i]->psi().mesh().moving();
        topoChanging |= mxs[i]->psi().mesh().topoChanging();
        word regionName = mxs[i]->psi().db().name();
        if (regionName == word::null)
        {
            regionName = mxs[i]->psi().mesh().name();
        }
        fieldNames.append(regionName);
    }

    Info<< "Regions: " << flatOutput(fieldNames) << ":" << endl;

    solverPerformance perf =
        monolithicSolve<Type>::solve
        (
            eqns,
            globalFvSolution_->solverDict(solverDictName),
            moving,
            topoChanging
        );
}




bool Foam::foamSolve::allConverged(label regioni)
{
    bool allConverged =
        (!solnControls_[regioni].residualControl().empty());
    const dictionary& perfDict =
        solutionRegistries_[regioni].mesh().solverPerformanceDict();

    bool foundAtLeastOneField = false;
    forAllConstIter(dictionary, perfDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldi = solnControls_[regioni].applyToField(variableName);
        if (fieldi != -1)
        {
            foundAtLeastOneField=true;
            solutionControl::fieldData& resControl =
                solnControls_[regioni].residualControl()[fieldi];

            scalar residual;
            scalar initialRes =
                solnControls_[regioni].maxResidual
                (
                    iter(),
                    residual
                );
            resControl.initialResidual = initialRes+ROOTVSMALL;
            allConverged &=
            (
                residual < resControl.absTol
                || residual/resControl.initialResidual
                < resControl.relTol
            );
        }
    }
    if (!foundAtLeastOneField)
    {
        //Residual Control Dict is not empty but does not contain any fields of this region
        return false;
    }
    return allConverged;
}


void Foam::foamSolve::printSolutionSchedule
(
    correctorInfo& corr, label iterLevel, Ostream& stream
)
{
    // Output the solution schedule
    forAll(corr.schedule, i)
    {
        if (corr.schedule[i].isCorrector())
        {
            for (label k = 0; k < iterLevel+1; k++)
            {
                stream << "  ";
            }
            correctorInfo& newCorr =
                dynamic_cast<correctorInfo&>(corr.schedule[i]);
            stream<< "Corrector " << newCorr.correctorName << endl;
            printSolutionSchedule(newCorr, iterLevel+1, stream);
        }
        else
        {
            solveInfo& info = dynamic_cast<solveInfo&>(corr.schedule[i]);
            for (label j = 0; j < iterLevel+1; j++)
            {
                stream << "  ";
            }

            if (info.nTeamThreads != 1)
            {
                stream << "|-";
            }

            if
            (
                info.readyToSolve
             && info.coupledGroup.size() <= 1
            )
            {
                stream << "Solve ";
            }
            else
            {
                stream << "Assemble ";
            }
            stream << info.solveName << " in region "
                << info.regionName;
            if (debug)
            {
                stream << " [" << flatOutput(info.dependencies)
                    << "," << info.regionOrder
                    << "," << info.solverOrder
                    << "," << info.fieldOrder
                    << "]";
            }
            stream << endl;
            if
            (
                info.readyToSolve
             && info.coupledGroup.size() > 1
            )
            {
                for (label j = 0; j < iterLevel+1; j++)
                {
                    stream << "  ";
                }
                stream << "Solve coupled: ( ";
                forAll(info.coupledGroup, j)
                {
                    stream << info.coupledGroup[j] << " ";
                }
                stream << ")" << endl;
            }
        }
    }
}


void Foam::foamSolve::printRegionName
(
    const label regioni
) const
{
    if (regionNames_.size()>1)
    {
        Info<< "Region: " << regionNames_[regioni] << ":" << endl;
    }
}


Foam::label Foam::foamSolve::getMaxTeamLen(correctorInfo& corr)
{
    label maxTeamLen = 0;
    for (auto& si : corr.schedule)
    {
        if (si.isCorrector())
        {
            maxTeamLen =
                max
                (
                    maxTeamLen,
                    getMaxTeamLen(dynamic_cast<correctorInfo&>(si))
                );
        }
        else
        {
            maxTeamLen =
                max(maxTeamLen, dynamic_cast<solveInfo&>(si).nTeamThreads);
        }
    }
    return maxTeamLen;
}

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

Foam::foamSolve::foamSolve()
:
    topLevelCorrector_(solverObject::timeLoopCorrectorName)
{}

bool Foam::foamSolve::isPostProcessActive(argList& args)
{
    bool postProcess = false;
    if (args.optionFound(argList::postProcessOptionName))
    {
        postProcess = true;
        if (args.optionFound("list"))
        {
            functionObjectList::list();
            return 0;
        }
    }
    return postProcess;
}

void Foam::foamSolve::createTime(const argList& args)
{
    Info<< "Creating time\n" << endl;
    time_.set(new Time(Time::controlDictName, args));
}

Foam::dictionary Foam::foamSolve::getSolverOptionConfiguration
(
    const word& configName
)
{
    Info<< "Reading default options for solver configuration: "
        << configName << nl << endl;
    Time& runTime(time_());
    fileName etcFileName =
        "caseDicts/solvers/foamSolve"/configName/"options.cfg";
    fileName configFileName = findEtcFile(etcFileName);
    dictionary configOptionsDict;
    if (!configFileName.empty())
    {
        configFileName.expand();
        configFileName.toAbsolute();
        IOdictionary iod
        (
            IOobject
            (
                configFileName,
                runTime,
                IOobject::MUST_READ
            )
        );
        configOptionsDict = iod;
    }
    else if (configName != "foamSolve")
    {
        FatalErrorInFunction
            << "foamSolve configuration " << configName
            << " was not found" << nl
            << "at " << "$FOAM_ETC"/etcFileName << nl
            << exit(FatalError);
    }
    fileName defaultEtcFileName =
        "caseDicts/solvers/foamSolve/options.cfg";
    fileName defaultConfigFileName = findEtcFile(defaultEtcFileName);
    dictionary defaultOptionsDict;
    if (!defaultConfigFileName.empty())
    {
        defaultConfigFileName.expand();
        defaultConfigFileName.toAbsolute();
        IOdictionary iod
        (
            IOobject
            (
                defaultConfigFileName,
                runTime,
                IOobject::MUST_READ
            )
        );
        defaultOptionsDict = iod;
    }
    configOptionsDict.merge(defaultOptionsDict, false);
    return configOptionsDict;
}
// * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * * //
void Foam::foamSolve::runPostProcess
(
    argList& args,
    const word& configName
)
{
    Time& runTime(time_());
    dictionary configOptionsDict = getSolverOptionConfiguration(configName);
    const bool lowMem = args.optionFound("lowMemory");

    Foam::instantList timeDirs = Foam::timeSelector::select0
    (
        runTime, args
    );

    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary functionsDict;

    HashSet<word> selectedFields;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;
        while (scheduler_->loop())
        {
            runTime.setTime
            (
                runTime.startTime().value(),
                runTime.startTimeIndex()
            );
            solutionRegistries_.clear();
            const PtrList<fvSolutionRegistry>& regList =
                scheduler_->getActiveSolutionRegistryList(runTime);
            solutionRegistries_.resize
            (
                scheduler_->getCurrentInstance().regionNames().size()
            );
            forAll(scheduler_->getCurrentInstance().regionNames(), i)
            {
                solutionRegistries_.set(i, &regList[i]);
            }

            // Get regions
            getRegions();


            FatalIOError.throwExceptions();

            try
            {
                forAll(solutionRegistries_,i)
                {
                    auto& m = solutionRegistries_[i].mesh();
                    m.readUpdate();
                }

                createOptions(configOptionsDict);

                // Re-generate functionObjects
                autoPtr<functionObjectList> functionsPtr =
                    functionObjectList::New
                    (
                        args,
                        runTime,
                        functionsDict,
                        selectedFields
                    );

                functionsPtr->execute(lowMem);

                // Execute the functionObject 'end()'
                // function for the last time
                if (!lowMem && timei == timeDirs.size()-1)
                {
                    functionsPtr->end();
                }
                // Destroy everything (except meshes) to force re-reading
                // from disk at next time
                fvOptions_.clear();

                // Additionally remove any objects owned by the object
                // registry to force re-reading
                forAll(solutionRegistries_, meshi)
                {
                    bool removedSomething;
                    do
                    {
                        removedSomething = false;
                        HashTable<regIOobject*> objs =
                            solutionRegistries_[meshi].mesh().lookupClass
                            <
                                regIOobject
                            >();
                        forAllIter(HashTable<regIOobject*>, objs, iter)
                        {
                            if (iter()->ownedByRegistry())
                            {
                                iter()->checkOut();
                                // Must re-start the loop as this may
                                // indirectly remove other objects in the
                                // hash table
                                removedSomething = true;
                                break;
                            }
                        }
                    } while (removedSomething);
                }
            }
            catch (IOerror& err)
            {
                Warning<< err << endl;
            }
            scheduler_->operator++();
            Info<< endl;
        }
    }
    Info<< "End\n" << endl;
}

int Foam::foamSolve::runSolver(argList& args, const word& configName)
{
    bool postProcess = isPostProcessActive(args);

    createTime(args);
    Time& runTime(time_());

    dictionary configOptionsDict = getSolverOptionConfiguration(configName);

    globalFvSolution_.set(new fvSolution(runTime));

    Info<<"Creating solution schedule "<<endl;
    scheduler_.set(new solutionScheduler(args, runTime));

    if (postProcess)
    {
        runPostProcess(args, configName);
        return 0;
    }

    while (scheduler_->loop())
    {
        runTime.setTime
        (
            runTime.startTime().value(),
            runTime.startTimeIndex()
        );
        solutionRegistries_.clear();
        const PtrList<fvSolutionRegistry>& regList =
            scheduler_->getActiveSolutionRegistryList(runTime);
        solutionRegistries_.resize
        (
            scheduler_->getCurrentInstance().regionNames().size()
        );
        forAll(scheduler_->getCurrentInstance().regionNames(), i)
        {
            solutionRegistries_.set(i, &regList[i]);
        }

        // Get regions
        getRegions();

        createOptions(configOptionsDict);
        // Look up solutionControls
        solnControls_.resize(solutionRegistries_.size());
        forAll(solutionRegistries_, i)
        {
            solnControls_.set
            (
                i,
                &solutionRegistries_[i].registry().
                    lookupObjectRef<solutionControl>
                    (
                        solutionControl::typeName
                    )
            );
        }

        // Gather all dependency info
        gatherDependencyInfo();

        // Perform sanity checks on the inputted data
        checkSolveInfo();

        // Recursively substitute any corrector names found in the
        // correctorMembers lists with the full contents of that corrector
        forAll(correctorMembers_, correctori)
        {
            DynamicList<label> correctorsBeingProcessed;
            substituteCorrectorNames(correctori, correctorsBeingProcessed);
        }

        // Collect groupings of coupled fields - note that coupled groups
        // are treated as one unit from a dependency point of view, but can span
        // multiple correctors
        forAllConstIter(solveInfoTable, allSolves_, node)
        {
            node()->flag = 0;
        }
        forAllConstIter(solveInfoTable, allSolves_, node)
        {
            if (node()->flag == 0)
            {
                gatherCoupledGroups(allSolves_, node.key());
            }
        }

        // Move the solves into correctors and resolve them into nested levels
        nestCorrectors();

        // Add dependencies to schedule the assembly of different members of
        // coupled groups correctly
        addCoupledGroupDependencies(topLevelCorrector_);
        addCoupledGroupNestingDependencies(topLevelCorrector_);

        // Calculate solution schedule based on assigned dependencies
        // Acts recursively on nested correctors
        constructSolveSchedule(topLevelCorrector_);
        if (topLevelCorrector_.localDependencies.size())
        {
            FatalErrorInFunction
                << "Unable to compute solution schedule - the following "
                << "dependencies could not be resolved:" << nl
                << topLevelCorrector_.localDependencies << nl
                << "Solution schedule:" << nl;
            printSolutionSchedule(topLevelCorrector_, 0, FatalErrorInFunction);
            FatalErrorInFunction << exit(FatalError);
        }
        topLevelCorrector_.scheduleOrder = 0;
        setScheduleOrder(topLevelCorrector_);
        processSolveSchedule(topLevelCorrector_);

        Info<< nl << "Solution schedule:" << nl;
        printSolutionSchedule(topLevelCorrector_, 0, Info);
        Info<< endl;


        // Create thread pool
        label maxTeamLen = getMaxTeamLen(topLevelCorrector_);
        if (maxTeamLen > 1)
        {
            Info<< "Creating pool of " << maxTeamLen
                << " sequential threads" << endl;
            serialThreads::create(maxTeamLen);
        }

        Info<< "\nStarting time loop\n" << endl;

        // Run the time steps
        runSchedule(topLevelCorrector_);

        topLevelCorrector_.clear();
        if (maxTeamLen > 1)
        {
            // Shut down threads
            serialThreads::end();
        }
        scheduler_->operator++();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
