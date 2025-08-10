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

Class
    Foam::foamSolve

Description
    Solver object scheduler and single-/multi-region assembler and solver

SourceFiles
    foamSolve.C

\*---------------------------------------------------------------------------*/

#ifndef foamSolve_H
#define foamSolve_H

#include "db/Time/Time.H"
#include "global/argList/argList.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "primitives/strings/lists/hashedWordList.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "primitives/Pair/PairKey/PairKey.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "solverObject/solverObjectList.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "solutionRegions/solutionScheduler.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
struct scheduleItem;
struct solveInfo;

// Use PairKey as key for (solve, region) [for solves] and
// (corrector, word::null) [for correctors]
typedef HashPtrTable<solveInfo, solveID, solveID::hash>
    solveInfoTable;
typedef HashTable<scheduleItem*, solveID, solveID::hash>
    scheduleUTable;
typedef HashSet<solveID,solveID::hash> solveHashSet;

// Gathers the information for each schedule item - which is
// either a solve or a group of solves in a corrector
struct scheduleItem
{
    // Functor to sort the solves in order to achieve the preferred ordering
    // - ordering to satisfy dependencies is imposed separately
    static class compareFunc
    {
        public:
            bool operator() (const scheduleItem* ps1, const scheduleItem* ps2);
    } compare;

    // Flag for marking as visited/done during tree traversal
    label flag;

    // (solve, region) of dependencies within this corrector
    DynamicList<solveID> localDependencies;

    // Final computed order in the solve schedule
    label scheduleOrder;

    scheduleItem()
    :
        flag(0),
        scheduleOrder(-1)
    {}

    virtual ~scheduleItem() {};

    virtual bool isCorrector() const = 0;

    virtual const solveInfo& firstSolveItem() const = 0;

    virtual const solveInfo& lastSolveItem() const = 0;
};


// Struct which gathers the information for each solve obtained from the solve
// graph of the solver objects
struct solveInfo
:
    public scheduleItem
{
    // Solve name and region of this solve
    word solveName;
    word regionName;
    // Solve name and any extra fields solved for
    wordHashSet solveNameAndFields;
    // Name of the region group and whether it is a consolidated group
    word groupName;
    // Original ordering of region in regionProperties
    label regionOrder;
    // Original ordering of this solver in the fv::option list for this region
    label solverOrder;
    // Original ordering of this field within the fv::option's list of fields
    label fieldOrder;
    // (solve, region) of dependencies (solves that we depend on)
    DynamicList<solveID> dependencies;
    // (field, region) for implicit solver couplings
    List<solveID> coupledFields;
    // (field, region) for the whole group of contiguously coupled fields
    // connected to this one
    List<solveID> coupledGroup;
    // Name of corrector that this solve is part of
    word parentCorrectorName;
    // Whether all coupled matrices have been assembled at this point in the
    // schedule
    bool readyToSolve;
    // Number of threads in this team, if in a consolidated
    // region group and scheduled consecutively
    label nTeamThreads;

    // Storage for matrices while assembling
    tmp<fvScalarMatrix> smx;
    tmp<fvVectorMatrix> vmx;
    tmp<fvSymmTensorMatrix> stmx;
    tmp<fvTensorMatrix> tmx;
    tmp<fvBlockMatrix<vector4>> v4mx;

    // Was matrix stored at the final iter of the iterLevel?
    bool finalIter;

    // Name to use to select solver dict
    word solverDictName;
    // Whether to select final solver dict
    bool solverDictFinal;

    solveInfo()
    :
        scheduleItem(),
        regionOrder(0),
        solverOrder(0),
        fieldOrder(0),
        readyToSolve(false),
        nTeamThreads(1),
        finalIter(false),
        solverDictFinal(false)
    {}

    virtual bool isCorrector() const override { return false; };

    virtual const solveInfo& firstSolveItem() const override
    {
        return *this;
    }

    virtual const solveInfo& lastSolveItem() const override
    {
        return *this;
    }
};


struct correctorInfo
:
    public scheduleItem
{
    word correctorName;

    //- Dependency tree
    scheduleUTable items;
    //- Solution schedule (ordered pointers into 'items')
    UPtrList<scheduleItem> schedule;

    correctorInfo(const word corrName)
    :
        correctorName(corrName)
    {}

    virtual bool isCorrector() const override { return true; };

    // Recursively find the schedule item - returns true if found and sets
    // foundItem to the name an item in this corrector if found (either the
    // item name itself or the name of a sub-corrector that contains it)
    bool find(const solveID& findItem, solveID& foundItem) const
    {
        if (items.found(findItem))
        {
            foundItem = findItem;
            return true;
        }
        forAllConstIter(scheduleUTable, items, item)
        {
            if (item()->isCorrector())
            {
                correctorInfo& corrItem = dynamic_cast<correctorInfo&>(*item());
                if (corrItem.find(findItem, foundItem))
                {
                    foundItem = item.key();
                    return true;
                }
            }
        }
        return false;
    }

    virtual const solveInfo& firstSolveItem() const override
    {
        if (!schedule.size())
        {
            FatalErrorInFunction
                << "Unexpected empty corrector" << nl << exit(FatalError);
        }
        return schedule[0].firstSolveItem();
    }

    virtual const solveInfo& lastSolveItem() const override
    {
        if (!schedule.size())
        {
            FatalErrorInFunction
                << "Unexpected empty corrector" << nl << exit(FatalError);
        }
        return schedule[schedule.size()-1].lastSolveItem();
    }

    virtual void clear()
    {
        schedule.clear();
        items.clear();
    }
};


class foamSolve
{

protected:

    //- Protected data

    autoPtr<Time> time_;

    //- List of all region names
    hashedWordList regionNames_;
    //- List of all group names
    wordList groupNames_;
    //- For a given region, the group number
    labelList regionGroup_;
    //- For a given group, whether it is consolidated
    boolList groupConsolidated_;

    //- Solution registry for each region
    UPtrList<const fvSolutionRegistry> solutionRegistries_;

    //- Option list for each region
    autoPtr<fv::options> globalFvOptions_;

    //- Option list for each region
    UPtrList<fv::options> fvOptions_;
    //- solverObject list for each region
    PtrList<solverObjectList> solverObjects_;
    //- Solution controls for each region
    UPtrList<solutionControl> solnControls_;

    //- Top-level fvSolution
    autoPtr<fvSolution> globalFvSolution_;

    //- Central storage of all the solve information
    solveInfoTable allSolves_;

    //- Map from solve and derived field names to their primary solve name
    HashTable<solveID, solveID, solveID::hash>
        fieldToSolveName_;

    //- Storage of all the correctors
    HashPtrTable<correctorInfo, word, word::hash> allCorrectors_;

    //- Solution schedule at the top level (time steps) - other correctors are
    //  nested
    correctorInfo topLevelCorrector_;

    //- Consolidated list of names of all correctors from all solver objects
    DynamicList<word> correctorNames_;
    //- For each corrector name, list of solves as (field, region)
    DynamicList<solveList> correctorMembers_;

    autoPtr<solutionScheduler> scheduler_;

protected:

    void getRegions();

    void createOptions(const dictionary& configOptionsDict);

    void orderFvOptions(fv::optionList& options);

    void gatherDependencyInfo();

    //- Perform sanity checks on input data
    void checkSolveInfo();

    //- Recursively substitute any corrector names found in the correctorMembers
    //  lists with the full contents of that corrector
    void substituteCorrectorNames
    (
        const label correctori,
        DynamicList<label>& correctorsBeingProcessed
    );

    template<class Type>
    void getCoupledFields
    (
        const fvMesh& mesh,
        const objectRegistry& solReg,
        const word& fieldName,
        List<solveID>& coupledFields
    );

    void gatherCoupledGroups
    (
        solveInfoTable& tree,
        const solveID& nodeKey
    );

    void nestCorrectors();

    void addCoupledGroupDependencies(correctorInfo& corr);
    void addCoupledGroupNestingDependencies(correctorInfo& corr);

    void constructSolveSchedule(correctorInfo& corr);
    void setScheduleOrder(correctorInfo& corr);
    void processSolveSchedule(correctorInfo& corr);

    void runSchedule(correctorInfo& corr);

    void assemble(correctorInfo& corr, const label stepi, const bool finalIter);

    void correctNext
    (
        correctorInfo& corr, label& groupIdx, const solveInfo& info
    );

    bool collectCoupledMatrices
    (
        correctorInfo& corr,
        solveInfo& info,
        DynamicList<tmp<fvScalarMatrix>>& smxs,
        DynamicList<tmp<fvVectorMatrix>>& vmxs,
        DynamicList<tmp<fvSymmTensorMatrix>>& stmxs,
        DynamicList<tmp<fvTensorMatrix>>& tmxs,
        DynamicList<tmp<fvBlockMatrix<vector4>>>& v4mxs,
        word& solverDictName
    );

    template<class Type>
    void callCoupledSolver
    (
        UList<tmp<fvMatrix<Type>>> mxs,
        const word& solverDictName
    );

    bool allConverged(label regioni);

    bool decideTransient(const fvMesh& mesh);

    void printSolutionSchedule
    (
        correctorInfo& corr, label iterLevel, Ostream& stream
    );

    void printRegionName(const label) const;

    label getMaxTeamLen(correctorInfo& corr);

    bool isPostProcessActive(argList& args);

    void createTime(const argList& args);

    dictionary getSolverOptionConfiguration(const word& configName);

public:

    //- Runtime type information
    ClassName("foamSolve");

    foamSolve();

    foamSolve(foamSolve&) = delete;

    void operator=(const foamSolve&) = delete;

    int runSolver(argList& args, const word& configName);

    void runPostProcess(argList& args, const word& configName);

    void updateSolControls()
    {
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
    }

    const HashTable<solveID, solveID, solveID::hash>&
    fieldToSolveName() const
    {
        return fieldToSolveName_;
    }
};


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
