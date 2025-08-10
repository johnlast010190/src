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
    (at your solverObject) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2019-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solverObjectList.H"
#include "solverOption/solverOption.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solverObjectList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solverObjectList::solverObjectList
(
    fv::optionList& fvOptions,
    bool legacySolver
)
:
    UPtrList<solverObject>(),
    solverOptions_(),
    fvOptions_(fvOptions)
{
    forAll(fvOptions, i)
    {
        fv::option& option = fvOptions[i];
        if (isA<fv::solverOption>(option))
        {
            solverObject& solverObj =
                refCast<fv::solverOption>(option).solverObject();
            this->resize(this->size()+1u);
            this->set
            (
                this->size()-1,
                &solverObj
            );
            solverObj.setLegacy(legacySolver);
            solverOptions_.resize(this->size());
            solverOptions_.set
            (
                this->size()-1,
                &refCast<fv::solverOption>(option)
            );
        }
    }
    correctorNames_.clear();
    correctorNames_.resize(this->size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solverObjectList::getMaxTimeStep()
{
    scalar maxTimeStep = GREAT;
    forAll(*this, i)
    {
        solverObject& objectI = this->operator[](i);
        maxTimeStep = min(maxTimeStep, objectI.getMaxTimeStep());
    }
    return maxTimeStep;
}


void Foam::solverObjectList::getSolveGraph
(
    DynamicList<word>& solveNames,
    HashTable<wordList>& derivedFields,
    SolveTable<List<Tuple2<solveID, bool>>>& solveDependencies,
    DynamicList<label>& solverOrder,
    HashTable<wordList>& correctorMembers
)
{
    // Record the originating option name per solve name
    DynamicList<word> optionNames;

    forAll(*this, i)
    {
        solverObject& objectI = this->operator[](i);

        DynamicList<word> sn;
        HashTable<wordList> sf;
        SolveTable<solveList> rds;
        SolveTable<solveList> ods;
        HashTable<wordList> cm;
        objectI.getSolveGraphBase(sn, sf, rds, ods, cm);

        // Merge solveNames
        for (word s : sn)
        {
            if (solveNames.found(s))
            {
                FatalErrorInFunction
                    << "Solve name '" << s << "' in solver object '"
                    << this->operator[](i).name() << "' duplicates another "
                    << "in solver object '" << optionNames[solveNames.find(s)]
                    << "' in region " << this->operator[](i).mesh().name()
                    << exit(FatalError);
            }
        }
        solveNames.append(sn);
        optionNames.append(wordList(sn.size(), this->operator[](i).name()));

        // Merge derivedFields
        forAllConstIters(sf, fields)
        {
            DynamicList<word> fieldNames
            (
                derivedFields.lookup(fields.key(), wordList())
            );
            for (word s : fields())
            {
                if (!fieldNames.found(s))
                {
                    fieldNames.append(s);
                }
            }
            derivedFields.set(fields.key(), fieldNames);
        }

        // Merge dependencies
        forAllConstIters(rds, dep)
        {
            DynamicList<Tuple2<solveID, bool>> depSolves
            (
                solveDependencies.lookup
                (
                    dep.key(), List<Tuple2<solveID, bool>>()
                )
            );
            for (const solveID& s : dep())
            {
                Tuple2<solveID, bool> d(s, true);
                if (!depSolves.found(d))
                {
                    depSolves.append(d);
                }
            }
            solveDependencies.set(dep.key(), depSolves);
        }
        forAllConstIters(ods, dep)
        {
            DynamicList<Tuple2<solveID, bool>> depSolves
            (
                solveDependencies.lookup
                (
                    dep.key(), List<Tuple2<solveID, bool>>()
                )
            );
            for (const solveID& s : dep())
            {
                Tuple2<solveID, bool> d(s, false);
                if (!depSolves.found(d))
                {
                    depSolves.append(d);
                }
            }
            solveDependencies.set(dep.key(), depSolves);
        }
        solverOrder.resize(solveNames.size(), i);

        // Merge correctors
        forAllConstIters(cm, corr)
        {
            DynamicList<word> corrMembers
            (
                correctorMembers.lookup(corr.key(), wordList())
            );
            for (word s : corr())
            {
                if (!corrMembers.found(s))
                {
                    corrMembers.append(s);
                }
            }
            correctorMembers.set(corr.key(), corrMembers);
        }

        // Store corrector names for each solver object to avoid calling into
        // solver objects for correctors that they are not part of
        correctorNames_[i].clear();
        correctorNames_[i].insert
        (
            solverObject::timeLoopCorrectorName
        );
        correctorNames_[i].set(cm.toc());
    }

    // Add any dependencies of the fvOption sources
    SolveTable<solveList> depList;
    fvOptions_.addSourceDependencies(depList);
    // Create if missing and append
    forAllConstIters(depList, iter)
    {
        solveDependencies.insert(iter.key(), List<Tuple2<solveID, bool>>());
        List<Tuple2<solveID, bool>>& deps = solveDependencies[iter.key()];
        forAll(iter(), j)
        {
            deps.append(Tuple2<solveID, bool>(iter()[j], false));
        }
    }
}


bool Foam::solverObjectList::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    bool done = true;
    forAll(*this, i)
    {
        if (correctorNames_[i].found(correctorName))
        {
            solverObject& objectI = this->operator[](i);
            done &= objectI.isFinalCorrectorBase(corrector, correctorName);
        }
    }
    return done;
}


void Foam::solverObjectList::beginIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    forAll(*this, i)
    {
        if (correctorNames_[i].found(correctorName))
        {
            solverObject& objectI = this->operator[](i);
            objectI.beginIterationBase(corrector, correctorName, finalIter);
        }
    }
}


void Foam::solverObjectList::endIteration
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    forAll(*this, i)
    {
        if (correctorNames_[i].found(correctorName))
        {
            solverObject& objectI = this->operator[](i);
            objectI.endIterationBase(corrector, correctorName, finalIter);
        }
    }
}

void Foam::solverObjectList::updateTimings()
{
    forAll(*this, i)
    {
        solverObject& objectI = this->operator[](i);
        objectI.updateSolverObjectTimes();
    }
}

// ************************************************************************* //
