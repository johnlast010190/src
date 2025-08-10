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

#include "solverObject.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "staticFvMesh/staticFvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solverObject::lookupInactiveRegions()
{
    if (obr_.time().found("regionProperties"))
    {
        const IOdictionary& regionProperties =
            obr_.time().lookupObject<IOdictionary>("regionProperties");
        const wordList inactiveRegions
        (
            regionProperties.lookupOrDefault
            (
                "inactiveRegions",
                wordList::null()
            )
        );
        if (inactiveRegions.found(mesh_.name()))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solverObject, 0);

    const word solverObject::timeLoopCorrectorName = "timeLoop";
    const word solverObject::outerCorrectorName = "outerCorrector";
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::solverObject::active() const
{
    if (inactiveRegion_)
    {
        return false;
    }
    label idx = nCorrectorsSinceLastActive_;
    const label solveInterval =
        round(solveInterval_->value(mesh_.time().value()));
    bool validIter = (solveInterval > 0 && (idx % solveInterval  == 0));
    if
    (
        solveInitial_
     && correctorNumber_.lookup(solverObject::outerCorrectorName, 0) == 0
    )
    {
        validIter = true;
    }
    if
    (
        solveFinal_
     && finalIter_.lookup(solverObject::outerCorrectorName, false)
    )
    {
        validIter = true;
    }
    if
    (
        mesh_.time().timeOutputValue() < timeStart_
     || mesh_.time().timeOutputValue() > timeEnd_
    )
    {
        validIter = false;
    }
    if (validIter)
    {
        nCorrectorsSinceLastActive_ = 0;
    }
    return validIter;
}


Foam::fvMesh& Foam::solverObject::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        const fvMesh& mesh = dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
        return const_cast<fvMesh&>(mesh);
    }
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr);
    return const_cast<fvMesh&>(mesh);
}


bool Foam::solverObject::isStatic() const
{
    if (!mesh().hasChangers())
    {
        return (!isA<dynamicFvMesh>(mesh()) || (isA<staticFvMesh>(mesh())));
    }
    else
    {
        return !mesh().dynamic();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solverObject::solverObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, solverObject>
    (
        getMesh(obr),
        name
    ),
    name_(name),
    obr_(obr),
    mesh_(getMesh(obr)),
    regionName_(obr.name() == word::null ? mesh_.name() : obr.name()),
    dict_(dict),
    coeffs_(dict.optionalSubDict(dict.dictName() + "Coeffs")),
    phaseName_(dict.lookupOrDefault("phaseName", word::null)),
    correctorName_(word::null),
    correctorNumber_(),
    finalIter_(),
    timeStart_(-VGREAT),
    timeEnd_(VGREAT),
    solveInterval_(new Function1Types::Constant<scalar>("solveInterval", 1)),
    nCorrectorsSinceLastActive_(0),
    fvOptionsPtr_(nullptr),
    legacySolver_(true),
    legacyOuterCorrector_(0),
    oneExecPerIter_(false),
    lastTimeIndex_(-1),
    verbose_(false),
    inactiveRegion_(lookupInactiveRegions()),
    solveController_(nullptr),
    adjoint_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solverObject::~solverObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::solverObject::setTimeStart(const scalar& timeStart)
{
    timeStart_ = timeStart;
}

void Foam::solverObject::setTimeEnd(const scalar& timeEnd)
{
    timeEnd_ = timeEnd;
}

void Foam::solverObject::readBase(const dictionary& dict)
{
    dict_ = dict;
    coeffs_ = dict.optionalSubDict(dict.dictName() + "Coeffs");
    timeStart_ = dict.lookupOrDefault("timeStart", -VGREAT);
    timeEnd_ = dict.lookupOrDefault("timeEnd", VGREAT);
    solveInitial_ = dict.lookupOrDefault("solveInitial", false);
    solveFinal_ = dict.lookupOrDefault("solveFinal", false);

    if (dict.found("solveInterval"))
    {
        solveInterval_.reset(Function1<scalar>::New("solveInterval", dict));
    }
    else
    {
        solveInterval_.reset
        (
            new Function1Types::Constant<scalar>
            (
                "solveInterval",
                (solveInitial_ || solveFinal_) ? 0 : 1
            )
        );
    }
    oneExecPerIter_ = dict.lookupOrDefault("oneExecPerIter", false);
    verbose_ = dict.lookupOrDefault("verbose", false);

    this->read(dict);
}


bool Foam::solverObject::initialiseBase()
{
    return this->initialise();
}


void Foam::solverObject::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& derivedFields,
    HashTable<solveList>& requiredDependencies,
    HashTable<solveList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Default implementation just adds the current region name to every
    // dependency - override to specify different regions

    HashTable<wordList> rds;
    HashTable<wordList> ods;
    this->getSolveGraph
    (
        solveNames,
        derivedFields,
        rds,
        ods,
        correctorMembers
    );

    forAllConstIters(rds, iter)
    {
        const wordList& rdsi = iter();
        requiredDependencies.insert
        (
            iter.key(), solveList(rdsi.size())
        );
        solveList& requiredDependenciesi
        (
            requiredDependencies[iter.key()]
        );
        forAll(rdsi, j)
        {
            requiredDependenciesi[j] =
                solveID(rdsi[j], regionName_);
        }
    }

    forAllConstIters(ods, iter)
    {
        const wordList& odsi = iter();
        optionalDependencies.insert
        (
            iter.key(), solveList(odsi.size())
        );
        solveList& optionalDependenciesi
        (
            optionalDependencies[iter.key()]
        );
        forAll(odsi, j)
        {
            optionalDependenciesi[j] = solveID(odsi[j], regionName_);
        }
    }
}


void Foam::solverObject::getSolveGraphBase
(
    wordList& solveNames,
    HashTable<wordList>& derivedFields,
    SolveTable<solveList>& requiredDependencies,
    SolveTable<solveList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    if
    (
        solveNames.size() != 0
     || derivedFields.size() != 0
     || requiredDependencies.size() != 0
     || optionalDependencies.size() != 0
     || correctorMembers.size() != 0
    )
    {
        FatalErrorInFunction
            << "Unexpected parameters received."
            << abort(FatalError);
    }

    HashTable<solveList> rDep;
    HashTable<solveList> oDep;
    this->getSolveGraph
    (
        solveNames,
        derivedFields,
        rDep,
        oDep,
        correctorMembers
    );

    forAllConstIters(rDep, iter)
    {
        solveID key(iter.key(), regionName_);
        requiredDependencies.insert(key, iter());
    }
    forAllConstIters(oDep, iter)
    {
        solveID key(iter.key(), regionName_);
        optionalDependencies.insert(key, iter());
    }

    // Add any dependencies specified in user input
    HashTable<wordList> inputDeps
    (
        coeffs_.lookupOrDefault("dependencies", HashTable<wordList>())
    );
    forAllConstIters(inputDeps, dep)
    {
        solveID key(dep.key(), regionName_);
        solveList depSolves(requiredDependencies.lookup(key, solveList()));
        const wordList& newDepSolves = dep();
        forAll(newDepSolves, j)
        {
            if (newDepSolves[j].size() && newDepSolves[j][0] == '~')
            {
                // Removal
                word solveName = newDepSolves[j].substr(1);
                label k =
                    depSolves.find(solveID(solveName, regionName_));
                if (k >= 0)
                {
                    for (; k < depSolves.size()-1; k++)
                    {
                        depSolves[k] = depSolves[k+1];
                    }
                    depSolves.resize(depSolves.size()-1);
                }
            }
            else
            {
                solveID newDep(newDepSolves[j], regionName_);
                label k = depSolves.find(newDep);
                if (k < 0)
                {
                    depSolves.append(newDep);
                }
            }
        }
        requiredDependencies.set(key, depSolves);
    }

    // Add any correctors specified in user input
    HashTable<wordList> inputCorrs
    (
        coeffs_.lookupOrDefault("correctors", HashTable<wordList>())
    );
    forAllConstIters(inputCorrs, corr)
    {
        wordList corrSolves = correctorMembers.lookup(corr.key(), wordList());
        const wordList& newCorrSolves = corr();
        forAll(newCorrSolves, j)
        {
            if (newCorrSolves[j].size() && newCorrSolves[j][0] == '~')
            {
                // Removal
                word solveName = newCorrSolves[j].substr(1);
                label k = corrSolves.find(solveName);
                if (k >= 0)
                {
                    for (; k < corrSolves.size()-1; k++)
                    {
                        corrSolves[k] = corrSolves[k+1];
                    }
                    corrSolves.resize(corrSolves.size()-1);
                }
            }
            else
            {
                label k = corrSolves.find(newCorrSolves[j]);
                if (k < 0)
                {
                    corrSolves.append(newCorrSolves[j]);
                }
            }
        }
        correctorMembers.set(corr.key(), corrSolves);
    }
}


bool Foam::solverObject::isFinalCorrectorBase
(
    const label corrector,
    const word& correctorName
)
{
    correctorName_ = correctorName;
    correctorNumber_.set(correctorName, corrector);

    if (correctorName_ == solverObject::outerCorrectorName)
    {
        nCorrectorsSinceLastActive_++;
    }

    if
    (
        mesh_.time().timeOutputValue() >= timeStart_
     && mesh_.time().timeOutputValue() <= timeEnd_
    )
    {
        return this->isFinalCorrector(corrector, correctorName);
    }
    else
    {
        return true;
    }

}


void Foam::solverObject::beginIterationBase
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    finalIter_.set(correctorName, finalIter);

    if (active())
    {
        this->beginIteration(corrector, correctorName, finalIter);
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::solverObject::assembleScalarMatrixBase
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    if (active())
    {
        return this->assembleScalarMatrix(solveName, finalSolve, dictName);
    }
    else
    {
        return tmp<fvScalarMatrix>();
    }
}


Foam::tmp<Foam::fvVectorMatrix> Foam::solverObject::assembleVectorMatrixBase
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    if (active())
    {
        return this->assembleVectorMatrix(solveName, finalSolve, dictName);
    }
    else
    {
        return tmp<fvVectorMatrix>();
    }
}


Foam::tmp<Foam::fvSymmTensorMatrix>
Foam::solverObject::assembleSymmTensorMatrixBase
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    if (active())
    {
        return this->assembleSymmTensorMatrix(solveName, finalSolve, dictName);
    }
    else
    {
        return tmp<fvSymmTensorMatrix>();
    }
}


Foam::tmp<Foam::fvTensorMatrix> Foam::solverObject::assembleTensorMatrixBase
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    if (active())
    {
        return this->assembleTensorMatrix(solveName, finalSolve, dictName);
    }
    else
    {
        return tmp<fvTensorMatrix>();
    }
}


Foam::tmp<Foam::fvBlockMatrix<Foam::vector4>>
Foam::solverObject::assembleVector4MatrixBase
(
    const word& solveName,
    bool& finalSolve,
    word& dictName
)
{
    if (active())
    {
        return this->assembleVector4Matrix(solveName, finalSolve, dictName);
    }
    else
    {
        return tmp<fvBlockMatrix<vector4>>();
    }
}


void Foam::solverObject::correctBase
(
    const word& solveName,
    const word& regionName
)
{
    if (active())
    {
        this->correct(solveName, regionName);
    }
}


void Foam::solverObject::endIterationBase
(
    const label corrector,
    const word& correctorName,
    const bool finalIter
)
{
    if (active())
    {
        this->endIteration(corrector, correctorName, finalIter);
    }
}


void Foam::solverObject::solve()
{
    if (!legacySolver_)
    {
        return;
    }

    // If running from a legacy solver, the assembly and solving is not done
    // externally, so we do it here in a naive, partitioned way

    if (mesh().time().timeIndex() == lastTimeIndex_)
    {
        legacyOuterCorrector_++;
        if (oneExecPerIter_)
        {
            return;
        }
    }
    else
    {
        lastTimeIndex_ = mesh().time().timeIndex();
        legacyOuterCorrector_ = 0;
    }

    DynamicList<word> solveNames;
    HashTable<wordList> derivedFields;
    SolveTable<solveList> requiredDependencies;
    SolveTable<solveList> optionalDependencies;
    HashTable<wordList> correctorMembers;
    this->getSolveGraphBase
    (
        solveNames,
        derivedFields,
        requiredDependencies,
        optionalDependencies,
        correctorMembers
    );
    correctorMembers.set(timeLoopCorrectorName, solveNames);

    // Simulate the process that foamSolve would use so that the object
    // is kept in a consistent state; however, run through the solves in
    // order returned - no scheduling is done.
    // Iterate on solves at each corrector level
    DynamicList<word> activeCorrectors;
    HashTable<label> correctorIter;
    HashTable<bool> finalIter;
    word currentCorrector = timeLoopCorrectorName;
    correctorIter.set(timeLoopCorrectorName, lastTimeIndex_);
    do   // Repeated iterations for correctors
    {
        activeCorrectors.clear();
        forAll(solveNames, solvei)
        {
            // Look for starting correctors
            forAllConstIters(correctorMembers, corr)
            {
                if (corr().found(solveNames[solvei]))
                {
                    const word& correctorName = corr.key();
                    if (!activeCorrectors.found(correctorName))
                    {
                        activeCorrectors.append(correctorName);
                        if (correctorIter.lookup(correctorName, 0) == 0)
                        {
                            currentCorrector = correctorName;
                            if (currentCorrector == outerCorrectorName)
                            {
                                correctorIter.set
                                (
                                    currentCorrector,
                                    legacyOuterCorrector_
                                );
                            }
                            else
                            {
                                correctorIter.set(currentCorrector, 0);
                            }

                        }
                        finalIter.set
                        (
                            correctorName,
                            this->isFinalCorrectorBase
                            (
                                correctorIter[correctorName],
                                correctorName
                            )
                        );
                        this->beginIterationBase
                        (
                            correctorIter[correctorName],
                            correctorName,
                            finalIter[correctorName]
                        );
                    }
                }
            }
            if (currentCorrector == activeCorrectors.last())
            {
                // Whether the 'final' solvers get selected
                bool finalSolver =
                    mesh().data::template lookupOrDefault<bool>
                    (
                        "finalIteration",
                        false
                    );

                const word& solveName = solveNames[solvei];
                word dictName = word::null;
                tmp<fvScalarMatrix> sm =
                    this->assembleScalarMatrixBase
                    (
                        solveName, finalSolver, dictName
                    );
                tmp<fvVectorMatrix> vm =
                    this->assembleVectorMatrixBase
                    (
                        solveName, finalSolver, dictName
                    );
                tmp<fvSymmTensorMatrix> stm =
                    this->assembleSymmTensorMatrixBase
                    (
                        solveName, finalSolver, dictName
                    );
                tmp<fvTensorMatrix> tm =
                    this->assembleTensorMatrixBase
                    (
                        solveName, finalSolver, dictName
                    );

                tmp<fvBlockMatrix<vector4>> vm4 =
                    this->assembleVector4MatrixBase
                    (
                        solveName, finalSolver, dictName
                    );
                if (finalSolver)
                {
                    dictName += "Final";
                }
                if (sm.valid())
                {
                    sm.ref().solve
                    (
                        mesh().solution().solverDict
                        (
                            dictName == word::null
                          ? sm->psi().select(finalSolver)
                          : dictName
                        )
                    );
                }
                if (vm.valid())
                {
                    vm.ref().solve
                    (
                        mesh().solution().solverDict
                        (
                            dictName == word::null
                          ? vm->psi().select(finalSolver)
                          : dictName
                        )
                    );
                }
                if (stm.valid())
                {
                    stm.ref().solve
                    (
                        mesh().solution().solverDict
                        (
                            dictName == word::null
                          ? stm->psi().select(finalSolver)
                          : dictName
                        )
                    );
                }
                if (tm.valid())
                {
                    tm.ref().solve
                    (
                        mesh().solution().solverDict
                        (
                            dictName == word::null
                          ? tm->psi().select(finalSolver)
                          : dictName
                        )
                    );
                }
                if (vm4.valid())
                {
                    vm4.ref().solve();
                }
                this->correctBase(solveName, "");
            }

            // Check for ending iteration levels
            const word nextSolveName =
                solvei < solveNames.size()-1
              ? solveNames[solvei+1]
              : word::null;
            wordList correctorNames = correctorMembers.toc();
            forAllReverse(correctorNames, correctori)
            {
                const word& correctorName = correctorNames[correctori];
                if (activeCorrectors.found(correctorName))
                {
                    if (!correctorMembers[correctorName].found(nextSolveName))
                    {
                        const label pos = activeCorrectors.find(correctorName);
                        if (pos >= 0)
                        {
                            if (pos != activeCorrectors.size()-1)
                            {
                                FatalErrorInFunction
                                    << "Error running solverObject "
                                    << this->name() << " in legacy mode: "
                                    << "correctors not properly nested. For "
                                    << "legacy mode, they need to be listed in "
                                    << "order of nesting."
                                    << exit(FatalError);
                            }
                            this->endIterationBase
                            (
                                correctorIter[activeCorrectors[pos]],
                                activeCorrectors[pos],
                                finalIter[activeCorrectors[pos]]
                            );
                            if
                            (
                                currentCorrector == activeCorrectors[pos]
                             && (
                                    finalIter[currentCorrector]
                                 || currentCorrector == timeLoopCorrectorName
                                 || currentCorrector == outerCorrectorName
                                )
                            )
                            {
                                if (pos > 0)
                                {
                                    currentCorrector = activeCorrectors[pos-1];
                                }
                                else
                                {
                                    currentCorrector = word::null;
                                }
                            }
                            activeCorrectors.resize(pos);
                        }
                    }
                }
            }
        }
        if (currentCorrector != word::null)
        {
            correctorIter.set
            (
                currentCorrector,
                correctorIter.lookup(currentCorrector, 0)+1
            );
        }
    } while (currentCorrector != word::null);
}


// ************************************************************************* //
