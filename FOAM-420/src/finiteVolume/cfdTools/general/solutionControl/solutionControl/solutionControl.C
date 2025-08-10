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
    (c) 2017-2021 Esi Ltd.
    (c) 2017 OpenCFD Ltd

\*---------------------------------------------------------------------------*/

#include "solutionControl.H"
#include "VectorN/primitives/vector4/vector4.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/solutionControl/pisoControl/pisoControl.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solutionControl, 0);
    defineRunTimeSelectionTable(solutionControl, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::solutionControl::readControls(const bool absTolOnly)
{
    const dictionary& solutionDict = this->dict();

    // Read solution controls
    nNonOrthCorr_ =
        solutionDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);
    momentumPredictor_ =
        solutionDict.lookupOrDefault("momentumPredictor", true);
    transonic_ = solutionDict.lookupOrDefault("transonic", false);
    consistent_ = solutionDict.lookupOrDefault("consistent", false);
    modifiedMomentumInterp_ =
        solutionDict.lookupOrDefault("modifiedMomentumInterp", false);
    frozenFlow_ = solutionDict.lookupOrDefault("frozenFlow", false);

    // Read residual information
    const dictionary residualDict
    (
        solutionDict.subOrEmptyDict("residualControl")
    );

    DynamicList<fieldData> data(residualControl_);

    forAllConstIter(dictionary, residualDict, iter)
    {
        const word& fName = iter().keyword();
        const label fieldi = applyToField(fName, false);
        if (fieldi == -1)
        {
            fieldData fd;
            fd.name = fName.c_str();

            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(fName));
                fd.relTol = -1;
                fd.initialResidual = -1;
            }
            else
            {
                if (iter().isDict())
                {
                    const dictionary& fieldDict(iter().dict());
                    fd.absTol = readScalar(fieldDict.lookup("tolerance"));
                    fd.relTol = readScalar(fieldDict.lookup("relTol"));
                    fd.initialResidual = 0.0;
                }
                else
                {
                    FatalErrorInFunction
                        << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }

            data.append(fd);
        }
        else
        {
            fieldData& fd = data[fieldi];
            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(fName));
            }
            else
            {
                if (iter().isDict())
                {
                    const dictionary& fieldDict(iter().dict());
                    fd.absTol = readScalar(fieldDict.lookup("tolerance"));
                    fd.relTol = readScalar(fieldDict.lookup("relTol"));
                }
                else
                {
                    FatalErrorInFunction
                        << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }
        }
    }

    residualControl_.transfer(data);

    if (debug)
    {
        forAll(residualControl_, i)
        {
            const fieldData& fd = residualControl_[i];
            Info<< "residualControl[" << i << "]:" << nl
                << "    name     : " << fd.name << nl
                << "    absTol   : " << fd.absTol << nl
                << "    relTol   : " << fd.relTol << nl
                << "    iniResid : " << fd.initialResidual << endl;
        }
    }
}


void Foam::solutionControl::readControls()
{
    readControls(false);
}


Foam::label Foam::solutionControl::applyToField
(
    const word& fieldName,
    const bool useRegEx
) const
{
    forAll(residualControl_, i)
    {
        if (useRegEx && residualControl_[i].name.match(fieldName))
        {
            return i;
        }
        else if (residualControl_[i].name == fieldName)
        {
            return i;
        }
    }

    return -1;
}

bool Foam::solutionControl::storeVars() const
{
    return false;
}


void Foam::solutionControl::storePrevIterFields() const
{
//    storePrevIter<label>();
    storePrevIter<scalar>();
    storePrevIter<vector>();
    storePrevIter<sphericalTensor>();
    storePrevIter<symmTensor>();
    storePrevIter<tensor>();
}

void Foam::solutionControl::setAdjointOptimization()
{
    adjointOptimization_ = true;
}

Foam::scalar Foam::solutionControl::maxResidual
(
    const entry& solverPerfDictEntry,
    scalar& lastRes
) const
{
    return maxResidual(obr_, solverPerfDictEntry, lastRes);
}


void Foam::solutionControl::setFirstIterFlag
(
    const bool check,
    const bool force
)
{
    if (debug)
    {
        Info
            << "solutionControl: force:" << force
            << " check: " << check
            << " corr: " << corr_
            << " corrNonOrtho:" << corrNonOrtho_
            << endl;
    }

    if (force || (check && corr_ <= 1 && corrNonOrtho_ == 0))
    {
        if (debug)
        {
            Info<< "solutionControl: set firstIteration flag" << endl;
        }
        mesh_.data::set("firstIteration", true);
    }
    else
    {
        if (debug)
        {
            Info<< "solutionControl: remove firstIteration flag" << endl;
        }
        mesh_.data::remove("firstIteration");
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::solutionControl::maxTypeResidual
(
    const objectRegistry& obr,
    const entry& solverPerfDictEntry,
    scalar& firstRes,
    scalar& lastRes
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    const word& fieldName = solverPerfDictEntry.keyword();

    if (obr.foundObject<fieldType>(fieldName))
    {
        const List<SolverPerformance<Type>> sp(solverPerfDictEntry.stream());
        firstRes = cmptMax(sp.first().initialResidual());
        lastRes = cmptMax(sp.last().initialResidual());
    }
}

template<class Type>
void Foam::solutionControl::maxTypeBlockResidual
(
    const objectRegistry& obr,
    const entry& solverPerfDictEntry,
    scalar& firstRes,
    scalar& lastRes
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    const word& fieldName = solverPerfDictEntry.keyword();

    if (obr.foundObject<fieldType>(fieldName))
    {
        const List<BlockSolverPerformance<Type>> sp(solverPerfDictEntry.stream());
        firstRes = cmptMax(sp.first().initialResidual());
        lastRes = cmptMax(sp.last().initialResidual());
    }
}


Foam::scalar Foam::solutionControl::maxResidual
(
    const objectRegistry& obr,
    const entry& solverPerfDictEntry,
    scalar& lastRes
)
{
    scalar firstRes = 0;

    maxTypeResidual<scalar>(obr, solverPerfDictEntry, firstRes, lastRes);
    maxTypeResidual<vector>(obr, solverPerfDictEntry, firstRes, lastRes);
    maxTypeResidual<sphericalTensor>(obr, solverPerfDictEntry, firstRes, lastRes);
    maxTypeResidual<symmTensor>(obr, solverPerfDictEntry, firstRes, lastRes);
    maxTypeResidual<tensor>(obr, solverPerfDictEntry, firstRes, lastRes);
    maxTypeBlockResidual<vector4>(obr, solverPerfDictEntry, firstRes, lastRes);

    return firstRes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionControl::solutionControl
(
    fvMesh& mesh,
    const objectRegistry& obr,
    const word& algorithmName
)
:
    regIOobject
    (
        IOobject
        (
            "solutionControl",
            obr.time().timeName(),
            obr
        )
    ),
    mesh_(mesh),
    obr_(obr),
    residualControl_(),
    algorithmName_(algorithmName),
    nNonOrthCorr_(0),
    momentumPredictor_(true),
    transonic_(false),
    consistent_(false),
    modifiedMomentumInterp_(false),
    frozenFlow_(false),
    adjointOptimization_(false),
    corr_(0),
    corrNonOrtho_(0)
{}

Foam::solutionControl::solutionControl
(
    fvMesh& mesh,
    const word& algorithmName
)
:
    solutionControl(mesh, mesh, algorithmName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solutionControl> Foam::solutionControl::New
(
    fvMesh& mesh,
    const objectRegistry& obr
)
{
    const fvSolution& fvSoln = mesh;

    word solnControlName = word::null;
    const dictionary* dictPtr =
        fvSoln.subDictPtr(solutionControl::typeName);
    if (dictPtr)
    {
        solnControlName = word(dictPtr->lookup("type"));
    }

    // Try old-style subdicts if solutionControl is not specified
    // Look for non-empty dicts first time around; failing that empty
    for (label i = 0; i < 2; i++)
    {
        if (solnControlName == word::null)
        {
            // Check all registered options
            for (const auto& p : dictionaryConstructorTable_()) {
                const word& type = p.first;
                const entry* entryPtr =
                    fvSoln.lookupEntryPtr(type, false, true);

                if (entryPtr && entryPtr->isDict())
                {
                    if (!entryPtr->dict().empty() || i == 1)
                    {
                        solnControlName = type;
                    }
                }
            }
        }
    }

    if (solnControlName == word::null)
    {
        FatalIOErrorInFunction(fvSoln)
            << "Please specify solution control type in solutionControl "
            << nl
            << "dictionary, since none of the following were found in "
            << fvSoln.name() << ":" << nl
            << exit(FatalIOError);
    }

    Info<< "Selecting solution control: " << solnControlName << endl;

    const auto ctor = ctorTableLookup("solution control type", dictionaryConstructorTable_(), solnControlName);
    return autoPtr<solutionControl>(ctor(mesh, obr));
}

Foam::autoPtr<Foam::solutionControl> Foam::solutionControl::NewSIMPLE
(
    fvMesh& mesh,
    const objectRegistry& obr
)
{
    word solnControlName = "SIMPLE";

    Info<< "Selecting solution control: " << solnControlName << endl;

    const auto ctor = ctorTableLookup("solution control type", dictionaryConstructorTable_(), solnControlName);
    return autoPtr<solutionControl>(ctor(mesh, obr));
}

Foam::solutionControl& Foam::solutionControl::lookupOrCreate
(
    fvMesh& mesh, const objectRegistry& obr
)
{
    if (!obr.foundObject<solutionControl>(solutionControl::typeName))
    {
        if (solutionControl::debug)
        {
            Pout<< "solutionControl::lookupOrCreate(const fvMesh&, const objectRegistry&) : "
                << "constructing solution control "
                << " for region " << obr.name() << endl;
        }
        regIOobject::store(New(mesh, obr).ptr());
    }

    return obr.lookupObjectRef<solutionControl>(solutionControl::typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionControl::~solutionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solutionControl::writeData(Ostream&) const
{
    return true;
}

// ************************************************************************* //
