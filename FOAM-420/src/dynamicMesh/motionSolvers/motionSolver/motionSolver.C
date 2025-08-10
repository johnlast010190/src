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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/motionSolverList/motionSolverList.H"
#include "twoDPointCorrector/twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolver::motionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    mesh_(mesh),
    coeffDict_(dict.optionalSubDict(type + "Coeffs"))
{
    initTransforms();
}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::clone() const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New
(
    const polyMesh& mesh,
    const dictionary& solverDict
)
{
    if (solverDict.found("solvers"))
    {
        return autoPtr<motionSolver>(new motionSolverList(mesh, solverDict));
    }
    else
    {
        const word solverName
        (
            solverDict.found("motionSolver")
          ? solverDict.lookup("motionSolver")
          : solverDict.lookup("solver")
        );

        Info<< "Selecting motion solver: " << solverName << endl;

        const_cast<Time&>(mesh.time()).libs().open
        (
            solverDict,
            "motionSolverLibs",
            dictionaryConstructorTable_()
        );

        const auto ctor = ctorTableLookup("solver type", dictionaryConstructorTable_(), solverName);
        return autoPtr<motionSolver>(ctor(mesh, solverDict));
    }
}


Foam::motionSolver::iNew::iNew(const polyMesh& mesh)
:
    mesh_(mesh)
{}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::iNew::operator()
(
    Istream& is
) const
{
    dictionaryEntry dict(dictionary::null, is);

    return motionSolver::New(mesh_, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolver::~motionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::motionSolver::initTransforms()
{
    isSolidBody_ = false;
    if
    (
        coeffDict().found("referenceFrame")
     || coeffDict().found(solidBodyMotionFunction::typeName)
    )
    {
        if (coeffDict().found(solidBodyMotionFunction::typeName))
        {
            SBMFs_.setSize(1);
            SBMFs_.set
            (
                0,
                solidBodyMotionFunction::New(coeffDict(), mesh().time())
            );
            isSolidBody_ = true;
        }
        else
        {
            refFrames_.setSize(1);
            refFrames_.set
            (
                0,
                coordinateFrame::lookupNew
                (
                    dynamic_cast<const fvMesh&>(mesh()),
                    coeffDict()
                )
            );
            if (!refFrames_[0].anyDynamic())
            {
                FatalErrorInFunction
                    << "coordinateFrame " << refFrames_[0].name()
                    << " is linked to dynamic mesh but it is not dynamic."
                    << abort(FatalError);
            }
        }
    }
    else
    {
        wordList names;

        if (coeffDict().found("referenceFrames"))
        {
            names = coeffDict().lookup<wordList>("referenceFrames");
            refFrames_.setSize(names.size());
        }
        else
        {
            forAllConstIter(dictionary, coeffDict(), iter)
            {
                const word& name = iter().keyword();
                if (iter().isDict() && name != word("solver"))
                {
                    if
                    (
                        coeffDict().subDict(name).found
                        (
                            solidBodyMotionFunction::typeName
                        )
                    )
                    {
                        names.append(name);
                    }
                }
            }
            SBMFs_.setSize(names.size());
            isSolidBody_ = true;
        }
        forAll(names, namei)
        {
            const word& name = names[namei];

            if (isSolidBody_)
            {
                SBMFs_.set
                (
                    namei,
                    solidBodyMotionFunction::New
                    (
                        coeffDict().subDict(name),
                        mesh().time()
                    )
                );
            }
            else
            {
                refFrames_.set
                (
                    namei,
                    &coordinateFrame::New
                    (
                        dynamic_cast<const fvMesh&>(mesh()),
                        name
                    )
                );
                if (!refFrames_[namei].anyDynamic())
                {
                    FatalErrorInFunction
                        << "coordinateFrame " << refFrames_[namei].name()
                        << " is linked to dynamic mesh but it is not dynamic."
                        << abort(FatalError);
                }
            }
        }
    }
}


bool Foam::motionSolver::isIncrementalMotion() const
{
    // Please note all motions have to be relative or absolute
    if (SBMFs_.size())
    {
        return SBMFs_[0].isIncrementalMotion();
    }
    else if (refFrames_.size())
    {
        return refFrames_[0].isIncrementalMotion();
    }
    return false;
}


Foam::septernion Foam::motionSolver::transformation
(
    const label transformationI,
    const wordList& nestedFrames
) const
{
    if (isSolidBody_)
    {
        return SBMFs_[transformationI].transformation();
    }

    // Necessary to support old solvers
    // (foamSolve has it's own updateStates)
    coordinateFrame::updateStates(mesh().thisDb());

    return refFrames_[transformationI].transformation();
}


Foam::tmp<Foam::pointField> Foam::motionSolver::newPoints()
{
    solve();
    return curPoints();
}


void Foam::motionSolver::twoDCorrectPoints(pointField& p) const
{
    twoDPointCorrector::New(mesh_).correctPoints(p);
}


void Foam::motionSolver::updateMesh(const mapPolyMesh& mpm)
{}


// ************************************************************************* //
