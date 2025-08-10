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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyGIBMotionSolver/sixDoFRigidBodyGIBMotionSolver.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "forces/forces.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyGIBMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        sixDoFRigidBodyGIBMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyGIBMotionSolver::sixDoFRigidBodyGIBMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    motionPtr_(),
    patches_(wordReList(coeffDict().lookup("patches"))),
    test_(coeffDict().lookupOrDefault<Switch>("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().lookupOrDefault<word>("rho", "rho")),
    curTimeIndex_(-1),
    triName_(coeffDict().lookup("triSurfaceName")),
    initStlPointsPtr_(nullptr),
    stlPointsDisPtr_(nullptr)
{
    if (rhoName_ == "rhoInf")
    {
        coeffDict().lookup("rhoInf") >> rhoInf_;
    }
    triSurfaceMesh stl
    (
        IOobject
        (
            triName_,
            mesh.time().constant(),
            "triSurface",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    if
    (
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
    )
    {
        motionPtr_.reset
        (
            new sixDoFRigidBodyMotion
            (
                coeffDict(),
                IOdictionary
                (
                    IOobject
                    (
                        "sixDoFRigidBodyMotionState",
                        mesh.time().timeName(),
                        "uniform",
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            )
        );
    }
    else if (coeffDict().found("rhoSolid"))
    {
        motionPtr_.reset
        (
            new sixDoFRigidBodyMotion
            (
                coeffDict(),
                coeffDict(),
                stl
            )
        );
    }
    else
    {
        motionPtr_.reset
        (
            new sixDoFRigidBodyMotion
            (
                coeffDict(),
                coeffDict()
            )
        );
    }
    initStlPointsPtr_ = new pointField(stl.points());
    const pointField& initStlPoints = *initStlPointsPtr_;
    /*
    stlPointsDisPtr_ = new pointField
    (
        initStlPointsPtr_->size(), vector::zero
    );
    */
    stlPointsDisPtr_ = new pointField
    (
        motionPtr_->transform(initStlPoints) - initStlPoints
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyGIBMotionSolver::~sixDoFRigidBodyGIBMotionSolver()
{
    delete initStlPointsPtr_;
    delete stlPointsDisPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyGIBMotionSolver::movePoints(const pointField& p)
{
    // No local data that needs adapting.
}

Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyGIBMotionSolver::curPoints() const
{
    const pointField& initStlPoints = *initStlPointsPtr_;
    const pointField& stlPointsDis = *stlPointsDisPtr_;
    return initStlPoints + stlPointsDis;
}


void Foam::sixDoFRigidBodyGIBMotionSolver::solve()
{
    const Time& t = mesh().time();

    // Store the motion state at the beginning of the time-stepbool
    bool firstIter = false;
    if (curTimeIndex_ != t.timeIndex())
    {
        motionPtr_->newTime();
        curTimeIndex_ = t.timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = mesh().lookupObject<uniformDimensionedVectorField>("g");
    }
    else if (coeffDict().found("g"))
    {
        coeffDict().lookup("g") >> g;
    }

    // scalar ramp = min(max((t.value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    if (test_)
    {
        motionPtr_->update
        (
            firstIter,
            ramp*(motionPtr_->mass()*g.value()),
            ramp*(motionPtr_->mass()*(motionPtr_->momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    else
    {
        dictionary forcesDict;
        forcesDict.add("type", functionObjects::forces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("rho", rhoName_);
        forcesDict.add("CofR", motionPtr_->centreOfRotation());

        functionObjects::forces f("forces", t, forcesDict);

        f.calcForcesMoment();

        motionPtr_->update
        (
            firstIter,
            ramp*(f.forceEff() + motionPtr_->mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + motionPtr_->mass()*(motionPtr_->momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    // Update the displacements
    pointField& initStlPoints = *initStlPointsPtr_;
    pointField& stlPointsDis_ = *stlPointsDisPtr_;

    stlPointsDis_ =
        motionPtr_->transform(initStlPoints) - initStlPoints;
}


bool Foam::sixDoFRigidBodyGIBMotionSolver::write() const
{
    IOdictionary dict
    (
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    motionPtr_->state().write(dict);

    return dict.regIOobject::write() && motionSolver::write();
}


// ************************************************************************* //
