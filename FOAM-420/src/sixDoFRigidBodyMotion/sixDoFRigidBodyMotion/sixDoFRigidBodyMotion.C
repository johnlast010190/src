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
    (c) 2015 OpenCFD Ltd.
    (c) 2017-2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyMotion/sixDoFRigidBodyMotion.H"
#include "sixDoFSolvers/sixDoFSolver/sixDoFSolver.H"
#include "primitives/septernion/septernion.H"
#include "momentOfInertia/momentOfInertia.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::updateCoeffsBasedonDensity
(
    const dictionary& dict,
    const triSurfaceMesh& stl
)
{
    scalar rhoSolid = readScalar(dict.lookup("rhoSolid"));
    scalar mass = 0.0;
    vector centreOfMass = Zero;
    tensor momentOfInertia = Zero;

    momentOfInertia::massPropertiesSolid
    (
        stl,
        rhoSolid,
        mass,
        centreOfMass,
        momentOfInertia
    );

    vector eVal = eigenValues(momentOfInertia);
    tensor eVec = eigenVectors(momentOfInertia);

    label pertI = 0;

    Random rand(57373);

    //- Double coding. Code exists at surfaceInertia utility.
    //  encapsulation is needed
    //  toDo: extend it for 2D cases because now works for the whole stl

    while ((magSqr(eVal) < VSMALL) && pertI < 10)
    {
        WarningInFunction
            << "No eigenValues found, shape may have symmetry, "
            << "perturbing inertia tensor diagonal" << endl;

        momentOfInertia.xx() *= 1.0 + SMALL*rand.sample01<scalar>();
        momentOfInertia.yy() *= 1.0 + SMALL*rand.sample01<scalar>();
        momentOfInertia.zz() *= 1.0 + SMALL*rand.sample01<scalar>();

        eVal = eigenValues(momentOfInertia);

        eVec = eigenVectors(momentOfInertia);

        pertI++;
    }

    if
    (
        (mag(eVec.x() ^ eVec.y()) > (1.0 - SMALL))
     && (mag(eVec.y() ^ eVec.z()) > (1.0 - SMALL))
     && (mag(eVec.z() ^ eVec.x()) > (1.0 - SMALL))
    )
    {
        // Make the eigenvectors a right handed orthogonal triplet
        eVec = tensor
        (
            eVec.x(),
            eVec.y(),
            eVec.z() * sign((eVec.x() ^ eVec.y()) & eVec.z())
        );

        // Finding the most natural transformation.  Using Lists
        // rather than tensors to allow indexed permutation.

        // Cartesian basis vectors - right handed orthogonal triplet
        List<vector> cartesian(3);

        cartesian[0] = vector(1, 0, 0);
        cartesian[1] = vector(0, 1, 0);
        cartesian[2] = vector(0, 0, 1);

        // Principal axis basis vectors - right handed orthogonal
        // triplet
        List<vector> principal(3);

        principal[0] = eVec.x();
        principal[1] = eVec.y();
        principal[2] = eVec.z();

        scalar maxMagDotProduct = -GREAT;

        // Matching axis indices, first: cartesian, second:principal

        Pair<label> match(-1, -1);

        forAll(cartesian, cI)
        {
            forAll(principal, pI)
            {
                scalar magDotProduct = mag(cartesian[cI] & principal[pI]);

                if (magDotProduct > maxMagDotProduct)
                {
                    maxMagDotProduct = magDotProduct;

                    match.first() = cI;

                    match.second() = pI;
                }
            }
        }

        scalar sense = sign
        (
            cartesian[match.first()] & principal[match.second()]
        );

        if (sense < 0)
        {
            // Invert the best match direction and swap the order of
            // the other two vectors

            List<vector> tPrincipal = principal;

            tPrincipal[match.second()] *= -1;

            tPrincipal[(match.second() + 1) % 3] =
                principal[(match.second() + 2) % 3];

            tPrincipal[(match.second() + 2) % 3] =
                principal[(match.second() + 1) % 3];

            principal = tPrincipal;

            vector tEVal = eVal;

            tEVal[(match.second() + 1) % 3] = eVal[(match.second() + 2) % 3];

            tEVal[(match.second() + 2) % 3] = eVal[(match.second() + 1) % 3];

            eVal = tEVal;
        }

        label permutationDelta = match.second() - match.first();

        if (permutationDelta != 0)
        {
            // Add 3 to the permutationDelta to avoid negative indices

            permutationDelta += 3;

            List<vector> tPrincipal = principal;

            vector tEVal = eVal;

            for (label i = 0; i < 3; i++)
            {
                tPrincipal[i] = principal[(i + permutationDelta) % 3];

                tEVal[i] = eVal[(i + permutationDelta) % 3];
            }

            principal = tPrincipal;

            eVal = tEVal;
        }

        label matchedAlready = match.first();

        match =Pair<label>(-1, -1);

        maxMagDotProduct = -GREAT;

        forAll(cartesian, cI)
        {
            if (cI == matchedAlready)
            {
                continue;
            }

            forAll(principal, pI)
            {
                if (pI == matchedAlready)
                {
                    continue;
                }

                scalar magDotProduct = mag(cartesian[cI] & principal[pI]);

                if (magDotProduct > maxMagDotProduct)
                {
                    maxMagDotProduct = magDotProduct;

                    match.first() = cI;

                    match.second() = pI;
                }
            }
        }

        sense = sign
        (
            cartesian[match.first()] & principal[match.second()]
        );

        if (sense < 0 || (match.second() - match.first()) != 0)
        {
            principal[match.second()] *= -1;

            List<vector> tPrincipal = principal;

            tPrincipal[(matchedAlready + 1) % 3] =
                principal[(matchedAlready + 2) % 3]*-sense;

            tPrincipal[(matchedAlready + 2) % 3] =
                principal[(matchedAlready + 1) % 3]*-sense;

            principal = tPrincipal;

            vector tEVal = eVal;

            tEVal[(matchedAlready + 1) % 3] = eVal[(matchedAlready + 2) % 3];

            tEVal[(matchedAlready + 2) % 3] = eVal[(matchedAlready + 1) % 3];

            eVal = tEVal;
        }

        eVec = tensor(principal[0], principal[1], principal[2]);

    }
    else
    {
        WarningInFunction
            << "Non-unique eigenvectors, cannot compute transformation "
            << "from Cartesian axes" << endl;
    }

    initialCentreOfMass_ = centreOfMass;
    initialCentreOfRotation_ = initialCentreOfMass_;
    mass_ = mass;

    momentOfInertia_.xx() = eVal.x();
    momentOfInertia_.yy() = eVal.y();
    momentOfInertia_.zz() = eVal.z();

    initialQ_ = eVec.T();

    Info<< endl;
    Info<< "motion coeffs are computed via density and stl:" << endl;
    Info<< "stl surface information:" << endl;
    Info<< tab << "density: " << rhoSolid << endl;
    Info<< tab << "mass: " << mass << endl;
    Info<< tab << "momentOfInertia: " << momentOfInertia_ << endl;
    Info<< tab << "centreOfMass: " << initialCentreOfMass_ << endl;
    Info<< tab << "orientation: " << initialQ_ << endl;
}

void Foam::sixDoFRigidBodyMotion::initialise
(
    const dictionary& dict,
    const dictionary& stateDict
)
{
    addRestraints(dict);

    // Set constraints and initial centre of rotation
    // if different to the centre of mass
    addConstraints(dict);

    // If the centres of mass and rotation are different ...
    vector R(initialCentreOfMass_ - initialCentreOfRotation_);
    if (magSqr(R) > VSMALL)
    {
        // ... correct the moment of inertia tensor using parallel axes theorem
        momentOfInertia_ += mass_*diag(I*magSqr(R) - sqr(R));

        // ... and if the centre of rotation is not specified for motion state
        // update it
        if (!stateDict.found("centreOfRotation"))
        {
            motionState_.centreOfRotation() = initialCentreOfRotation_;
        }
    }

    // Save the old-time motion state
    motionState0_ = motionState_;
}

void Foam::sixDoFRigidBodyMotion::applyRestraints()
{
    if (restraints_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        forAll(restraints_, rI)
        {
            if (report_)
            {
                Info<< "Restraint " << restraints_[rI].name() << ": ";
            }

            // Restraint position
            point rP = Zero;

            // Restraint force
            vector rF = Zero;

            // Restraint moment
            vector rM = Zero;

            // Accumulate the restraints
            restraints_[rI].restrain(*this, rP, rF, rM);

            // Update the acceleration
            a() += rF/mass_;

            // Moments are returned in global axes, transforming to
            // body local to add to torque.
            tau() += Q().T() & (rM + ((rP - centreOfRotation()) ^ rF));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion()
:
    motionState_(),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_(Zero),
    initialCentreOfRotation_(Zero),
    initialQ_(I),
    mass_(VSMALL),
    momentOfInertia_(diagTensor::one*VSMALL),
    aRelax_(1.0),
    aDamp_(1.0),
    report_(false),
    solver_(nullptr)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const dictionary& dict,
    const dictionary& stateDict
)
:
    motionState_(stateDict),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_
    (
        dict.lookupOrDefault
        (
            "initialCentreOfMass",
            vector(dict.lookup("centreOfMass"))
        )
    ),
    initialCentreOfRotation_(initialCentreOfMass_),
    initialQ_
    (
        dict.lookupOrDefault
        (
            "initialOrientation",
            dict.lookupOrDefault("orientation", tensor::I)
        )
    ),
    mass_(readScalar(dict.lookup("mass"))),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    aRelax_(dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0)),
    aDamp_(dict.lookupOrDefault<scalar>("accelerationDamping", 1.0)),
    report_(dict.lookupOrDefault<Switch>("report", false)),
    solver_(sixDoFSolver::New(dict.subDict("solver"), *this))
{
    addRestraints(dict);
    initialise(dict, stateDict);
}

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const dictionary& dict,
    const dictionary& stateDict,
    const triSurfaceMesh& stl
)
:
    motionState_(stateDict),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_
    (
        dict.lookupOrDefault
        (
            "initialCentreOfMass",
            vector(dict.lookup("centreOfMass"))
        )
    ),
    initialCentreOfRotation_(initialCentreOfMass_),
    initialQ_
    (
        dict.lookupOrDefault
        (
            "initialOrientation",
            dict.lookupOrDefault("orientation", tensor::I)
        )
    ),
    mass_(readScalar(dict.lookup("mass"))),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    aRelax_(dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0)),
    aDamp_(dict.lookupOrDefault<scalar>("accelerationDamping", 1.0)),
    report_(dict.lookupOrDefault<Switch>("report", false)),
    solver_(sixDoFSolver::New(dict.subDict("solver"), *this))
{
    updateCoeffsBasedonDensity(dict, stl);
    initialise(dict, stateDict);
}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const sixDoFRigidBodyMotion& sDoFRBM
)
:
    motionState_(sDoFRBM.motionState_),
    motionState0_(sDoFRBM.motionState0_),
    restraints_(sDoFRBM.restraints_),
    constraints_(sDoFRBM.constraints_),
    tConstraints_(sDoFRBM.tConstraints_),
    rConstraints_(sDoFRBM.rConstraints_),
    initialCentreOfMass_(sDoFRBM.initialCentreOfMass_),
    initialCentreOfRotation_(sDoFRBM.initialCentreOfRotation_),
    initialQ_(sDoFRBM.initialQ_),
    mass_(sDoFRBM.mass_),
    momentOfInertia_(sDoFRBM.momentOfInertia_),
    aRelax_(sDoFRBM.aRelax_),
    aDamp_(sDoFRBM.aDamp_),
    report_(sDoFRBM.report_),
    solver_(sDoFRBM.solver_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::~sixDoFRigidBodyMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        forAllConstIter(IDLList<entry>, restraintDict, iter)
        {
            if (iter().isDict())
            {
                restraints_.set
                (
                    i++,
                    sixDoFRigidBodyMotionRestraint::New
                    (
                        iter().keyword(),
                        iter().dict()
                    )
                );
            }
        }

        restraints_.setSize(i);
    }
}


void Foam::sixDoFRigidBodyMotion::addConstraints
(
    const dictionary& dict
)
{
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        pointConstraint pct;
        pointConstraint pcr;

        forAllConstIter(IDLList<entry>, constraintDict, iter)
        {
            if (iter().isDict())
            {
                constraints_.set
                (
                    i,
                    sixDoFRigidBodyMotionConstraint::New
                    (
                        iter().keyword(),
                        iter().dict(),
                        *this
                    )
                );

                constraints_[i].setCentreOfRotation(initialCentreOfRotation_);
                constraints_[i].constrainTranslation(pct);
                constraints_[i].constrainRotation(pcr);

                i++;
            }
        }

        constraints_.setSize(i);

        tConstraints_ = pct.constraintTransformation();
        rConstraints_ = pcr.constraintTransformation();

        Info<< "Translational constraint tensor " << tConstraints_ << nl
            << "Rotational constraint tensor " << rConstraints_ << endl;
    }
}


void Foam::sixDoFRigidBodyMotion::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    static bool first = false;

    // Save the previous iteration accelerations for relaxation
    vector aPrevIter = a();
    vector tauPrevIter = tau();

    // Calculate new accelerations
    a() = fGlobal/mass_;
    tau() = (Q().T() & tauGlobal);
    applyRestraints();

    // Relax accelerations on all but first iteration
    if (!first)
    {
        a() = aRelax_*a() + (1 - aRelax_)*aPrevIter;
        tau() = aRelax_*tau() + (1 - aRelax_)*tauPrevIter;
    }

    first = false;
}


void Foam::sixDoFRigidBodyMotion::update
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    if (Pstream::master())
    {
        solver_->solve(firstIter, fGlobal, tauGlobal, deltaT, deltaT0);

        if (report_)
        {
            status();
        }
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::status() const
{
    vector angle = vector::zero;

    angle.x() = 180*atan(Q().zy()/Q().yy()) / constant::mathematical::pi;
    angle.y() = 180*atan(Q().xz()/Q().zz()) / constant::mathematical::pi;
    angle.z() = 180*atan(Q().yx()/Q().xx()) / constant::mathematical::pi;

    Info<< "6-DoF rigid body motion" << nl
        << "    Centre of rotation: " << centreOfRotation() << nl
        << "    Centre of mass: " << centreOfMass() << nl
        << "    Orientation: " << orientation() << nl
        << "    Angular displacement: " << angle << " [deg]" << nl
        << "    Linear velocity: " << v() << nl
        << "    Angular velocity: " << omega()
        << endl;
}


Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::transform
(
    const pointField& initialPoints
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoints - initialCentreOfRotation_))
    );
}


Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::transform
(
    const pointField& initialPoints,
    const scalarField& scale
) const
{
    // Calculate the transformation septerion from the initial state
    septernion s
    (
        centreOfRotation() - initialCentreOfRotation(),
        quaternion(Q().T() & initialQ())
    );

    tmp<pointField> tpoints(new pointField(initialPoints));
    pointField& points = tpoints.ref();

    forAll(points, pointi)
    {
        // Move non-stationary points
        if (scale[pointi] > SMALL)
        {
            // Use solid-body motion where scale = 1
            if (scale[pointi] > 1 - SMALL)
            {
                points[pointi] = transform(initialPoints[pointi]);
            }
            // Slerp septernion interpolation
            else
            {
                septernion ss(slerp(septernion::I, s, scale[pointi]));

                points[pointi] =
                    initialCentreOfRotation()
                  + ss.invTransformPoint
                    (
                        initialPoints[pointi]
                      - initialCentreOfRotation()
                    );
            }
        }
    }

    return tpoints;
}


// ************************************************************************* //
