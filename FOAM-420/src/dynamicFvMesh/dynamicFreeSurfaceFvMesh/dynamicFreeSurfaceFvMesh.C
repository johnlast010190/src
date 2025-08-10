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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicFreeSurfaceFvMesh/dynamicFreeSurfaceFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFreeSurfaceFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicFreeSurfaceFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicFreeSurfaceFvMesh::dynamicFreeSurfaceFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_(dynamicMeshDict().optionalSubDict(typeName + "Coeffs")),
    velocityProfile_(Function1<scalar>::New("velocity", dynamicMeshCoeffs_)),
    origin_(dynamicMeshCoeffs_.lookup("origin")),
    motionDir_(dynamicMeshCoeffs_.lookup("direction")),
    heightRef_(readScalar(dynamicMeshCoeffs_.lookup("heightRef"))),
    heightMax_(0.0),
    initialDisplacement_(readScalar(dynamicMeshCoeffs_.lookup("initialDisplacement"))),
    stationaryPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
     ),
    identifier_(stationaryPoints_.size(), 0.0),
    ratio_(stationaryPoints_.size(), 0.0)

{
    // normalise axis
    if (mag(motionDir_) < SMALL)
    {
        FatalErrorInFunction
            << "Badly motion direction: zero magnitude: " << motionDir_
            << abort(FatalError);
    }
    motionDir_ /= mag(motionDir_);

    scalar maxHeight = -GREAT;
    forAll(stationaryPoints_, pI)
    {
        scalar parallel = ((stationaryPoints_[pI] - origin_) & motionDir_);
        maxHeight = max(maxHeight, parallel);
    }
    reduce(maxHeight, maxOp<scalar>());

    heightMax_=maxHeight;

    // mark points to be moved (above heightRef_)
    forAll(identifier_, pI)
    {
        identifier_[pI] = pos0
            (((stationaryPoints_[pI] - origin_) & motionDir_)-(heightRef_+SMALL));
    }

    forAll(ratio_, pI)
    {
        if (identifier_[pI]  == 1.0)
            ratio_[pI] =
                (((stationaryPoints_[pI] - origin_) & motionDir_) - (heightRef_+SMALL))
                /(heightMax_-heightRef_);
    }

    stationaryPoints_.rename("stationaryPoints");

    Info<< "Performing a dynamic mesh calculation: " << endl
        << " origin: " << origin_
        << " direction: " << motionDir_
        << " heightRef: " << heightRef_
        << " heightMax: " << heightMax_<<endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicFreeSurfaceFvMesh::~dynamicFreeSurfaceFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicFreeSurfaceFvMesh::update()
{
    static bool hasWarned = false;

    pointField newPoints = this->points();

    // apply initial displacement
    if (this->time().timeIndex() == 1)
    {
        if (initialDisplacement_ != 0.0)
        {
            scalar dL =
                (heightMax_ - heightRef_ - initialDisplacement_);

            forAll(identifier_, pI)
            {
                if (identifier_[pI] == 1.0)
                {
                    scalar d = motionDir_& (newPoints[pI] - origin_);
                    scalar d0 = heightRef_;
                    scalar dMax = heightMax_;

                    scalar frac = (d-d0)/(dMax-d0);

                    newPoints[pI] -= dL*frac*motionDir_;
                }
            }
        }
    }

    scalar Un = velocityProfile_->value(time().value());

    forAll(identifier_, pI)
    {
        if (identifier_[pI] == 1.0)
        {
            // move points along motionDir_
            // using ratio to achive uniform spacing
            newPoints[pI] += ratio_[pI]*Un*time().deltaT().value()*motionDir_;
        }
    }

    fvMesh::movePoints(newPoints);

    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //
