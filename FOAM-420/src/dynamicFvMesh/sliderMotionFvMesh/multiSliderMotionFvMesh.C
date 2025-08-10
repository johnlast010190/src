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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sliderMotionFvMesh/multiSliderMotionFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/Fields/transformField/transformField.H"
#include "meshes/polyMesh/zones/ZoneMesh/cellZoneMesh.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSliderMotionFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        multiSliderMotionFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSliderMotionFvMesh::multiSliderMotionFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_(dynamicMeshDict().subDict(typeName + "Coeffs")),
    displacement_(points().size(), vector::zero)
{
    zoneIDs_.setSize(dynamicMeshCoeffs_.size());
    pointIDs_.setSize(dynamicMeshCoeffs_.size());
    movingIDs_.setSize(dynamicMeshCoeffs_.size());
    motionProfile_.setSize(dynamicMeshCoeffs_.size());
    velocityProfile_.setSize(dynamicMeshCoeffs_.size());
    initialDisplacement_.setSize(dynamicMeshCoeffs_.size());
    bufferSize_.setSize(dynamicMeshCoeffs_.size());
    blendType_.setSize(dynamicMeshCoeffs_.size());

    label zoneI = 0;

    const Time& t = time();
    scalar time = t.timeOutputValue();

    forAllConstIter(dictionary, dynamicMeshCoeffs_, iter)
    {
        if (iter().isDict())
        {
            zoneIDs_[zoneI] = cellZones().findZoneID(iter().keyword());

            if (zoneIDs_[zoneI] == -1)
            {
                FatalIOErrorIn
                (
                    "multiSliderMotionFvMesh::"
                    "multiSliderMotionFvMesh(const IOobject&)",
                    dynamicMeshCoeffs_
                )   << "Cannot find cellZone named " << iter().keyword()
                    << ". Valid zones are " << cellZones().names()
                    << exit(FatalIOError);
            }

            const dictionary& subDict = iter().dict();

            //List of moving patch id's to define fixed mesh displacement region
            wordList patchNames = wordList(subDict.lookup("movingPatches"));
            HashSet<word> patchNameSet(patchNames);

            DynamicList<label> mIDs(patchNames.size());
            forAll(boundary(), patchI)
            {
                word name = boundary()[patchI].name();
                if (patchNameSet.found(name))
                {
                    mIDs.append(patchI);
                }
            }
            mIDs.shrink();
            movingIDs_[zoneI].transfer(mIDs);

            //Velocity DataEntry for each zone
            if (subDict.found("velocity"))
            {
                velocityProfile_[zoneI] = true;
                motionProfile_.set
                (
                    zoneI,
                    Function1<vector>::New("velocity", subDict)
                );
            }
            else if (subDict.found("displacement"))
            {
                velocityProfile_[zoneI] = false;
                motionProfile_.set
                (
                    zoneI,
                    Function1<vector>::New("displacement", subDict)
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Need displacement or velocity profile for zone : "
                    << iter().keyword()
                    << exit(FatalError);
            }

            if (time < SMALL)
            {
                //Define initial mesh displacement if present
                initialDisplacement_[zoneI] = subDict.lookupOrDefault<vector>
                (
                    "initialDisplacement",
                    vector::zero
                 );
            }
            else
            {
                initialDisplacement_[zoneI] = vector::zero;
            }

            //Define buffer region size and blending function
            bufferSize_[zoneI] =
                subDict.lookupOrDefault<scalar>("bufferSize",0.01);
            blendType_[zoneI] =
                subDict.lookupOrDefault<word>("blendType","linear");

            // Collect points of cell zone.
            const cellZone& cz = cellZones()[zoneIDs_[zoneI]];

            boolList movePts(nPoints(), false);

            forAll(cz, i)
            {
                label cellI = cz[i];
                const cell& c = cells()[cellI];
                forAll(c, j)
                {
                    const face& f = faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointI = f[k];
                        movePts[pointI] = true;
                    }
                }
            }

            syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zoneI].transfer(ptIDs);

            Info<< "Applying solid body motion to "
                << pointIDs_[zoneI].size() << " points of cellZone "
                << iter().keyword() << endl;

            zoneI++;
        }
    }
    zoneIDs_.setSize(zoneI);
    pointIDs_.setSize(zoneI);
    movingIDs_.setSize(zoneI);
    motionProfile_.setSize(zoneI);
    velocityProfile_.setSize(zoneI);
    initialDisplacement_.setSize(zoneI);
    bufferSize_.setSize(zoneI);
    blendType_.setSize(zoneI);

    //move points according if initial displacement defined

    slidePoints(initialDisplacement_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSliderMotionFvMesh::~multiSliderMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::multiSliderMotionFvMesh::makePatch
(
    const labelList& patchIDs
)
{
    const polyBoundaryMesh& patches = boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces.
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFaceI++;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(faces(), addressing),
            points()
        )
    );
}

void Foam::multiSliderMotionFvMesh::slidePoints
(
    const vectorField& zoneDisplacements
)
{
    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];
        if (mag(zoneDisplacements[i]) > SMALL)
        {
            vector dir = zoneDisplacements[i] /
                mag(zoneDisplacements[i]);
            const labelList& movingPatchIds = movingIDs_[i];
            autoPtr<indirectPrimitivePatch> pp(makePatch(movingPatchIds));

            scalar minRegionDP = GREAT;
            scalar maxRegionDP = -GREAT;

            forAll(pp().localPoints(), pointI)
            {
                scalar dProd = (pp().localPoints()[pointI] & dir);
                minRegionDP =  min(minRegionDP, dProd);
                maxRegionDP =  max(maxRegionDP, dProd);
            }

            scalar minZoneDP = GREAT;
            scalar maxZoneDP = -GREAT;

            pointField zPts(points(), zonePoints);
            scalarField zDP(zPts.size(), 0.);

            forAll(zPts, pointI)
            {
                zDP[pointI] = (zPts[pointI] & dir);
                minZoneDP =  min(minZoneDP, zDP[pointI]);
                maxZoneDP =  max(maxZoneDP, zDP[pointI]);
            }

            Foam::reduce(
                std::tie(minRegionDP, maxRegionDP, minZoneDP, maxZoneDP),
                ParallelOp<minOp<scalar>, maxOp<scalar>, minOp<scalar>, maxOp<scalar>>{},
                comm()
            );

            scalarField dispRatio(zPts.size(), 0.0);

            scalar totalLength = maxZoneDP - minZoneDP;

            scalar bufferLength = bufferSize_[i] * totalLength;

            if
            (
                bufferLength > 0.25*(maxZoneDP-maxRegionDP)
                || bufferLength > 0.25*(minRegionDP-minZoneDP)
            )
            {
                bufferLength = 0.0;
            }

            forAll(zPts, pointI)
            {
                scalar dn = zDP[pointI];

                if
                (
                    dn > minZoneDP + bufferLength
                    && dn < maxZoneDP - bufferLength
                )
                {
                    if
                    (
                        dn >= minRegionDP - bufferLength
                        && dn <= maxRegionDP + bufferLength
                    )
                    {
                        dispRatio[pointI] = 1.0;
                    }
                    else
                    {
                        if (dn > maxRegionDP + bufferLength)
                        {
                            scalar num = dn - (maxRegionDP + bufferLength);
                            scalar den = maxZoneDP - maxRegionDP
                                - 2*bufferLength;
                            dispRatio[pointI] = 1.0 - num/den;
                        }
                        else
                        {
                            scalar num = dn - (minZoneDP + bufferLength);
                            scalar den = minRegionDP - minZoneDP
                                - 2*bufferLength;
                            dispRatio[pointI] = num/den;
                        }
                    }
                }
            }

            if (blendType_[i] == "linear")
            {
                UIndirectList<point>(displacement_, zonePoints) =
                    dispRatio*zoneDisplacements[i];
            }
            else if (blendType_[i] == "cos")
            {
                UIndirectList<point>(displacement_, zonePoints) = 0.5*
                    (1.0-Foam::cos(dispRatio*Foam::constant::mathematical::pi))
                    *zoneDisplacements[i];
            }
            else if (blendType_[i] == "smoothStep")
            {
                UIndirectList<point>(displacement_, zonePoints) =
                    (3.*pow(dispRatio,2)-2.*pow(dispRatio,3))
                    *zoneDisplacements[i];
            }
            else
            {
                FatalErrorInFunction
                    <<"Blending scheme "<<blendType_[i]<<" not available"
                    << exit(FatalError);
            }
        }
        else
        {
            UIndirectList<point>(displacement_, zonePoints) = vector::zero;
        }
    }

    syncTools::syncPointList
    (
        *this,
        displacement_,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value
    );

    fvMesh::movePoints(points()+displacement_);
}


bool Foam::multiSliderMotionFvMesh::update()
{
    static bool hasWarned = false;

    const Time& t = time();
    vectorField zoneDisplacements(zoneIDs_.size());
    forAll(zoneIDs_, i)
    {
        scalar t0 = t.timeOutputValue() - t.deltaTValue();
        if (velocityProfile_[i])
        {
            vector vel = motionProfile_[i].value(t0);
            zoneDisplacements[i] = vel * t.deltaTValue();
        }
        else
        {
            scalar t1 = t.timeOutputValue();
            vector disp0 = motionProfile_[i].value(t0);
            vector disp1 = motionProfile_[i].value(t1);
            zoneDisplacements[i] = disp1 - disp0;
        }
    }

    slidePoints(zoneDisplacements);

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
