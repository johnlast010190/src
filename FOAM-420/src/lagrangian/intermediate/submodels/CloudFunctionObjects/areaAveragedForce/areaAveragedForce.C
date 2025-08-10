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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/CloudFunctionObjects/areaAveragedForce/areaAveragedForce.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "sampledSurface/writers/surfaceWriter.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::areaAveragedForce<CloudType>::makeLogFile
(
    const word& zoneName,
    const label zoneI,
    const label nFaces,
    const scalar totArea
)
{
    // Create the output file if not already created
    if (log_)
    {
        if (debug)
        {
            Info<< "Creating output file." << endl;
        }

        if (Pstream::master())
        {
            // Create directory if does not exist
            mkDir(this->writeTimeDir());

            // Open new file at start up
            outputFilePtr_.set
            (
                zoneI,
                new OFstream
                (
                    this->writeTimeDir()/(type() + '_' + zoneName + ".dat")
                )
            );

            outputFilePtr_[zoneI]
                << "# Source    : " << type() << nl
                << "# Face zone : " << zoneName << nl
                << "# Faces     : " << nFaces << nl
                << "# Area      : " << totArea << nl
                << "# Time" << tab << "Pwt" << tab << "Pnt"<<endl;
        }
    }
}


template<class CloudType>
void Foam::areaAveragedForce<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();
    const faceZoneMesh& fzm = mesh.faceZones();

    forAll(faceZoneIDs_, zoneI)
    {
        scalar PwtP = 0;
        scalar PntP = 0;

        if (totalArea_[zoneI])
        {
            PwtP = mag(Fwt_[zoneI])/totalArea_[zoneI];
            PntP = mag(Fnt_[zoneI])/totalArea_[zoneI];
        }

        reduce(PwtP, sumOp<scalar>());
        reduce(PntP, sumOp<scalar>());

        Pwt_[zoneI] = PwtP;
        Pnt_[zoneI] = PntP;
    }

    const label procI = Pstream::myProcNo();

    Info<< type() << " output:" << nl;

    forAll(faceZoneIDs_, zoneI)
    {
        const word& zoneName = fzm[faceZoneIDs_[zoneI]].name();

        Info<< "    " << zoneName
            << ": Tangential Pressure = " << Pwt_[zoneI]
            << "; Normal Pressure = " << Pnt_[zoneI]
            << nl;

        if (outputFilePtr_.set(zoneI))
        {
            OFstream& os = outputFilePtr_[zoneI];
            os  << time.timeName() << token::TAB << Pwt_[zoneI] << token::TAB
                << Pnt_[zoneI] << endl;
        }
    }

    Info<< endl;


    if (surfaceFormat_ != "none")
    {
        forAll(faceZoneIDs_, zoneI)
        {
            const faceZone& fZone = fzm[faceZoneIDs_[zoneI]];

            labelList pointToGlobal;
            labelList uniqueMeshPointLabels;
            autoPtr<globalIndex> globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    fZone().meshPoints(),
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

            pointField uniquePoints(mesh.points(), uniqueMeshPointLabels);
            List<pointField> allProcPoints(Pstream::nProcs());
            allProcPoints[procI] = uniquePoints;
            Pstream::gatherList(allProcPoints);

            faceList faces(fZone().localFaces());
            forAll(faces, i)
            {
                inplaceRenumber(pointToGlobal, faces[i]);
            }
            List<faceList> allProcFaces(Pstream::nProcs());
            allProcFaces[procI] = faces;
            Pstream::gatherList(allProcFaces);

            if (Pstream::master())
            {
                pointField allPoints
                (
                    ListListOps::combine<pointField>
                    (
                        allProcPoints, accessOp<pointField>()
                    )
                );

                faceList allFaces
                (
                    ListListOps::combine<faceList>
                    (
                        allProcFaces, accessOp<faceList>()
                    )
                );

                autoPtr<surfaceWriter> writer
                (
                    surfaceWriter::New
                    (
                        surfaceFormat_,
                        this->coeffDict().subOrEmptyDict("formatOptions").
                            subOrEmptyDict(surfaceFormat_)
                    )
                );

                writer->write
                (
                    this->writeTimeDir(),
                    fZone.name(),
                    meshedSurfRef(allPoints, allFaces),
                    "tangentialForce",
                    Fw_[zoneI],
                    false
                );

                writer->write
                (
                    this->writeTimeDir(),
                    fZone.name(),
                    meshedSurfRef(allPoints, allFaces),
                    "normalForce",
                    Fn_[zoneI],
                    false
                );
            }
        }
    }

    forAll(Fnt_, zoneI)
    {
        Fnt_[zoneI] = vector::zero;
        Fwt_[zoneI] = vector::zero;
    }

    // writeProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::areaAveragedForce<CloudType>::areaAveragedForce
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    faceZoneIDs_(),
    surfaceFormat_(this->coeffDict().template lookupOrDefault<word>("surfaceFormat", "none")),
    Fw_(),
    Fwt_(),
    Pwt_(),
    Fn_(),
    Fnt_(),
    Pnt_(),
    totalArea_(),
    bulkDir_(this->coeffDict().template lookupOrDefault<vector>("bulkDirection",vector::zero)),
    log_(this->coeffDict().lookup("log")),
    outputFilePtr_(),
    timeOld_(owner.mesh().time().value())
{

    wordList faceZoneNames(this->coeffDict().lookup("faceZones"));

    Fw_.setSize(faceZoneNames.size());//added
    Fn_.setSize(faceZoneNames.size());//added
    Fwt_.setSize(faceZoneNames.size(), Foam::vector(0,0,0));//added
    Pwt_.setSize(faceZoneNames.size(), 0.0);//added
    Fnt_.setSize(faceZoneNames.size(), Foam::vector(0,0,0));//added
    Pnt_.setSize(faceZoneNames.size(), 0.0);//added
    totalArea_.setSize(faceZoneNames.size(), 0.0);//added

    outputFilePtr_.setSize(faceZoneNames.size());

    DynamicList<label> zoneIDs;
    const faceZoneMesh& fzm = owner.mesh().faceZones();
    const surfaceScalarField& magSf = owner.mesh().magSf();
    const polyBoundaryMesh& pbm = owner.mesh().boundaryMesh();
    label nZonesFound = 0;
    forAll(faceZoneNames, i)
    {
        const word& zoneName = faceZoneNames[i];
        label zoneI = fzm.findZoneID(zoneName);
        if (zoneI != -1)
        {
            zoneIDs.append(zoneI);
            const faceZone& fz = fzm[zoneI];

            Fw_[nZonesFound].setSize(fz.size(), vector(0,0,0));//added
            Fn_[nZonesFound].setSize(fz.size(), vector(0,0,0));//added

            label nFaces = returnReduce(fz.size(), sumOp<label>());
            Info<< "        " << zoneName << " faces: " << nFaces << nl;

            scalar totArea = 0.0;
            totalArea_[i] = 0;

            forAll(fz, j)
            {
                label faceI = fz[j];
                if (faceI < owner.mesh().nInternalFaces())
                {
                    totArea += magSf[fz[j]];
                }
                else
                {
                    label bFaceI = faceI - owner.mesh().nInternalFaces();
                    label patchI = pbm.patchID()[bFaceI];
                    const polyPatch& pp = pbm[patchI];

                    if
                    (
                        !magSf.boundaryField()[patchI].coupled()
                     || refCast<const coupledPolyPatch>(pp).owner()
                    )
                    {
                        label localFaceI = pp.whichFace(faceI);
                        totArea += magSf.boundaryField()[patchI][localFaceI];
                    }
                }
            }
            totArea = returnReduce(totArea, sumOp<scalar>());
            totalArea_[i]=totArea;

            makeLogFile(zoneName, nZonesFound, nFaces, totArea);
            nZonesFound++;
        }
        else
        {
            WarningInFunction
                << "Could not find faceZone " << zoneName << endl;
        }
    }

    faceZoneIDs_.transfer(zoneIDs);

    // readProperties(); AND initialise mass... fields
}


template<class CloudType>
Foam::areaAveragedForce<CloudType>::areaAveragedForce
(
    const areaAveragedForce<CloudType>& pff
)
:
    CloudFunctionObject<CloudType>(pff),
    faceZoneIDs_(pff.faceZoneIDs_),
    surfaceFormat_(pff.surfaceFormat_),
    Fw_(pff.Fw_),
    Fwt_(pff.Fwt_),
    Pwt_(pff.Pwt_),
    Fn_(pff.Fn_),
    Fnt_(pff.Fnt_),
    Pnt_(pff.Pnt_),
    totalArea_(pff.totalArea_),
    bulkDir_(pff.bulkDir_),
    log_(pff.log_),
    outputFilePtr_(),
    timeOld_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::areaAveragedForce<CloudType>::~areaAveragedForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::areaAveragedForce<CloudType>::postFace
(
    const parcelType& p,
    const label faceI,
    bool&
)
{
    if
    (
        this->owner().solution().output()
     || this->owner().solution().transient()
    )
    {
        const faceZoneMesh& fzm = this->owner().mesh().faceZones();

        forAll(faceZoneIDs_, i)
        {
            const faceZone& fz = fzm[faceZoneIDs_[i]];

            label faceId = -1;
            forAll(fz, j)
            {
                if (fz[j] == faceI)
                {
                    faceId = j;
                    break;
                }
            }

            if (faceId != -1)
            {

            //- determine normal direction and normal force
                vector faceNormal = this->owner().mesh().Sf()[fz[faceId]];
                //- Sf points outward but we define it as inward in analysis
                vector nHat = -faceNormal/mag(faceNormal);
                // store normal force
                Fn_[i][faceId] = nHat*(p.parcelForces() & nHat);

            //- determine tangential direction and tangential force
                //- initialize the tangential direction
                vector faceTangential = vector::zero;

                //- if a user defines a bulk direction us that for
                //  tangential direction determination
                if (mag(bulkDir_))
                {
                    faceTangential = faceNormal ^ bulkDir_;
                } else
                {
                    faceTangential = faceNormal ^ p.U();
                }

                //- if the bulk direction or velocity is not
                //  parallel to the face normal, calculate the tangential forces
                if (mag(faceTangential))
                {
                    vector tHat = faceTangential/mag(faceTangential);
                    Fw_[i][faceId] = tHat*(p.parcelForces() & tHat);

                } else
                {
                    Fw_[i][faceId] = vector::zero;
                }

            //- sum the normal and tangential forces in the faceZone
                Fnt_[i] += Fn_[i][faceId];
                Fwt_[i] += Fw_[i][faceId];
            }
        }
    }
}


// ************************************************************************* //
