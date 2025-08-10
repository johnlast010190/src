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
    (c) 2019 Esi Ltd.
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fieldCalcs/surfaceFieldCalc/surfaceFieldCalc.H"
#include "fvMesh/fvMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.H"
#include "sampledSurface/sampledSurface/sampledSurface.H"
#include "meshes/meshTools/mergePoints.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldCalcs
{
    defineTypeNameAndDebug(surfaceFieldCalc, 0);
    addToRunTimeSelectionTable(fieldCalc, surfaceFieldCalc, dictionary);
    addToRunTimeSelectionTable(functionObject, surfaceFieldCalc, dictionary);
}
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldCalcs::surfaceFieldCalc::regionTypes,
    4
>::names[] =
{
    "faceZone",
    "patch",
    "surface",
    "sampledSurface"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldCalcs::surfaceFieldCalc::regionTypes,
    4
> Foam::functionObjects::fieldCalcs::surfaceFieldCalc::regionTypeNames_;


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldCalcs::surfaceFieldCalc::calcTypes,
    4
>::names[] =
{
    "field",
    "swirlAngle",
    "bulkSwirlAngle",
    "axialDeviation"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldCalcs::surfaceFieldCalc::calcTypes,
    4
> Foam::functionObjects::fieldCalcs::surfaceFieldCalc::calcTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::obr() const
{
    if (regionType_ == stSurface)
    {
        return mesh_.lookupObject<objectRegistry>(regionName_);
    }
    else
    {
        return mesh_;
    }
}



void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::setProbeIndex
(
    const pointField& Cf
)
{
    probeIndex_.setSize(probeLocations_.size(), -1);

    List<List<scalar>> faceDist(Pstream::nProcs());

    faceDist[Pstream::myProcNo()].setSize(probeLocations_.size());
    faceDist[Pstream::myProcNo()] = GREAT;

    List<scalar>& myFaceDist = faceDist[Pstream::myProcNo()];

    forAll(probeLocations_, probei)
    {
        const vector& location = probeLocations_[probei];
        scalar minDist = GREAT;
        label minFace = -1;
        forAll(Cf, facei)
        {
            scalar dist = mag(location-Cf[facei]);
            if (dist < minDist)
            {
                minDist = dist;
                minFace = facei;
            }
        }
        probeIndex_[probei] = minFace;
        if (probeIndex_[probei] != -1)
        {
            myFaceDist[probei] = minDist;
        }
    }

    Pstream::allGatherList(faceDist);

    forAll(probeLocations_, probei)
    {
        label procWithClosestFace = -1;
        scalar closestFaceDist = GREAT;

        forAll(faceDist, proci)
        {
            if (faceDist[proci][probei] < closestFaceDist)
            {
                closestFaceDist = faceDist[proci][probei];
                procWithClosestFace = proci;
            }
            else if
            (
                faceDist[proci][probei] == closestFaceDist
                && proci < procWithClosestFace
            )
            {
                closestFaceDist = faceDist[proci][probei];
                procWithClosestFace = proci;
            }
        }

        if (Pstream::myProcNo() !=  procWithClosestFace)
        {
            probeIndex_[probei] = -1;
        }
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::setFaceZoneFaces()
{
    const label zoneId = mesh_.faceZones().findZoneID(regionName_);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown face zone name: " << regionName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneId];

    DynamicList<label> faceIds(fZone.size());
    DynamicList<label> facePatchIds(fZone.size());
    DynamicList<bool> faceFlip(fZone.size());

    forAll(fZone, i)
    {
        const label facei = fZone[i];

        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(facei);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = facei - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            faceIds.append(faceId);
            facePatchIds.append(facePatchId);
            faceFlip.append(fZone.flipMap()[i] ? true : false);
        }
    }

    faceId_.transfer(faceIds);
    facePatchId_.transfer(facePatchIds);
    faceFlip_.transfer(faceFlip);


    bool geomFlip(dict_.lookupOrDefault<Switch>("faceZoneFlip",false));

    if (geomFlip)
    {
        vector vec0 = vector::zero;
        forAll(faceId_, i)
        {
            label facei = faceId_[i];
            label patchi = facePatchId_[i];
            if (patchi == -1)
            {
                vec0 = mesh_.Sf()[facei];
                break;
            }
            else
            {
                vec0 = mesh_.Sf().boundaryField()[patchi][facei];
                break;
            }
        }
        reduce(vec0, maxMagSqrOp<vector>());

        forAll(faceId_, i)
        {
            label facei = faceId_[i];
            label patchi = facePatchId_[i];
            vector fv = vector::zero;
            if (patchi == -1)
            {
                fv = mesh_.Sf()[facei];
            }
            else
            {
                fv = mesh_.Sf().boundaryField()[patchi][facei];
            }
            faceFlip_[i] = (fv&vec0) < 0 ? true : false;
        }
    }

    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    if (debug)
    {
        Pout<< "Original face zone size = " << fZone.size()
            << ", new size = " << faceId_.size() << endl;
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::setPatchFaces()
{
    const label patchid = mesh_.boundaryMesh().findPatchID(regionName_);

    if (patchid < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown patch name: " << regionName_
            << ". Valid patch names are: "
            << mesh_.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const polyPatch& pp = mesh_.boundaryMesh()[patchid];

    label nFaces = pp.size();
    if (isA<emptyPolyPatch>(pp))
    {
        nFaces = 0;
    }

    faceId_.setSize(nFaces);
    facePatchId_.setSize(nFaces);
    faceFlip_.setSize(nFaces);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    forAll(faceId_, facei)
    {
        faceId_[facei] = facei;
        facePatchId_[facei] = patchid;
        faceFlip_[facei] = false;
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::combineMeshGeometry
(
    faceList& faces,
    pointField& points
) const
{
    List<faceList> allFaces(Pstream::nProcs());
    List<pointField> allPoints(Pstream::nProcs());

    labelList globalFacesIs(faceId_);
    forAll(globalFacesIs, i)
    {
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            globalFacesIs[i] += mesh_.boundaryMesh()[patchi].start();
        }
    }

    // Add local faces and points to the all* lists
    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh_.faces(), globalFacesIs),
        mesh_.points()
    );
    allFaces[Pstream::myProcNo()] = pp.localFaces();
    allPoints[Pstream::myProcNo()] = pp.localPoints();

    Pstream::gatherList(allFaces);
    Pstream::gatherList(allPoints);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, proci)
    {
        nFaces += allFaces[proci].size();
        nPoints += allPoints[proci].size();
    }

    faces.setSize(nFaces);
    points.setSize(nPoints);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const faceList& fcs = allFaces[Pstream::myProcNo()];
        forAll(fcs, i)
        {
            const face& f = fcs[i];
            face& newF = faces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        forAll(pts, i)
        {
            points[nPoints++] = pts[i];
        }
    }

    // Other proc data follows
    forAll(allFaces, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const faceList& fcs = allFaces[proci];
            forAll(fcs, i)
            {
                const face& f = fcs[i];
                face& newF = faces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[proci];
            forAll(pts, i)
            {
                points[nPoints++] = pts[i];
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newPoints;
    bool hasMerged = mergePoints
    (
        points,
        SMALL,
        false,
        oldToNew,
        newPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << points.size()
                << " down to " << newPoints.size() << " points" << endl;
        }

        points.transfer(newPoints);
        forAll(faces, i)
        {
            inplaceRenumber(oldToNew, faces[i]);
        }
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::
combineSurfaceGeometry
(
    faceList& faces,
    pointField& points
) const
{
    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());

        if (Pstream::parRun())
        {
            // Dimension as fraction of surface
            const scalar mergeDim = 1e-10*boundBox(s.points(), true).mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
    else if (surfacePtr_.valid())
    {
        const sampledSurface& s = surfacePtr_();

        if (Pstream::parRun())
        {
            // Dimension as fraction of mesh bounding box
            const scalar mergeDim = 1e-10*mesh_.bounds().mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
}


Foam::scalar
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::totalArea() const
{
    scalar totalArea;

    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());

        totalArea = gSum(s.magSf());
    }
    else if (surfacePtr_.valid())
    {
        totalArea = gSum(surfacePtr_().magSf());
    }
    else
    {
        totalArea = gSum(filterField(mesh_.magSf()));
    }

    return totalArea;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fieldCalcs::surfaceFieldCalc::needsGeom() const
{
    // Many operations use the Sf or Cf field
    return true;
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::initialise
(
    const dictionary& dict
)
{
    dict.lookup("name") >> regionName_;

    switch (regionType_)
    {
        case stFaceZone:
        {
            setFaceZoneFaces();
            surfacePtr_.clear();
            break;
        }
        case stPatch:
        {
            setPatchFaces();
            surfacePtr_.clear();
            break;
        }
        case stSurface:
        {
            const surfMesh& s = dynamicCast<const surfMesh>(obr());
            nFaces_ = returnReduce(s.size(), sumOp<label>());

            faceId_.clear();
            facePatchId_.clear();
            faceFlip_.clear();
            surfacePtr_.clear();
            break;
        }
        case stSampledSurface:
        {
            faceId_.clear();
            facePatchId_.clear();
            faceFlip_.clear();

            surfacePtr_ = sampledSurface::New
            (
                name(),
                mesh_,
                dict.subDict("sampledSurfaceDict")
            );
            surfacePtr_().update();
            nFaces_ =
                returnReduce(surfacePtr_().faces().size(), sumOp<label>());
            break;
        }
        default:
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << int(regionType_) << "(" << regionName_ << "):"
                << nl << "    Unknown region type. Valid region types are:"
                << regionTypeNames_ << nl
                << exit(FatalError);
        }
    }

    if (nFaces_ == 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Region has no faces" << exit(FatalError);
    }

    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    vectorField Cf;
    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());
        Cf = s.Cf();
    }
    else if (surfacePtr_.valid())
    {
        Cf = surfacePtr_().Cf();
    }
    else
    {
        Cf = filterField(mesh_.Cf());
    }

    if (probeLocations_.size() > 0)
    {
        setProbeIndex(Cf);
    }

    totalArea_ = totalArea();

    operations_.clear();
    operations_ = PtrList<dictionary>(dict.lookup("operations"));

    operationsWriter_.clear();
    operationsWriter_.setSize(operations_.size());
    forAll(operations_, opI)
    {
        const dictionary& opDict(operations_[opI]);
        word fieldName(opDict.lookup("fieldName"));
        word operation(opDict.lookup("operation"));

        word outputName(operation + "_" + fieldName);
        operationsWriter_.set
        (
            opI,
            new writeFile(obr_, name(), outputName, dict)
        );
    }

    probesWriter_.setSize(operations_.size());
    forAll(operations_, opI)
    {
        if (probeLocations_.size() > 0)
        {
            const dictionary& opDict(operations_[opI]);
            word fieldName = word(opDict.lookup("fieldName"))+word("_probes");
            word operation(opDict.lookup("operation"));

            word outputName(operation + "_" + fieldName);

            probesWriter_[opI] = autoPtr<writeFile>
            (
                new writeFile(obr_, name(), fieldName, dict)
            );
        }
        else
        {
            probesWriter_[opI] = autoPtr<writeFile>();
        }
    }

    Info<< "    total faces   = " << nFaces_ << nl
        << "    total area    = " << totalArea_ << nl;

    surfaceWriterPtr_.clear();
    const word surfaceFormat(dict.lookup("surfaceFormat"));

    if (surfaceFormat != "none")
    {
        surfaceWriterPtr_.reset
        (
            surfaceWriter::New
            (
                surfaceFormat,
                dict.subOrEmptyDict("formatOptions").
                subOrEmptyDict(surfaceFormat)
            ).ptr()
         );

        Info<< "    surfaceFormat = " << surfaceFormat << nl;
    }

    Info<< nl << endl;
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::writeFileHeader
(
    writeFile& writer
) const
{
    writer.writeCommented(writer.file(), "Region type : ");
    writer.file() << regionTypeNames_[regionType_] << " "
                  << regionName_ << endl;

    writer.writeHeaderValue(writer.file(), "Faces", nFaces_);
    writer.writeHeaderValue(writer.file(), "Area", totalArea_);

    writer.writeCommented(writer.file(), "Time");
    writer.file()  << "Min, Max, Ave(Mag), Ave, StdDev, Integral" <<endl;
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::writeProbesFileHeader
(
    writeFile& writer
) const
{
    writer.writeCommented(writer.file(), "Probes for Region type : ");
    writer.file() << regionTypeNames_[regionType_] << " "
                  << regionName_ << endl;

    writer.writeCommented(writer.file(), "Time");

    forAll(probeLocations_,probei)
    {
        writer.file()  << tab << probeLocations_[probei];
    }
    writer.file()<< endl;
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::calcGeomFilterFaces
(
    const vectorField& Cf,
    const PtrList<dictionary>& geomFilters
)
{
    boolList markedFaces(Cf.size(), true);
    forAll(geomFilters, filteri)
    {
        const dictionary& filterDict(geomFilters[filteri]);
        word type(filterDict.lookup("type"));
        bool inside(filterDict.lookupOrDefault<Switch>("inside",true));

        if (type == "sphere")
        {
            point sc(filterDict.lookup("centre"));
            scalar sr(readScalar(filterDict.lookup("radius")));

            forAll(Cf, facei)
            {
                const point centre = Cf[facei];
                scalar offset = mag(sc - centre);
                if (offset <= sr)
                {
                    markedFaces[facei] = inside ? true : false;
                }
                else
                {
                    markedFaces[facei] = inside ? false : true;
                }
            }
        }
        else if (type == "box")
        {
            treeBoundBox box(filterDict.lookup("boundBox"));
            forAll(Cf, facei)
            {
                const point centre = Cf[facei];
                if (box.contains(centre))
                {
                    markedFaces[facei] = inside ? true : false;
                }
                else
                {
                    markedFaces[facei] = inside ? false : true;
                }
            }
        }
    }
    filterFaces_.clearStorage();
    filterFaces_.reserve(Cf.size());

    forAll(markedFaces, facei)
    {
        if (markedFaces[facei])
        {
            filterFaces_.append(facei);
        }
    }

    filterFaces_.shrink();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCalcs::surfaceFieldCalc::surfaceFieldCalc
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldCalc(name, runTime, dict, typeName),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operations_(dict.lookup("operations")),
    probeLocations_(dict.lookupOrDefault("probeLocations",pointField(0))),
    probeIndex_(probeLocations_.size(), -1),
    nFaces_(0),
    filterFaces_(),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);

    forAll(operationsWriter_, opI)
    {
        writeFile& writer(operationsWriter_[opI]);
        writeFileHeader(writer);
    }

    forAll(probesWriter_, opI)
    {
        if (probesWriter_[opI].valid())
        {
            writeFile& writer(probesWriter_[opI]());
            writeProbesFileHeader(writer);
        }
    }
}


Foam::functionObjects::fieldCalcs::surfaceFieldCalc::surfaceFieldCalc
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldCalc(name, obr, dict, typeName),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operations_(dict.lookup("operations")),
    nFaces_(0),
    filterFaces_(),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);

    forAll(operationsWriter_, opI)
    {
        writeFile& writer(operationsWriter_[opI]);
        writeFileHeader(writer);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCalcs::surfaceFieldCalc::~surfaceFieldCalc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldCalcs::surfaceFieldCalc::read
(
    const dictionary& dict
)
{
    fieldCalc::read(dict);
    initialise(dict);

    return true;
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::calculate
(
    const meshedSurfRef& surfToWrite
)
{
    // Many operations use the Sf field
    vectorField Sf,Cf;
    if (needsGeom())
    {
        if (regionType_ == stSurface)
        {
            const surfMesh& s = dynamicCast<const surfMesh>(obr());
            Sf = s.Sf();
            Cf = s.Cf();
        }
        else if (surfacePtr_.valid())
        {
            Sf = surfacePtr_().Sf();
            Cf = surfacePtr_().Cf();
        }
        else
        {
            Sf = filterField(mesh_.Sf());
            Cf = filterField(mesh_.Cf());
        }
    }

    forAll(operations_, opI)
    {
        const dictionary& opDict(operations_[opI]);
        word operation(opDict.lookup("operation"));
        writeFile& opWriter(operationsWriter_[opI]);
        autoPtr<writeFile>& prWriter = probesWriter_[opI];

        const PtrList<dictionary> geomFilters
        (
            opDict.found("subsetDict") ?
            opDict.lookup("subsetDict") :  PtrList<dictionary>(0)
        );
        calcGeomFilterFaces(Cf,geomFilters);

        switch(calcTypeNames_.read(IStringStream(operation)()))
        {
            case ctField:
                field(opDict,Sf,surfToWrite,opWriter,prWriter);
            break;
            case ctSwirlAngle:
                swirlAngle(opDict,Cf,Sf,surfToWrite,opWriter,prWriter);
            break;
            case ctBulkSwirlAngle:
                bulkSwirlAngle(opDict,Cf,Sf,surfToWrite,opWriter,prWriter);
            break;
            case ctAxialDeviation:
                axialDeviation(opDict,Sf,surfToWrite,opWriter,prWriter);
            break;
            default:
                WarningInFunction
                    << operation << " is an invalid operation." << endl;
            break;
        }
    }
}

bool Foam::functionObjects::fieldCalcs::surfaceFieldCalc::write()
{
    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    fieldCalc::write();

    forAll(operationsWriter_, opI)
    {
        Ostream& os = operationsWriter_[opI].file();
        operationsWriter_[opI].writeTime(os);
    }

    Log << "    total area = " << totalArea_ << endl;

    // Faces and points for surface format (if specified)
    faceList faces;
    pointField points;

    if (surfaceWriterPtr_.valid())
    {
        if (regionType_ == stSurface || surfacePtr_.valid())
        {
            combineSurfaceGeometry(faces, points);
        }
        else
        {
            combineMeshGeometry(faces, points);
        }
    }

    meshedSurfRef surfToWrite(points, faces);

    // Perform field calculations
    calculate(surfToWrite);

    Log << endl;

    return true;
}


// ************************************************************************* //
