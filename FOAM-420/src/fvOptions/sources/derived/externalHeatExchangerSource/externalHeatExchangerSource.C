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
    (c) 2010-2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "externalHeatExchangerSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "regionSplit/regionSplit.H"
#include "primitives/strings/stringOps/stringOps.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"
#include "primitives/strings/string/string.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(externalHeatExchangerSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        externalHeatExchangerSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fv::externalHeatExchangerSource::zoneBoundaryFaces(const labelList& regionCells)
{

    DynamicList<label> zoneBoundFaces(0);

    boolList zoneCells(mesh_.nCells(), false);

    forAll(regionCells, i)
    {
        zoneCells[regionCells[i]] = true;
    }

    const labelList& nei = mesh_.neighbour();
    const labelList& own = mesh_.faceOwner();

    forAll(nei, i)
    {
        if (zoneCells[own[i]] && !zoneCells[nei[i]])
        {
            zoneBoundFaces.append(i);
        }
        else if (zoneCells[nei[i]] && !zoneCells[own[i]])
        {
            zoneBoundFaces.append(i);
        }
    }

    // check on processor patches
    boolList neiZoneID(mesh_.nFaces()-mesh_.nInternalFaces(),false);

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        if (zoneCells[own[faceI]])
        {
            neiZoneID[faceI-mesh_.nInternalFaces()] =  true;
        }
    }

    boolList swappedZoneID(neiZoneID);

    syncTools::swapBoundaryFaceList(mesh_, swappedZoneID);

    const polyPatchList& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
        const processorPolyPatch& procPatch
            = refCast<const processorPolyPatch>(patches[patchI]);

            forAll(procPatch, i)
            {
                label meshFaceI = procPatch.start() + i;
                label bFaceI =  meshFaceI - mesh_.nInternalFaces();

                if
                (
                    !swappedZoneID[bFaceI] && neiZoneID[bFaceI]
                )
                {
                    zoneBoundFaces.append(meshFaceI);
                }
            }
        }
        else
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            forAll(patch, i)
            {
                label meshFaceI = patch.start() + i;
                label bFaceI =  meshFaceI - mesh_.nInternalFaces();
                if (neiZoneID[bFaceI])
                {
                    zoneBoundFaces.append(meshFaceI);
                }
            }
        }
    }

    return labelList(zoneBoundFaces, true);
}

void Foam::fv::externalHeatExchangerSource::heConstruction()
{
    Info<<"External heat exchanger model autozoning:"<<endl;

    // 1) split separate regions OK!
    boolList boundaryFace(mesh_.faces().size(), false);

    // get faces bounding the cell zone
    const labelList& fList = zoneBoundaryFaces(cells_);

    forAll(fList, i)
    {
        const label& facei = fList[i];
        boundaryFace[facei] = true;
    }

    syncTools::syncFaceList(mesh_, boundaryFace, orEqOp<bool>());

    regionSplit cellMark(mesh_, boundaryFace);

    // 2) store tubes only in the list OK!

    scalar maxID = max(cellMark);
    reduce(maxID, maxOp<scalar>());

    scalar nRegions = maxID+1;

    List<labelList> allRegions_(nRegions);

    forAll(cells_, i)
    {
        label regionID = cellMark[cells_[i]];
        allRegions_[regionID].append(cells_[i]); // store cell index
    }

    // clean list from regions not belonging to tubes
    List<labelList> unsortedTubes_(nRows_*nCol_);

    label index = 0;
    forAll(allRegions_, i)
    {
        // parallel handling
        scalar gRegionSize = allRegions_[i].size();
        reduce(gRegionSize, sumOp<scalar>());

        // make sure it is a tube
        if (gRegionSize!=0)
        {
            forAll(allRegions_[i], j)
            {
                unsortedTubes_[index].append(allRegions_[i][j]);
            }
            index ++;
        }
    }

    // 3) calculate origin for each tube and store in a list

    originTubes_.setSize(nRows_*nCol_);

    forAll(unsortedTubes_, j)
    {
        labelList tubej = unsortedTubes_[j];
        originTubes_[j] = tubeOPoint(tubej);
    }

    // 4) sort by row based on local origin OK

    List<point> originTubesLocalCS_(nRows_*nCol_);

    //transform points to the HE global coordinate system to perform sorting
    forAll(originTubes_, pointI)
    {
        // point localPoint = coordSysPtr_->R().invTransform(originTubes_[pointI]);
        point localPoint = coordSysPtr_->localPosition(originTubes_[pointI]);
        originTubesLocalCS_[pointI] = localPoint;
    }

    // sort rows - e3 dir
    label lastIndex = originTubesLocalCS_.size()-1;
    bSortTubes(unsortedTubes_, originTubesLocalCS_, 2, 0, lastIndex);

    // sort columns - e1 dir
    for (label i = 0; i < originTubesLocalCS_.size(); i+=nCol_)
    {
        label iStart = i;
        label iEnd = i+3;
        bSortTubes(unsortedTubes_, originTubesLocalCS_, 0, iStart, iEnd);
    }

    sortedTubes_ = unsortedTubes_;

    Info<<"Sorted tubes in global coordinate system x y z: "<<endl;
    forAll(originTubes_, i)
    {
        Info<<originTubes_[i]<<endl;
    }

    // 5) store monitoring faces OK

    monitorTubes_.setSize(nRows_*nCol_);
    const faceZone& fZoneInlet = mesh_.faceZones()[inletID_];
    const faceZone& fZoneOutlet = mesh_.faceZones()[outletID_];

    forAll(sortedTubes_, j)
    {
        // get faces bounding the tube j
        labelList tubej = sortedTubes_[j];
        const labelList& fListTube = zoneBoundaryFaces(tubej);

        // store the face if it belongs to the inlet or outlet faceSet
        forAll(fListTube, i)
        {
            // get face
            const label& faceTubei = fListTube[i];

            // check if it belongs to inlet faceSet
            forAll(fZoneInlet, i)
            {
                label fZoneInleti = fZoneInlet[i];

                if (mesh_.isInternalFace(fZoneInleti))
                {
                    if (faceTubei == fZoneInleti)
                    {
                        monitorTubes_[j].append(faceTubei);
                    }
                }
                else
                {
                    label facePatchId = mesh_.boundaryMesh().whichPatch(fZoneInleti);
                    const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
                    if (isA<coupledPolyPatch>(pp))
                    {
                        if (refCast<const coupledPolyPatch>(pp).owner())
                        {
                            if (faceTubei == fZoneInleti)
                            {
                                monitorTubes_[j].append(faceTubei);
                            }
                        }
                    }
                }
            }

            // check if it belongs to outlet faceSet
            forAll(fZoneOutlet, i)
            {
                label fZoneOutleti = fZoneOutlet[i];

                if (mesh_.isInternalFace(fZoneOutleti))
                {
                    if (faceTubei == fZoneOutleti)
                    {
                        monitorTubes_[j].append(faceTubei);
                    }
                }
                else
                {
                    label facePatchId = mesh_.boundaryMesh().whichPatch(fZoneOutleti);
                    const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
                    if (isA<coupledPolyPatch>(pp))
                    {
                        if (refCast<const coupledPolyPatch>(pp).owner())
                        {
                            if (faceTubei == fZoneOutleti)
                            {
                                monitorTubes_[j].append(faceTubei);
                            }
                        }
                    }
                }
            }
        }
    }

    //6) store tube length and width OK

    lengthTubes_.setSize(nRows_*nCol_);
    widthTubes_.setSize(nRows_*nCol_);

    forAll(sortedTubes_, j)
    {
        labelList tubej = sortedTubes_[j];
        scalarList dimensions = calcTubeDim(tubej);

        lengthTubes_[j] = dimensions[0];
        widthTubes_[j] = dimensions[1];
    }

    // calculate row width pCB if not provided in dict
    if (pCB_ < 0)
    {
        scalarList zoneDim = calcTubeDim(cells_);
        pCB_ = zoneDim[1];

    }

}


Foam::point Foam::fv::externalHeatExchangerSource::tubeOPoint(const labelList& tubeCells)
{
    vector tubeLocalCS(3);
    tubeLocalCS[0] = GREAT;   //min
    tubeLocalCS[1] = GREAT;   //min
    tubeLocalCS[2] =-GREAT;   //max

    // loop over the cells of a single tube
    forAll(tubeCells, i)
    {
        // get label of each cell point
        const labelList& cellPoints = mesh_.cellPoints()[tubeCells[i]];

        // loop over the cell points
        forAll(cellPoints, cPtI)
        {
            const label pointI = cellPoints[cPtI];
            const point& cc = mesh_.points()[pointI];

            //transform point in HE coord system e1 e2 and e3
            // const point ccHe = coordSysPtr_->R().invTransform(cc);
            const point ccHe = coordSysPtr_->localPosition(cc);

            //min e1
            if (ccHe[0]<tubeLocalCS[0])
            {
                tubeLocalCS[0] = ccHe[0];
            }
            //min e2
            if (ccHe[1]<tubeLocalCS[1])
            {
                tubeLocalCS[1] = ccHe[1];
            }
            //max e3
            if (ccHe[2]>tubeLocalCS[2])
            {
                tubeLocalCS[2] = ccHe[2];
            }
        }
    }

    reduce( tubeLocalCS[0], minOp<scalar>());
    reduce( tubeLocalCS[1], minOp<scalar>());
    reduce( tubeLocalCS[2], maxOp<scalar>());

    // convert to global coordinate system x y z
    //point tubeCS = coordSysPtr_->R().transform(tubeLocalCS);
    point tubeCS = coordSysPtr_->globalPosition(tubeLocalCS);

    return tubeCS;
}

void Foam::fv::externalHeatExchangerSource::bSortTubes
(
List<labelList>& unsortedTubes_,
List<point>& originPoints,
label dir,
label iStart,
label iEnd
)
{
    for (label i = iStart; i < iEnd ; i++)
    {
        for (label j = iStart; j < iEnd; j++)
        {
            if (originPoints[j][dir] > originPoints[j+1][dir])
            {
                // sort points in HE coord system e1 e2 e3
                point tempLocalPoint = originPoints[j];
                originPoints[j] = originPoints[j+1];
                originPoints[j+1] = tempLocalPoint;

                // sort points in global coord system x y z
                point tempGlobalPoint = originTubes_[j];
                originTubes_[j] = originTubes_[j+1];
                originTubes_[j+1] = tempGlobalPoint;

                // sort cells
                labelList tempRow = unsortedTubes_[j];
                unsortedTubes_[j] = unsortedTubes_[j+1];
                unsortedTubes_[j+1] = tempRow;
            }
        }
    }
}

Foam::scalarList Foam::fv::externalHeatExchangerSource::calcTubeDim(const labelList& tubeCells)
{

    scalarList tubeDim(2);

    PackedBoolList sourcePts(mesh_.nPoints(), 0);

    forAll(tubeCells, i)
    {
        label celli = tubeCells[i];
        const labelList& cellPoints = mesh_.cellPoints()[celli];

        forAll(cellPoints, cPtI)
        {
            sourcePts.set(cellPoints[cPtI], 1);
        }
    }

    scalar minLength = GREAT;
    scalar maxLength = -GREAT;
    scalar minWidth = GREAT;
    scalar maxWidth = -GREAT;

    forAll(mesh_.points(), pointi)
    {
        if (sourcePts.get(pointi) == 1)
        {
            vector v(mesh_.points()[pointi]);

            //length
            scalar lparallel = (v & coordSysPtr_->e2());
            minLength = min(minLength, lparallel);
            maxLength = max(maxLength, lparallel);

            //width
            scalar wparallel = (v & coordSysPtr_->e1());
            minWidth = min(minWidth, wparallel);
            maxWidth = max(maxWidth, wparallel);

        }
    }

    reduce(minLength, minOp<scalar>());
    reduce(maxLength, maxOp<scalar>());
    reduce(minWidth, minOp<scalar>());
    reduce(maxWidth, maxOp<scalar>());

    tubeDim[0] = mag(maxLength - minLength);
    tubeDim[1] = mag(maxWidth - minWidth);

    return tubeDim;

}

Foam::tmp<Foam::volScalarField>
Foam::fv::externalHeatExchangerSource::getRho() const
{

    if (obr_.found("rho"))
    {
        return(obr_.lookupObject<volScalarField>("rho"));
    }
    else
    {
        dimensionedScalar rho
        (
            "rho",
            obr_.lookupObject<dictionary>("transportProperties")
                .lookup("rho")
        );

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    obr_
                ),
                mesh_,
                rho
            )
        );
    }
}

Foam::scalarList Foam::fv::externalHeatExchangerSource::tubeMFluxTemp
(
    const fvMesh& mesh,
    const labelList& faces
)
{
    scalar fsum = 0;
    const surfaceScalarField& phi
        = mesh.lookupObject<surfaceScalarField>(phiName_);

    surfaceScalarField rhof(fvc::interpolate(getRho()));

    scalar Tsum = 0;
    const volScalarField& T
        = mesh.lookupObject<volScalarField>(TName_);

    surfaceScalarField Tf(fvc::interpolate(T));

    bool incomp = false;
    if (phi.dimensions() == dimVelocity * dimArea)
    {
        incomp = true;
    }

    forAll(faces, i)
    {
        label meshFaceI = faces[i];
        if (meshFaceI < mesh.nInternalFaces())
        {
            if (incomp)
            {
                fsum += rhof[meshFaceI]*mag(phi[meshFaceI]);

                Tsum += rhof[meshFaceI]*mag(phi[meshFaceI])*Tf[meshFaceI];
            }
            else
            {
                fsum += mag(phi[meshFaceI]);

                Tsum += mag(phi[meshFaceI])*Tf[meshFaceI];
            }
        }
        else
        {
            label patchI = mesh.boundaryMesh().whichPatch(meshFaceI);
            const polyPatch& patch = mesh.boundaryMesh()[patchI];
            if (!isA<emptyPolyPatch>(patch))
            {
                label bfI = meshFaceI - patch.start();

                if (incomp)
                {
                    fsum += rhof.boundaryField()[patchI][bfI]
                        *mag(phi.boundaryField()[patchI][bfI]);

                    Tsum += rhof.boundaryField()[patchI][bfI]
                        *mag(phi.boundaryField()[patchI][bfI])
                        *Tf.boundaryField()[patchI][bfI];
                }
                else
                {
                    fsum += mag(phi.boundaryField()[patchI][bfI]);

                    Tsum += mag(phi.boundaryField()[patchI][bfI])
                        *Tf.boundaryField()[patchI][bfI];
                }
            }
        }
    }
    reduce(fsum, sumOp<scalar>());
    reduce(Tsum, sumOp<scalar>());

    scalarList monitoredQuantities(2);

    monitoredQuantities[0] = fsum; // mass flow rate
    monitoredQuantities[1] = Tsum/fsum-273.15; // flux-ave temperature in deg C

    return monitoredQuantities;
}


Foam::scalar Foam::fv::externalHeatExchangerSource::finDensity(label rowIndex)
{

    scalar pFD;

    if ((rowIndex - 1) <= 4)
    {
        pFD = 38.0;
    }
    else if ((rowIndex - 1) <= 7)
    {
        pFD = 36.0;
    }
    else if ((rowIndex - 1) <= 10)
    {
        pFD = 34.0;
    }
    else if ((rowIndex - 1) <= 20)
    {
        pFD = 30.0;
    }
    else if ((rowIndex - 1) <= 25)
    {
        pFD = 26.0;
    }
    else if ((rowIndex - 1) <= 30)
    {
        pFD = 24.0;
    }
    else if ((rowIndex - 1) <= 35)
    {
        pFD = 22.0;
    }
    else if ((rowIndex - 1) <= 40)
    {
        pFD = 20.0;
    }
    else if ((rowIndex - 1) <= 50)
    {
        pFD = 18.0;
    }
    else if ((rowIndex - 1) <= 60)
    {
        pFD = 16.0;
    }
    else if ((rowIndex - 1) <= 70)
    {
        pFD = 14.0;
    }
    else
    {
        pFD = 12.0;
    }

    return pFD;
}



void Foam::fv::externalHeatExchangerSource::readHRfile(pointField& rowPoints, scalarField& rowHR)
{
    DynamicList<point> dPoints;
    DynamicList<scalar> dHR;

    // read done by each processor to handle parallel
    // this can be also done from master only but a split of information
    // to different processors is required after reading

    fileName expandedFile(fName_);

    // wait till the file is created
    fileState
    (
        fileName(fName_),
        true
    );

    IFstream is(expandedFile.expand());

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot read file."
            << exit(FatalIOError);
    }

    Info<<"Reading .dat file ..."<<endl;

    while (is.good())
    {
        string line;
        is.getLine(line);

            // skip first line
        if (!stringOps::upper(line).find('X'))
        {
            continue;
        }

        label n = 0;
        std::size_t pos = 0;
        DynamicList<string> splitted;

        label nEntries = 3;

        while ((pos != std::string::npos) && (n <= nEntries))
        {
            std::size_t nPos = line.find(',', pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos = nPos;
                n++;
            }
            else
            {
                splitted.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
                n++;
            }
        }

        if (splitted.size() <= 1)
        {
            break;
        }

        scalar x = readScalar(IStringStream(splitted[0])());
        scalar y = readScalar(IStringStream(splitted[1])());
        scalar z = readScalar(IStringStream(splitted[2])());
        scalar Q = readScalar(IStringStream(splitted[3])());

        dPoints.append(point(x, y, z));
        dHR.append(Q);
    }

    dPoints.shrink();
    dHR.shrink();

    rowPoints.setSize(dPoints.size());
    rowPoints.transfer(dPoints);
    dPoints.clear();

    rowHR.setSize(dHR.size());
    rowHR.transfer(dHR);
    dHR.clear();

    Info<<"Reading completed."<<endl;

}

void Foam::fv::externalHeatExchangerSource::fileState(const fileName& f, bool shouldExist) const
{
    label tooLong = 1;
    label iters = 0;

    if (shouldExist)
    {
        while (!isFile(f))
        {
            sleep(1);
            iters++;
            if (iters > tooLong)
            {
                iters = 0;
                Pout<< "Waiting for file: "
                     << f << " to be created." << endl;
            }
        }
    }
    else
    {
        while (isFile(f))
        {
            sleep(1);
            iters++;
            if (iters > tooLong)
            {
                iters = 0;
                Pout<< "Waiting for file: "
                     << f << " to be deleted." << endl;
            }
        }

    }
}

void Foam::fv::externalHeatExchangerSource::mapRowHR
(
label firstTube,
label lastTube,
point csPoint,
pointField& rowPoints,
scalarField& rowHR
)
{
    Info<<"Row origin: "<<csPoint<<". Tubes: "<<firstTube<<" to "<<lastTube<<"."<<endl;

    // convert point cloud form local row coord system to x y x
    pointField rowPointsGCS(rowPoints.size());
    forAll(rowPoints, i)
    {
        //generic point + local tube origin (this is already in GCS)
        //rowPointsGCS[i] = coordSysPtr_->R().transform(rowPoints[i]) + csPoint;
        rowPointsGCS[i] = coordSysPtr_->globalPosition(rowPoints[i]) + csPoint;

        // Info<<rowPointsGCS[i][0]<<" , "<<rowPointsGCS[i][1]<<" , "<<rowPointsGCS[i][2]<<endl;
    }

    treeBoundBox bb(rowPointsGCS);
    Random rndGen(123456);
    bb = bb.extend(rndGen, 1e-4);

    indexedOctree<treeDataPoint> tree
    (
        treeDataPoint(rowPointsGCS),
        bb,         // overall search domain
        20,         // max levels
        100,        // maximum ratio of cubes v.s. cells
        100.0       // max
    );

    scalarField& Qhr = Qhr_->primitiveFieldRef();

    // map over 4 tubes
    for (label i = firstTube; i <= lastTube; i++)
    {
        labelList tubei = sortedTubes_[i];
        forAll(tubei, j)
        {
            label celli = tubei[j];
            const point& cc = mesh_.cellCentres()[celli];

            pointIndexHit info = tree.findNearest
            (
                cc,
                GREAT
            );

            if (info.hit())
            {
                label index = info.index();
                Qhr[celli] = rowHR[index];
            }
            else
            {
                Pout<<"Cannot find point: "<<cc<<endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::externalHeatExchangerSource::externalHeatExchangerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    TName_(coeffs_.lookupOrDefault<word>("TName", "T")),
    phiName_(coeffs_.lookupOrDefault<word>("phiName", "phi")),
    nRows_(coeffs_.lookupOrDefault<label>("rowsNumber", 80)),
    nCol_(4),
    inletZoneName_(coeffs_.lookup("tubeInletFaceZone")),
    outletZoneName_(coeffs_.lookup("tubeOutletFaceZone")),
    inletID_(mesh_.faceZones().findZoneID(inletZoneName_)),
    outletID_(mesh_.faceZones().findZoneID(outletZoneName_)),
    solutionFrequency_(coeffs_.lookupOrDefault<label>("solutionFrequency", 25)),
    firstIter_(true),
    linearised_(coeffs_.lookupOrDefault<Switch>("linearisedSource", true)),
    fName_(coeffs_.lookupOrDefault<fileName>("file", "csvHR.dat")),
    appPath_(coeffs_.lookupOrDefault<fileName>("appPath", "/home/user/app/HeatExchanger.sh")),
    mathPath_(coeffs_.lookupOrDefault<fileName>("mathPath", "/home/user/Matlab/v92")),
    pCB_(coeffs_.lookupOrDefault<scalar>("pCB", -1)),
    passes_(coeffs_.lookupOrDefault<label>("passes", 2)),
    pCC_(coeffs_.lookupOrDefault<scalar>("pCC", 5.048)),
    chargePressure_(coeffs_.lookupOrDefault<scalar>("chargePressure", 450000)),
    pFA_(coeffs_.lookupOrDefault<scalar>("pFA", 0.00216)),
    pFB_(coeffs_.lookupOrDefault<scalar>("pFB", 25)),
    pFC_(coeffs_.lookupOrDefault<scalar>("pFC", 0.00485)),
    ambientTemperature_(coeffs_.lookupOrDefault<scalar>("ambientTemperature", 30)),
    pCA_(coeffs_.lookupOrDefault<word>("pCA", "StandardReturn")),
    returnWidth_(coeffs_.lookupOrDefault<scalar>("returnWidth", 0.0)),
    pAA_(coeffs_.lookupOrDefault<scalar>("pAA", 0.0107)),
    pAB_(coeffs_.lookupOrDefault<scalar>("pAB", 0.0057)),
    pAC_(coeffs_.lookupOrDefault<scalar>("pAC", 0.002)),
    pAD_(coeffs_.lookupOrDefault<scalar>("pAD", 0.00135)),
    coordSysPtr_(),
    originTubes_(),
    sortedTubes_(),
    monitorTubes_(),
    lengthTubes_(),
    widthTubes_(),
    Qhr_()
{

    if (inletID_ < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << inletZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }
    if (outletID_ < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << outletZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::externalHeatExchangerSource::initialise()
{
    // set field name
    const surfaceScalarField& phi
        = obr_.lookupObject<surfaceScalarField>("phi");

    bool incomp = false;
    if (phi.dimensions() == dimVelocity * dimArea)
    {
        incomp = true;
    }

    if (!incomp)
    {
        const basicThermo& thermo =
            obr_.lookupObject<basicThermo>(basicThermo::dictName);
        fieldNames_.setSize(1, thermo.he().name());
    }
    else
    {
        fieldNames_.setSize(1,TName_);
    }

    applied_.setSize(1, false);

    read(this->coeffs());
    heConstruction();

    // init Qhr to 0
    Qhr_.reset
    (
        new volScalarField(
        IOobject
         (
            "Qhr",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
         mesh_,
        dimensionedScalar("Qhr",dimensionSet(1,-1,-3,0,0,0,0),0) // [W/m3]
     )
   );

    if (passes_ != 2)
    {
        passes_ = 2;
        WarningInFunction << " The model is designed to work with "
                          << " 2 cold tubes and 2 hot tubes per row."
                          << " If you require any extension,"
                          << " please contact ESI."<< endl;
    }

    return true;
}


void Foam::fv::externalHeatExchangerSource::correct()
{

    if (firstIter_ || (mesh_.time().timeIndex() % solutionFrequency_ == 0))
    {
        Info<<"***External Heat Exchanger Solver***"<<endl;

        // solve per row starting from last row
        label rowIndex = nRows_;
        for (label i = 0; i < originTubes_.size(); i+=nCol_)
        {
                // first and last tubes of each row
            label iStart = i;
            label iEnd = i+3; //nCol_-1;

                // local coord system per row i
            point csRowi = originTubes_[i];

            scalar pFD = finDensity(rowIndex);

            scalarList tubeCold1(2);
            scalarList tubeCold2(2);
            scalarList tubeHot1(2);
            scalarList tubeHot2(2);

            tubeCold1 = tubeMFluxTemp(mesh_, monitorTubes_[i+0]);
            tubeCold2 = tubeMFluxTemp(mesh_, monitorTubes_[i+1]);
            tubeHot1 = tubeMFluxTemp(mesh_, monitorTubes_[i+2]);
            tubeHot2 = tubeMFluxTemp(mesh_, monitorTubes_[i+3]);

            string length_s = Foam::name(lengthTubes_[i]);
            string pCB_s = Foam::name(pCB_);
            string passes_s = Foam::name(passes_);
            string pCC_s = Foam::name(pCC_);
            string chargePressure_s = Foam::name(chargePressure_);
            string chargeTemperature_s =
                "'[" + string(Foam::name(tubeCold1[1])) + ", "+string(Foam::name(tubeCold2[1]))
                +"; "+string(Foam::name(tubeHot1[1]))+", "+string(Foam::name(tubeHot2[1])) +"]'";
            string chargeFlow_s =
                "'[" + string(Foam::name(tubeCold1[0])) + ", "+string(Foam::name(tubeCold2[0]))
                +"; "+string(Foam::name(tubeHot1[0]))+", "+string(Foam::name(tubeHot2[0])) +"]'";
            string pFD_s = Foam::name(pFD);
            string pFA_s = Foam::name(pFA_);
            string pFB_s = Foam::name(pFB_);
            string pFC_s = Foam::name(pFC_);
            string pTA_s = Foam::name(ambientTemperature_);
            string pCA_s = pCA_;
            string returnWidth_s = Foam::name(returnWidth_);
            string pAA_s = Foam::name(pAA_);
            string pAB_s = Foam::name(pAB_);
            string pAC_s = Foam::name(pAC_);
            string pAD_s = Foam::name(pAD_);

            // systemCall
            if (Pstream::master())
            {
                Info<<"System call - row "<<rowIndex<<"."<<endl;

                string inPar = length_s + ", " + pCB_s + ", " + passes_s + ", " + pCC_s + ", "
                    + chargePressure_s + ", " + chargeTemperature_s + ", " + chargeFlow_s + ", "
                    + pFD_s + ", " + pFA_s + ", " + pFB_s + ", " + pFC_s + ", " + pTA_s + ", "
                    + pCA_ + ", " + returnWidth_s + ", " + pAA_s + ", " + pAB_s + ", "
                    + pAC_s + ", " + pAD_s;

                string command = appPath_ + " " + mathPath_ + " " + inPar;

                Info<<"Call arguments: "<<inPar<<endl;
                Foam::system(command);

            }
            // .dat file is ready and can be read

            pointField rowPoints;
            scalarField rowHR;

            label readSync = 0;
            reduce(readSync, sumOp<label>());

            readHRfile(rowPoints,rowHR);

            label readSync2 = 0;
            reduce(readSync2, sumOp<label>());

            fileName targetFile(("HR_Tube_" + Foam::name(rowIndex) +".cvs"));

            if (Pstream::master())
            {
                mv(fName_,targetFile);
            }

            Info<<"Mapping row "<<rowIndex<<endl;
            mapRowHR(iStart,iEnd,csRowi,rowPoints,rowHR);
            Info<<"Mapping row "<<rowIndex<<" completed."<<endl;

            // update row index for fin density calculation
            rowIndex--;
        }

        firstIter_ = false;
    }

}


void Foam::fv::externalHeatExchangerSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label
)
{

    if (linearised_)
    {
        const volScalarField rho( getRho() );

        const incompressible::turbulenceModel & turbModel =
            obr_.lookupObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        volScalarField Cp( turbModel.transport().Cp() );
        const volScalarField& Qhr = Qhr_;
        const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

        dimensionedScalar Tamb
        (
            "Tamb",
            dimTemperature,
            ambientTemperature_+273.15
        );

        volScalarField HTC
        (
            IOobject
            (
                "HTC",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(1,-1,-3,-1,0,0,0), 0),
            "calculated"
        );

        forAll(cells_, i)
        {
            label celli = cells_[i];

            scalar dTmax = -0.1;
            scalar dT = Tamb.value() - T[celli];

            if (dT >= dTmax) // always negative
            {
                dT = dTmax;
            }

            HTC[celli] = max(Qhr[celli]/dT, 0.0);
        }

        eqn +=HTC/(rho*Cp)*Tamb;

        eqn +=fvm::Sp(-HTC/(rho*Cp),eqn.psi());

    }
    else
    {
        const volScalarField rho( getRho() );

        const incompressible::turbulenceModel & turbModel =
            obr_.lookupObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        volScalarField Cp( turbModel.transport().Cp() );

        const volScalarField& Qhr = Qhr_;

        eqn +=Qhr/(rho*Cp);
    }

}

void Foam::fv::externalHeatExchangerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label
)
{
    if (linearised_)
    {
        const volScalarField& Qhr = Qhr_;
        const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

        const basicThermo& thermo =
            obr_.lookupObject<basicThermo>(basicThermo::dictName);

        volScalarField Cp( thermo.Cp() );

        dimensionedScalar Tamb
        (
            "Tamb",
            dimTemperature,
            ambientTemperature_+273.15
        );

        volScalarField HTC
        (
            IOobject
            (
                "HTC",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(1,-1,-3,-1,0,0,0), 0),
            "calculated"
        );

        forAll(cells_, i)
        {
            label celli = cells_[i];

            scalar dTmax = -0.1;
            scalar dT = Tamb.value() - T[celli];

            if (dT >= dTmax) // always negative
            {
                dT = dTmax;
            }

            HTC[celli] = max(Qhr[celli]/dT, 0.0);
        }

        eqn +=HTC*Tamb;

        eqn +=fvm::Sp(-HTC/Cp,eqn.psi());
    }
    else
    {
        const volScalarField& Qhr = Qhr_;

        eqn +=Qhr;
    }

}

bool Foam::fv::externalHeatExchangerSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coordSysPtr_ = coordinateSystem::New(mesh_, coeffs_);

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
