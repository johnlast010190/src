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
    (c) 2014-1014 Esi Ltd.
\*---------------------------------------------------------------------------*/

#include "hexReport/hexReport.H"
#include "meshes/meshShapes/cellMatcher/hexMatcher.H"
#include "meshes/meshShapes/cellMatcher/wedgeMatcher.H"
#include "meshes/meshShapes/cellMatcher/prismMatcher.H"
#include "meshes/meshShapes/cellMatcher/pyrMatcher.H"
#include "meshes/meshShapes/cellMatcher/tetWedgeMatcher.H"
#include "meshes/meshShapes/cellMatcher/tetMatcher.H"
#include "regionSplit/regionSplit.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(hexReport, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hexReport::write()
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<<nl;
    Info<<"Final Mesh Summary Report"<<nl;
    Info<<"-------------------------"<<nl;
    Info<<endl;

    Info<< "Mesh stats" << nl
        << "----------" << nl
        << nl
        << "    points:              "
        << returnReduce(mesh.points().size(), sumOp<label>()) << nl;

    Info<< "    edges:               "
        << returnReduce(mesh.nEdges(), sumOp<label>()) << nl;

    label nFaces = returnReduce(mesh.faces().size(), sumOp<label>());
    label nIntFaces = returnReduce(mesh.faceNeighbour().size(), sumOp<label>());
    label nCells = returnReduce(mesh.cells().size(), sumOp<label>());

    label nPatches = 0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        if (!isA<processorPolyPatch>(pp))
        {
            nPatches++;
        }
    }
    regionSplit rs(mesh);

    scalar totalVolume = 0;
    forAll(mesh.cells(), cellI)
    {
        totalVolume += mesh.cellVolumes()[cellI];
    }
    totalVolume = returnReduce(totalVolume, sumOp<scalar>());

    //Calculate mesh bounding box
    const boundBox& globalBb = mesh.bounds();

    Info<< "    faces:               " << nFaces << nl
        << "    internal faces:      " << nIntFaces << nl
        << "    cells:               " << nCells << nl
        << "    layers cells:        " << numLayerCells_ << nl
        << "    layer coverage:      " << layerCoverage_ <<" %"<< nl
        << "    boundary patches:    " << nPatches << nl
        << "    point zones:         " << mesh.pointZones().size() << nl
        << "    face zones:          " << mesh.faceZones().size() << nl
        << "    cell zones:          " << mesh.cellZones().size() << nl
        << "    mesh regions:        " << rs.nRegions() << nl
        << "    domain volume:       " << setprecision(5) << totalVolume << nl
        << "    domain bounding box: "
        << globalBb.min() << " " << globalBb.max() << nl
        << endl;

    //Calculate zone statistics
    mesh.checkZones();

    // Construct shape recognizers
    hexMatcher hex;
    prismMatcher prism;
    wedgeMatcher wedge;
    pyrMatcher pyr;
    tetWedgeMatcher tetWedge;
    tetMatcher tet;

    // Counters for different cell types
    label nHex = 0;
    label nWedge = 0;
    label nPrism = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTetWedge = 0;
    label nUnknown = 0;

    Map<label> polyhedralFaces;

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (hex.isA(mesh, cellI))
        {
            nHex++;
        }
        else if (tet.isA(mesh, cellI))
        {
            nTet++;
        }
        else if (pyr.isA(mesh, cellI))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, cellI))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, cellI))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, cellI))
        {
            nTetWedge++;
        }
        else
        {
            nUnknown++;
            polyhedralFaces(mesh.cells()[cellI].size())++;
        }
    }

    reduce(
        std::tie(
            nHex,
            nPrism,
            nWedge,
            nPyr,
            nTetWedge,
            nTet,
            nUnknown
        ),
        UniformParallelOp<sumOp<label>, 7>{}
    );

    Info<< "Overall number of cells of each type" << nl
        << "------------------------------------"<< nl
        << endl;
    Info<< "    hexahedra:  " << nHex
        << " (" <<100*scalar(nHex)/nCells << setprecision(3)<< " %)"<<endl;
    Info<< "    prisms:     " << nPrism
        << " (" <<100*scalar(nPrism)/nCells << setprecision(3)<< " %)"<<endl;
    Info<< "    wedges:     " << nWedge
        << " (" <<100*scalar(nWedge)/nCells << setprecision(3)<< " %)"<<endl;
    Info<< "    pyramids:   " << nPyr
        << " (" <<100*scalar(nPyr)/nCells << setprecision(3)<< " %)"<<endl;
    Info<< "    tet wedges: " << nTetWedge
        << " (" <<100*scalar(nTetWedge)/nCells << setprecision(3)<< " %)"<<endl;
    Info<< "    tetrahedra: " << nTet
        << " (" <<100*scalar(nTet)/nCells << setprecision(3)<< " %)"<<endl;
    Info<< "    polyhedra:  " << nUnknown
        << " (" <<100*scalar(nUnknown)/nCells << setprecision(3)<< " %)"<<endl;

    Info<< endl;

    //cell level information
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    label nLevels = gMax(cellLevel);

    if (nLevels > -1)
    {
        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
        labelList nCells(nLevels+1, 0);

        forAll(cellLevel, cellI)
        {
            nCells[cellLevel[cellI]]++;
        }

        Pstream::listCombineReduce(nCells, plusOp<label>());

        Info<< "Cells per refinement level" << endl;
        Info<< "--------------------------" << endl;
        Info<<endl;

        forAll(nCells, levelI)
        {
            scalar edgeLen = 1000*edge0Len / pow(2., levelI);
            Info<< "    " << levelI << '\t' << nCells[levelI]
                <<" (spacing: "<<edgeLen<<" mm)"<< nl;
        }
        Info<<endl;
    }


    if (patchLayerInfo_.size() > 0)
    {
        label maxPatchNameLen = 0;
        forAll(patchLayerInfo_, patchI)
        {
            word patchName = patchLayerInfo_[patchI].first();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }

        Info<< "Patch layer info" << endl;
        Info<< "----------------" << endl;
        Info<< nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "Key: " <<nl
            << setf(ios_base::left) << setw(maxPatchNameLen)
            << "First cell height (fch)" <<nl
            << setf(ios_base::left) << setw(maxPatchNameLen)
            << "Final layer height (flh)" <<nl
            << setf(ios_base::left) << setw(maxPatchNameLen)
            << "Total layer height (tlh)" <<nl
            << setf(ios_base::left) << setw(maxPatchNameLen)
            << "Expansion ratio (exp)" <<nl;

        Info<< nl
             << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
             << setw(0) << " faces    layers fch (m)  flh (m) "
             << " tlh (m)  exp      coverage (%)" << nl
             << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
             << setw(0) << " -----    ------ -------  -------  -------  "
             << "---      ------------" << endl;

        forAll(patchLayerInfo_, patchI)
        {
            Info<<endl;
            const scalarField& pLInfo = patchLayerInfo_[patchI].second();

            Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                << patchLayerInfo_[patchI].first() << setprecision(3)
                << " " << setw(8)
                << label(pLInfo[0])
                << " " << setw(6) << pLInfo[5]
                << " " << setw(8) << pLInfo[1]
                << " " << setw(8) << pLInfo[2]
                << " " << setw(8) << pLInfo[3]
                << " " << setw(8) << pLInfo[4]
                << "  " << setw(8) <<  pLInfo[6]
                << endl;
        }
        Info<<endl;
    }

    Info<< "Mesh quality" << endl;
    Info<< "------------" << endl;
    Info<<endl;

    bool report = true;
    label noFailedChecks = 0;
    if (mesh.checkPoints(report)) noFailedChecks++;
    if (mesh.checkUpperTriangular(report)) noFailedChecks++;
    if (mesh.checkCellsZipUp(report)) noFailedChecks++;
    if (mesh.checkFaceVertices(report)) noFailedChecks++;

    Info<<endl;

    if (mesh.checkClosedBoundary(report)) noFailedChecks++;
    if (mesh.checkClosedCells(report)) noFailedChecks++;
    if (mesh.checkFaceAreas(report)) noFailedChecks++;
    if (mesh.checkCellVolumes(report)) noFailedChecks++;
    if (mesh.checkFaceOrthogonality(report)) noFailedChecks++;
    if (mesh.checkFacePyramids(report)) noFailedChecks++;

    Info<<endl;

    if (noFailedChecks == 0)
    {
        Info<< "Mesh OK." << endl;
    }
    else
    {
        Info<< "Warning in " << noFailedChecks
            << " mesh checks." << endl;
    }
    Info<<endl;

    Info<< "Meshing time" << endl;
    Info<< "------------" << endl;
    Info<<endl;

    scalar totalTime = initialisationTime_ + refinementTime_ + dualTime_
        + snappingTime_ + layerTime_;

    Info<<"   Overall Meshing Time "<< setprecision(6)
        << totalTime << setprecision(3)
        <<" sec ("<<totalTime/3600<<" hours)."<<nl;
    Info<<endl;

    if (totalTime < SMALL)
    {
        return;
    }

    Info<<"   Performed mesh initialisation in "<< setprecision(6)
        << initialisationTime_ << setprecision(3)
        <<" sec (" << 100.0*initialisationTime_/totalTime << " %)."<<nl;

    Info<<"   Performed mesh refinement in " << setprecision(6)
        << refinementTime_ << setprecision(3)
        <<" sec (" << 100.0*refinementTime_/totalTime << " %)."<<nl;

    if (dualTime_ > SMALL)
    {
        if (meshRefiner_.controller().algorithm() == meshControl::DUAL)
        {
            Info<<"   Performed dual setup in " << setprecision(6)
                << dualTime_ << setprecision(3)
                <<" sec (" << 100.0*dualTime_/totalTime << " %)."<<nl;
        }
        else
        {
            Info<<"   Performed extrude setup in " << setprecision(6)
                << dualTime_ << setprecision(3)
                <<" sec (" << 100.0*dualTime_/totalTime << " %)."<<nl;
        }
    }

    Info<<"   Performed mesh snapping in " << setprecision(6)
        << snappingTime_ << setprecision(3)
        <<" sec (" << 100.0*snappingTime_/totalTime << " %)."<<nl;

    Info<<"   Performed mesh layers in " << setprecision(6)
        << layerTime_ << setprecision(3)
        <<" sec (" << 100.0*layerTime_/totalTime << " %)."<<nl;
    Info<<endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexReport::hexReport
(
    const meshRefinement& meshRefiner
)
:
    meshRefiner_(meshRefiner),
    initialisationTime_(0),
    refinementTime_(0),
    dualTime_(0),
    snappingTime_(0),
    layerTime_(0),
    numLayerCells_(0),
    layerCoverage_(0),
    patchLayerInfo_(0)
{}

// ************************************************************************* //
