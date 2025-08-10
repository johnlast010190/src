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
    (c) 2015-2016 OpenCFD Ltd.

Application
    surfaceRedistributePar

Group
    grpSurfaceUtilities

Description
    (Re)distribution of triSurface. Either takes an undecomposed surface
    or an already decomposed surface and redistributes it so that each
    processor has all triangles that overlap its mesh.

Note
    - best decomposition option is hierarchGeomDecomp since it
      guarantees square decompositions.
    - triangles might be present on multiple processors.
    - merging uses geometric tolerance so take care with writing precision.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "distributedTriSurfaceMesh/distributedTriSurfaceMesh.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistribute.H"
#include "db/IOobjects/IOdictionary/localIOdictionary.H"
#include "decompositionModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Print on master all the per-processor surface stats.
void writeProcStats
(
    const triSurface& s,
    const List<List<treeBoundBox>>& meshBb
)
{
    // Determine surface bounding boxes, faces, points
    List<treeBoundBox> surfBb(Pstream::nProcs());
    {
        surfBb[Pstream::myProcNo()] = treeBoundBox(s.points());
        Pstream::allGatherList(surfBb);
    }

    labelList nPoints(Pstream::nProcs());
    nPoints[Pstream::myProcNo()] = s.points().size();
    Pstream::allGatherList(nPoints);

    labelList nFaces(Pstream::nProcs());
    nFaces[Pstream::myProcNo()] = s.size();
    Pstream::allGatherList(nFaces);

    forAll(surfBb, proci)
    {
        Info<< "processor" << proci << nl;

        const List<treeBoundBox>& bbs = meshBb[proci];
        if (bbs.size())
        {
            Info<< "\tMesh bounds          : " << bbs[0] << nl;
            for (label i = 1; i < bbs.size(); i++)
            {
                Info<< "\t                       " << bbs[i]<< nl;
            }
        }
        Info<< "\tSurface bounding box : " << surfBb[proci] << nl
            << "\tTriangles            : " << nFaces[proci] << nl
            << "\tVertices             : " << nPoints[proci]
            << endl;
    }
    Info<< endl;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Redistribute a triSurface. "
        "The specified surface must be located in the constant/triSurface directory"
    );

    argList::validArgs.append("triSurfaceMesh");
    argList::validArgs.append("distributionType");
    argList::addBoolOption
    (
        "keepNonMapped",
        "Preserve surface outside of mesh bounds"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    runTime.functionObjects().off();

    const fileName surfFileName = args[1];
    const word distTypeName = args[2];
    const label distType =
        distributedTriSurfaceMesh::distributionTypeNames_[distTypeName];

    Info<< "Reading surface from " << surfFileName << nl << nl
        << "Using distribution method "
        << distTypeName << nl << endl;

    const bool keepNonMapped = args.optionFound("keepNonMapped");

    if (keepNonMapped)
    {
        Info<< "Preserving surface outside of mesh bounds." << nl << endl;
    }
    else
    {
        Info<< "Removing surface outside of mesh bounds." << nl << endl;
    }


    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Please run this program on the decomposed case."
            << " It will read surface " << surfFileName
            << " and decompose it such that it overlaps the mesh bounding box."
            << exit(FatalError);
    }


    Random rndGen(653213);

    // For independent decomposition, ensure that distributedTriSurfaceMesh
    // can find the alternative decomposeParDict specified via the
    // -decomposeParDict option.
    if (distType == distributedTriSurfaceMesh::INDEPENDENT)
    {
        fileName decompDictFile;
        args.optionReadIfPresent("decomposeParDict", decompDictFile);

        // A demand-driven decompositionMethod can have issues finding
        // an alternative decomposeParDict location.

        IOdictionary* dictPtr = new IOdictionary
        (
            decompositionModel::selectIO
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                decompDictFile
            )
        );

        // Store it on the object registry, but to be found it must also
        // have the expected "decomposeParDict" name.

        dictPtr->rename("decomposeParDict");
        runTime.store(dictPtr);
    }

    // Determine mesh bounding boxes:
    List<List<treeBoundBox>> meshBb(Pstream::nProcs());
    if (distType == distributedTriSurfaceMesh::FOLLOW)
    {
        #include "include/createPolyMesh.H"

        meshBb[Pstream::myProcNo()] = List<treeBoundBox>
        (
            1,
            treeBoundBox
            (
                boundBox(mesh.points(), false)
            ).extend(rndGen, 1e-3)
        );
        Pstream::allGatherList(meshBb);
    }

    IOobject io
    (
        surfFileName,         // name
        //runTime.findInstance("triSurface", surfFileName),   // instance
        runTime.constant(),   // instance
        "triSurface",         // local
        runTime,              // registry
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    // Look for file (using searchableSurface rules)
    const fileName actualPath(typeFilePath<searchableSurface>(io));
    fileName localPath(actualPath);
    localPath.replace(runTime.rootPath() + '/', "");


    autoPtr<distributedTriSurfaceMesh> surfMeshPtr;

    if (actualPath == io.objectPath())
    {
        Info<< "Loading local (decomposed) surface " << localPath << nl <<endl;
        surfMeshPtr.reset(new distributedTriSurfaceMesh(io));
    }
    else
    {
        Info<< "Loading undecomposed surface " << localPath
            << " on master only" << endl;

        triSurface s;
        List<treeBoundBox> bbs;
        if (Pstream::master())
        {
            // Actually load the surface
            triSurfaceMesh surf(io);
            s = surf;
            bbs = List<treeBoundBox>(1, treeBoundBox(boundBox::greatBox));
        }
        else
        {
            bbs = List<treeBoundBox>(1, treeBoundBox(boundBox::invertedBox));
        }

        dictionary dict;
        dict.add("distributionType", distTypeName);
        dict.add("mergeDistance", SMALL);
        dict.add("bounds", bbs);

        // Scatter patch information
        Pstream::scatter(s.patches());

        // Construct distributedTrisurfaceMesh from components
        IOobject notReadIO(io);
        notReadIO.readOpt() = IOobject::NO_READ;
        surfMeshPtr.reset(new distributedTriSurfaceMesh(notReadIO, s, dict));
    }

    distributedTriSurfaceMesh& surfMesh = surfMeshPtr();


    // Write per-processor stats
    Info<< "Before redistribution:" << endl;
    writeProcStats(surfMesh, meshBb);


    // Do redistribution
    Info<< "Redistributing surface" << nl << endl;
    autoPtr<mapDistribute> faceMap;
    autoPtr<mapDistribute> pointMap;
    surfMesh.distribute
    (
        meshBb[Pstream::myProcNo()],
        keepNonMapped,
        faceMap,
        pointMap
    );
    faceMap.clear();
    pointMap.clear();

    Info<< endl;


    // Write per-processor stats
    Info<< "After redistribution:" << endl;
    writeProcStats(surfMesh, meshBb);


    Info<< "Writing surface." << nl << endl;
    surfMesh.objectRegistry::write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
