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
    (c) 2010-2023 Esi Ltd.

Application
    frontalArea

Description
    Calculates the frontal area on a specified geometry

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "searchableSurfaces/searchableSphere/searchableSphere.H"
#include "searchableSurfaces/searchableBox/searchableBox.H"
#include "searchableSurfaces/searchableCylinder/searchableCylinder.H"
#include "searchableSurfaces/searchableRing/searchableRing.H"
#include "searchableSurfaces/searchableCone/searchableCone.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "coordinate/systems/coordinateSystem.H"
#include "meshes/meshTools/mergePoints.H"
#include "delabella.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

dictionary getCoordinateDict
(
    const Time& runTime,
    const dictionary& dict
)
{
    if (dict.found("coordinates"))
    {
        return dict.subDict("coordinates");
    }
    else if (dict.found("referenceFrame"))
    {
        const word frameName = dict.lookup("referenceFrame");
        IOdictionary meshObjects
        (
            IOobject
            (
                "meshObjects",
                runTime.time().system(),
                runTime.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (!meshObjects.found(frameName))
        {
            FatalErrorInFunction
                << "specified referenceFrame " << frameName
                << " not found in system/meshObjects"
                << abort(FatalError);
        }

        return meshObjects.subDict(frameName).subDict("coordinateSystem");
    }
    else
    {
        dictionary defaultDict("coordinates");
        defaultDict.add("origin",vector::zero,true);
        defaultDict.add(word("coordinateRotation"), dictionary(), false);
        dictionary& rotDict = defaultDict.subDict("coordinateRotation");
        rotDict.add("type","axesRotation",true);
        rotDict.add("e1",vector(1,0,0),true);
        rotDict.add("e2",vector(0,1,0),true);

        return defaultDict;
    }
}


boundBox getBoundingBox
(
    const searchableSurfaces& allGeometry,
    const coordinateSystem& coord
)
{
    point minPt = point(GREAT, GREAT, GREAT);
    point maxPt = point(-GREAT, -GREAT, -GREAT);

    forAll(allGeometry, surfI)
    {
        const searchableSurface& s = allGeometry[surfI];
        if (isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& triMesh =
                refCast<const triSurfaceMesh>(s);
            triSurface& tri = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
            );
            const pointField& pts = tri.points();

            forAll(pts, pti)
            {
                point localPt = coord.invTransformPoint(pts[pti]);
                minPt = min(minPt, localPt);
                maxPt = max(maxPt, localPt);
            }
        }
        else if
        (
            isA<searchableSphere>(s) || isA<searchableBox>(s)
         || isA<searchableCylinder>(s) || isA<searchableRing>(s)
         || isA<searchableCone>(s)
        )
        {
            const boundBox& surfBb = s.bounds();
            const tmp<pointField> tpts = surfBb.points();
            const pointField& pts = tpts();
            forAll(pts, pti)
            {
                point localPt = coord.invTransformPoint(pts[pti]);
                minPt = min(minPt, localPt);
                maxPt = max(maxPt, localPt);
            }
        }
    }

    reduce
    (
       std::tie(minPt, maxPt),
       ParallelOp<minMagSqrOp<vector>,
       maxMagSqrOp<vector>>{}
    );

    return boundBox(minPt,maxPt);
}


void filterFaces
(
   const autoPtr<searchableSurfaces>& constraintSurfaces,
   const boundBox& constraintBB,
   const pointField& filterPts,
   const labelList& filterIndex,
   labelList& markedIndex
)
{
   Info<<"Filter outside faces"<<endl;

   if (constraintSurfaces.valid())
   {
      pointField endPoints(filterPts.size());
      forAll(constraintSurfaces(), surfi)
      {
         const searchableSurface& s =
            constraintSurfaces()[surfi];
         scalar dr = 0;
         point midBbPt = vector::zero;
         if
         (
            isA<triSurfaceMesh>(s) || isA<searchableSphere>(s)
            || isA<searchableBox>(s) || isA<searchableCylinder>(s)
            || isA<searchableRing>(s) || isA<searchableCone>(s)
         )
         {
            if (isA<triSurfaceMesh>(s) && !s.hasVolumeType())
            {
               WarningInFunction
               << "The following tri-surface is non-manifold which "
               << "might produce isssues for filter insideness check : "
               << s.name() <<endl;
            }
            const boundBox& surfBb = s.bounds();
            midBbPt = surfBb.midpoint();
            dr = 2*surfBb.mag();
         }
         else
         {
            WarningInFunction
               << "Ignoring filter surface as not of closed type  : "
               << s.name() <<endl;
            continue;
         }

         labelList nOut(filterPts.size(), 0);
         for (direction i = 0; i < vector::nComponents; i++)
         {
            endPoints = vector::zero;
            forAll(filterPts, filteri)
            {
               vector outsidePt = vector::zero;
               outsidePt[i] = dr;
               outsidePt += midBbPt;
               endPoints[filteri] = outsidePt;
            }
            List<List<pointIndexHit>> hitInfo(filterPts.size());
            constraintSurfaces()[surfi].findLineAll
            (
               filterPts,
               endPoints,
               hitInfo
            );

            forAll(hitInfo,hiti)
            {
               label nHits = hitInfo[hiti].size();
               //Filter out coincident hit points
               if (nHits > 1)
               {
                  for (int j = 1;j < nHits;j++)
                  {
                     scalar hitNbrDist = mag
                     (
                        hitInfo[hiti][j].hitPoint()
                        - hitInfo[hiti][j-1].hitPoint()
                     );
                     if (hitNbrDist < SMALL)
                     {
                        nHits--;
                     }
                  }
               }

               if ((nHits % 2) == 0)
               {
                  nOut[hiti]++;
               }
            }
         }

         forAll(nOut, hiti)
         {
            if (nOut[hiti] > 1)
            {
               label index = filterIndex[hiti];
               markedIndex[index] = -1;
            }
         }
      }
   }
   else
   {
      forAll(filterPts, filteri)
      {
         if (!constraintBB.contains(filterPts[filteri]))
         {
            label index = filterIndex[filteri];
            markedIndex[index] = -1;
         }
      }
   }

   return;
}


pointField getMasterPoints
(
    const Time& runTime,
    const dictionary& dict,
    const coordinateSystem& coord,
    const autoPtr<searchableSurfaces>& constraintSurfaces,
    const boundBox& constraintBB,
    point& origin
)
{
    const label coarsen =
        dict.lookupOrDefault<label>("coarsen", 1);

    DynamicList<point> projectedPts(0);
    if (dict.found("geometry"))
    {
        Info<<"Reading frontal area geometry :" <<endl;
        if (Pstream::master())
        {
            const dictionary& geometryDict = dict.subDict("geometry");
            autoPtr<searchableSurfaces> allGeometryPtr
            (
                new searchableSurfaces
                (
                    IOobject
                    (
                        "abc",              // dummy name
                        runTime.time().constant(),     // directory
                        "triSurface",    // instance
                        runTime.time(),                // registry
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    geometryDict
                 )
            );
            const searchableSurfaces& allGeometry = allGeometryPtr();
            label nPts = 0;
            forAll(allGeometry, surfI)
            {
                const searchableSurface& s = allGeometry[surfI];
                if (isA<triSurfaceMesh>(s))
                {
                    const triSurfaceMesh& triMesh =
                        refCast<const triSurfaceMesh>(s);
                    triSurface& tri = const_cast<triSurface&>
                    (
                        refCast<const triSurface>(triMesh)
                    );
                    const pointField& pts = tri.points();
                    nPts += pts.size();
                }
            }

            projectedPts.setSize(nPts);
            forAll(allGeometry, surfI)
            {
                const searchableSurface& s = allGeometry[surfI];
                if (isA<triSurfaceMesh>(s))
                {
                    const triSurfaceMesh& triMesh =
                    refCast<const triSurfaceMesh>(s);
                    triSurface& tri = const_cast<triSurface&>
                    (
                        refCast<const triSurface>(triMesh)
                     );
                    const pointField& pts = tri.points();

                    forAll(pts, i)
                    {
                        if ((i % coarsen) == 0)
                        {
                           const point& pt = pts[i];
                           projectedPts.append(pt);
                        }
                    }
                }
            }
        }
    }
    else if (dict.found("patches"))
    {
        const Foam::polyMesh mesh
        (
            Foam::IOobject
            (
                Foam::polyMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
             )
        );

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        labelHashSet sourceSet = patches.patchSet
        (
            wordReList(dict.lookup("patches")),
            false,
            true
        );
        const labelList& patchIDs = sourceSet.toc();
        label nFaces = 0;
        forAll(patchIDs, i)
        {
            const polyPatch& pp = patches[patchIDs[i]];
            nFaces += pp.size();
        }
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

        autoPtr<indirectPrimitivePatch> ppPtr
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), addressing),
                mesh.points()
             )
        );

        // Gather all unique points on master
        List<pointField> gatheredPoints(Pstream::nProcs());
        gatheredPoints[Pstream::myProcNo()] = ppPtr().localPoints();
        Pstream::gatherList(gatheredPoints);

        if (Pstream::master())
        {
            forAll(gatheredPoints, proci)
            {
                const pointField& procPts = gatheredPoints[proci];
                forAll(procPts, i)
                {
                    const point& pt = procPts[i];
                    projectedPts.append(pt);
                }
            }
            gatheredPoints.clear();
        }
    }

    labelList markedIndex(projectedPts.size(), 0);
    pointField filterPts(projectedPts);
    filterFaces
    (
       constraintSurfaces,
       constraintBB,
       filterPts,
       identity(projectedPts.size()),
       markedIndex
    );

    label nKept = 0;
    forAll(markedIndex, i)
    {
       if (markedIndex[i] > -1)
       {
          filterPts[nKept++] = projectedPts[i];
       }
    }
    filterPts.setSize(nKept);

    boundBox bbox(filterPts);
    if (Pstream::master() && filterPts.size())
    {
        origin = bbox.midpoint();
        const vector& dir = coord.e1();
        plane pl(origin, dir);
        forAll(filterPts, i)
        {
            point pt = pl.nearestPoint(filterPts[i]);
            filterPts[i] = pt;
        }

        //Merge coincdent points to sepcified tolerance
        labelList oldToNew;
        pointField uniquePts;
        mergePoints
        (
            filterPts,
            1e-8,
            false,
            oldToNew,
            uniquePts
        );
        return uniquePts;
    }

    return pointField(0);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates the frontal area"
    );

    argList::noCheckProcessorDirectories();

    argList::addBoolOption
    (
        "delaunay",
        "Use Delaunay triangulation method for calculation"
    );

    argList::addBoolOption
    (
        "writeFrontal",
        "Write frontal are information to file"
    );

#include "include/addDictOption.H"

#include "include/setRootCase.H"
#include "include/createTime.H"

    const word dictName("frontalAreaDict");
#include "include/setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    const IOdictionary dict(dictIO);

    const bool writeVTK =
        dict.lookupOrDefault<bool>("writeVTK", false);

    coordinateSystem coord
    (
        coordinateSystem::New
        (
            "cartesian",
            getCoordinateDict(runTime, dict)
         )
    );

    vector constraintMin(vector::min);
    vector constraintMax(vector::max);
    if (dict.found("boundBox"))
    {
        const dictionary bbDict = dict.subDict("boundBox");
        constraintMin = bbDict.lookup("min");
        constraintMax = bbDict.lookup("max");
    }
    boundBox constraintBB(constraintMin, constraintMax);

    autoPtr<searchableSurfaces> constraintSurfaces;
    if (dict.found("boundSurfaces"))
    {
        Info<<"Reading filter surfaces :" <<endl;
        const dictionary& boundingDict = dict.subDict("boundSurfaces");
        constraintSurfaces.set
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",              // dummy name
                    runTime.time().constant(),     // directory
                    "triSurface",    // instance
                    runTime.time(),                // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                 ),
                boundingDict
             )
        );
    }

    scalar frontalArea = 0.;
    if (args.optionFound("delaunay"))
    {
        Info<<"Using Delaunay frontal area method"<< nl <<endl;

        point origin = vector::zero;
        pointField cloudPts = getMasterPoints
        (
            runTime,
            dict,
            coord,
            constraintSurfaces,
            constraintBB,
            origin
        );

        if (Pstream::master() && cloudPts.size())
        {
            const scalar filterLength =
                dict.lookupOrDefault<scalar>("filterLength", 0.05);
            label nPts = cloudPts.size();
            struct MyPoint
            {
                char something;
                scalar x;
                scalar y;
            };
            MyPoint* cloud = new MyPoint[nPts];
            forAll(cloudPts, i)
            {
                point pt = cloudPts[i];
                vector lvec = pt-origin;
                cloud[i].x = (lvec&coord.e2());
                cloud[i].y = (lvec&coord.e3());
            }

            IDelaBella0* idb = IDelaBella0::Create();

            int verts = idb->Triangulate
            (
                nPts,
                &cloud->x,
                &cloud->y,
                sizeof(MyPoint)
             );

            label tris = verts / 3;
            DynamicList<labelledTri> faces(tris);
            if (verts > 0)
            {
                const DelaBella0_Triangle* dela =
                    idb->GetFirstDelaunayTriangle();
                for (label i = 0; i < tris; i++)
                {
                    label v0 = dela->v[0]->i;
                    label v1 = dela->v[1]->i;
                    label v2 = dela->v[2]->i;
                    point pt0(cloudPts[v0]);
                    point pt1(cloudPts[v1]);
                    point pt2(cloudPts[v2]);

                    if
                    (
                        mag(pt0-pt1) > filterLength
                        || mag(pt1-pt2) > filterLength
                        || mag(pt2-pt0) > filterLength
                    )
                    {
                        dela = dela->next;
                        continue;
                    }
                    else
                    {
                        labelledTri f = labelledTri
                        (
                            dela->v[0]->i,
                            dela->v[1]->i,
                            dela->v[2]->i,
                            0
                         );
                        faces.append(f);
                        dela = dela->next;
                    }
                }
            }

            delete[] cloud;
            idb->Destroy();

            triSurface tri(faces, cloudPts);
            const scalarField& sAreas = tri.magSf();
            forAll(sAreas, facei)
            {
                frontalArea += sAreas[facei];
            }

            if (writeVTK)
            {
                tri.write("delaunayFrontalArea.vtk");
            }
        }
    }
    else
    {
        Info<<"Using shadow frontal area method"<< nl <<endl;
        const scalar resolution =
            dict.lookupOrDefault<scalar>("resolution", 0.0005);
        Info<<"Reading frontal area geometry : " <<endl;
        const dictionary& geometryDict = dict.subDict("geometry");
        autoPtr<searchableSurfaces> allGeometryPtr
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",              // dummy name
                    runTime.time().constant(),     // directory
                    "triSurface",    // instance
                    runTime.time(),                // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                 ),
                geometryDict
             )
        );
        const searchableSurfaces& allGeometry = allGeometryPtr();
        boundBox bbox = getBoundingBox
        (
            allGeometry,
            coord
        );

        scalar dx = bbox.max().x() - bbox.min().x();
        scalar dy = bbox.max().y() - bbox.min().y();
        scalar dz = bbox.max().z() - bbox.min().z();

        scalar startx =  bbox.min().x() - 0.1*dx;
        scalar endx = bbox.max().x() + 0.1*dx;

        point minPt(startx,bbox.min().y(),bbox.min().z());

        label ny = max(dy/resolution,1) + 4;
        label nz = max(dz/resolution,1) + 4;

        label startY = 0;
        label endY = ny;

        if (Pstream::parRun())
        {
            label nPerProcY = ny/Pstream::nProcs();
            if (Pstream::myProcNo() == Pstream::nProcs() -1)
            {
                startY = nPerProcY*Pstream::myProcNo();
                endY = ny;
            }
            else
            {
                startY = nPerProcY*Pstream::myProcNo();
                endY = nPerProcY*(Pstream::myProcNo()+1)+1;
            }
        }

        label nPoints = (endY-startY)*nz;
        label nFaces = ((endY-startY)-1)*(nz-1);

        if (nPoints < 0)
        {
            FatalErrorIn(args.executable())
                <<"resolution in frontalAreaDict needs to be increased"
                << abort(FatalError);
        }
        else if (nPoints > 0.025*labelMax)
        {
            scalar estPts = nPoints/5.0e7;
            WarningIn(args.executable())
                <<"resolution in frontalAreaDict may be too small "
                <<"suggested minimum value: "<<Foam::sqrt(estPts)*resolution
                <<endl;
        }

        pointField points(nPoints);
        label pCount = 0;
        for (int i = startY; i < endY; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                point basePt = bbox.min();
                basePt.y() = basePt.y() + resolution*(i-1);
                basePt.z() = basePt.z() + resolution*(j-1);
                points[pCount++] = coord.transformPoint(basePt);
            }
        }

        List<face> faces(nFaces);
        pCount = 0;
        label fCount = 0;

        for (int i = startY; i < endY-1; i++)
        {
            for (int j = 0; j < nz-1; j++)
            {
                face f0(4);
                f0[0] = pCount;
                f0[1] = pCount+1;
                f0[2] = pCount+nz+1;
                f0[3] = pCount+nz;
                faces[fCount++] = f0;
                pCount++;
            }
            pCount++;
        }

        pointField start(nFaces);

        forAll(faces, faceI)
        {
            const face& f = faces[faceI];
            point centre = f.centre(points);
            start[faceI] = centre;
        }

        pointField end(start+(endx-startx)*coord.e1());

        labelList markedIndex(start.size());
        List<pointIndexHit> hitInfo(start.size());

        Info<<"Calculate intersections"<<endl;

        allGeometry.findAnyIntersection
        (
            start,
            end,
            markedIndex,
            hitInfo
        );

        pointField filterPts(start.size());
        labelList filterIndex(start.size());
        label nHits = 0;
        forAll(hitInfo, hiti)
        {
           if (hitInfo[hiti].hit())
           {
              filterPts[nHits] = hitInfo[hiti].hitPoint();
              filterIndex[nHits] = hiti;
              nHits++;
           }
        }
        filterPts.setSize(nHits);
        filterIndex.setSize(nHits);

        filterFaces
        (
           constraintSurfaces,
           constraintBB,
           filterPts,
           filterIndex,
           markedIndex
        );

        forAll(filterIndex, hiti)
        {
           label index = filterIndex[hiti];
           if (markedIndex[index] > -1)
           {
              const face& f = faces[index];
              scalar area = f.mag(points);
              frontalArea += area;
           }
        }
        reduce(frontalArea, sumOp<scalar>());

        if (writeVTK)
        {
            Info<<"Write intersections"<<endl;
            simpleVTKWriter frontalAreaVTK
            (
                faces,
                points
            );
            frontalAreaVTK.addFaceData("hitSurface", markedIndex);
            frontalAreaVTK.write("rayFrontalArea.vtk");
        }
    }

    if (args.optionFound("writeFrontal"))
    {
        fileName baseDir = runTime.time().path();
        if (Pstream::parRun())
        {
           baseDir = baseDir/".."/"postProcessing"/"frontalArea";
        }
        else
        {
           baseDir = baseDir/"postProcessing"/"frontalArea";
        }
        if (!exists(baseDir))
        {
           mkDir(baseDir);
        }
        const fileName fAreaFile
        (
           args.optionFound("delaunay")
           ? "delaunayFrontalArea.dat" : "frontalArea.dat"
        );
        OFstream ofs = OFstream(baseDir/fAreaFile);
        ofs << setw(1) << "#" << setw(1) << ' '
            << setf(ios_base::left) << "frontalArea" <<nl;
        ofs << frontalArea << nl;
    }

    Info<< nl << "Frontal area : " << frontalArea <<endl;

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
