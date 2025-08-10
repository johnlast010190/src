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
    (c) 2015 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "primitives/random/Random/Random.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "triSurface/fields/triSurfaceFields.H"
#include "db/Time/Time.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "primitives/Pair/Pair.H"
#include "primitives/quaternion/quaternion.H"
#include "global/constants/mathematical/mathematicalConstants.H"

#ifdef FOAM_USE_TBB
  #include <tbb/parallel_for.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMesh, 0);
    addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);
    word triSurfaceMesh::meshSubDir = "triSurface";
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::triSurfaceMesh::checkFile
(
    const regIOobject& io,
    const bool isGlobal
)
{
    const fileName fName
    (
        isGlobal
      ? io.globalFilePath(typeName)
      : io.localFilePath(typeName)
    );
    if (fName.empty())
    {
        FatalErrorInFunction
            << "Cannot find triSurfaceMesh starting from "
            << io.objectPath() << exit(FatalError);
    }

    return fName;
}


Foam::fileName Foam::triSurfaceMesh::checkFile
(
    const regIOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
{
    fileName fName;
    if (dict.readIfPresent("file", fName, false, false))
    {
        fName.expand();
        if (!fName.isAbsolute())
        {
            fName = io.objectPath().path()/fName;
        }
        if (!exists(fName))
        {
            FatalErrorInFunction
                << "Cannot find triSurfaceMesh at " << fName
                << exit(FatalError);
        }
    }
    else
    {
        fName =
        (
            isGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        if (!exists(fName))
        {
            FatalErrorInFunction
                << "Cannot find triSurfaceMesh starting from "
                << io.objectPath() << exit(FatalError);
        }
    }

    return fName;
}


bool Foam::triSurfaceMesh::addFaceToEdge
(
    const edge& e,
    EdgeMap<label>& facesPerEdge
)
{
    EdgeMap<label>::iterator eFnd = facesPerEdge.find(e);
    if (eFnd != facesPerEdge.end())
    {
        if (eFnd() == 2)
        {
            return false;
        }
        eFnd()++;
    }
    else
    {
        facesPerEdge.insert(e, 1);
    }
    return true;
}


bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::isSurfaceClosed:"
            << " determining closedness for surface with "
            << triSurface::size() << " triangles and "
            << triSurface::points().size() << " points" << endl;
    }

    const pointField& pts = triSurface::points();

    // Construct pointFaces. Let's hope surface has compact point
    // numbering ...
    labelListList pointFaces;
    invertManyToMany(pts.size(), *this, pointFaces);

    // Loop over all faces surrounding point. Count edges emanating from point.
    // Every edge should be used by two faces exactly.
    // To prevent doing work twice per edge only look at edges to higher
    // point
    EdgeMap<label> facesPerEdge(100);
    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        facesPerEdge.clear();
        forAll(pFaces, i)
        {
            const triSurface::FaceType& f = triSurface::operator[](pFaces[i]);
            label fp = findIndex(f, pointi);

            // Something weird: if I expand the code of addFaceToEdge in both
            // below instances it gives a segmentation violation on some
            // surfaces. Compiler (4.3.2) problem?


            // Forward edge
            label nextPointi = f[f.fcIndex(fp)];

            if (nextPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, nextPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    return false;
                }
            }
            // Reverse edge
            label prevPointi = f[f.rcIndex(fp)];

            if (prevPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, prevPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    return false;
                }
            }
        }

        // Check for any edges used only once.
        forAllConstIter(EdgeMap<label>, facesPerEdge, iter)
        {
            if (iter() != 2)
            {
                return false;
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io),
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(s),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts, false);
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        io,
        checkFile(static_cast<const searchableSurface&>(*this), true)
    ),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts, false);
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        io,
        checkFile(static_cast<const searchableSurface&>(*this), dict, true),
        dict
    ),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Reading from supplied file name instead of objectPath/filePath
    dict.readIfPresent("file", fName_, false, false);

    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts, false);

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const readAction r)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Check IO flags
    if (io.readOpt() != IOobject::NO_READ)
    {
        const bool searchGlobal(r == localOrGlobal || r == masterOnly);

        const fileName actualFile
        (
            searchGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        if (debug)
        {
            Pout<< "triSurfaceMesh(const IOobject& io) :"
                << " loading surface " << io.objectPath()
                << " local filePath:" << io.localFilePath(typeName)
                << " from:" << actualFile << endl;
        }

        if (searchGlobal && Pstream::parRun())
        {
            // Check where surface was found
            const fileName localFile(io.localFilePath(typeName));

            if (r == masterOnly && (actualFile != localFile))
            {
                // Found undecomposed surface. Load on master only
                if (Pstream::master())
                {
                    triSurface s2(actualFile);
                    triSurface::transfer(s2);
                }
                Pstream::scatter(triSurface::patches());
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
            else
            {
                // Read on all processors
                triSurface s2(actualFile);
                triSurface::transfer(s2);
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
        }
        else
        {
            // Read on all processors
            triSurface s2(actualFile);
            triSurface::transfer(s2);
            if (debug)
            {
                Pout<< "triSurfaceMesh(const IOobject& io) :"
                    << " loaded triangles:" << triSurface::size() << endl;
            }
        }
    }

    const pointField& pts = triSurface::points();
    bounds() = boundBox(pts, false);
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict,
    const readAction r
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    if (io.readOpt() != IOobject::NO_READ)
    {
        const bool searchGlobal(r == localOrGlobal || r == masterOnly);

        fileName actualFile
        (
            searchGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        // Reading from supplied file name instead of objectPath/filePath
        if (dict.readIfPresent("file", fName_, keyType::LITERAL))
        {
            //fName_ = triSurface::relativeFilePath
            //(
            //    static_cast<const searchableSurface&>(*this),
            //    fName_,
            //    searchGlobal
            //);
            actualFile = fName_;
        }

        if (debug)
        {
            Pout<< "triSurfaceMesh(const IOobject& io, const dictionary&) :"
                << " loading surface " << io.objectPath()
                << " local filePath:" << io.localFilePath(typeName)
                << " from:" << actualFile << endl;
        }

        if (searchGlobal && Pstream::parRun())
        {
            // Check where surface was found
            const fileName localFile(io.localFilePath(typeName));

            if (r == masterOnly && (actualFile != localFile))
            {
                // Surface not loaded from processor directories -> undecomposed
                // surface. Load on master only
                if (Pstream::master())
                {
                    //Read master (optional transforms)
                    triSurface s2(io,actualFile,dict);
                    triSurface::transfer(s2);
                }
                Pstream::scatter(triSurface::patches());
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
            else
            {
                // Read on all processors (optional transforms)
                triSurface s2(io,actualFile,dict);
                triSurface::transfer(s2);

                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
        }
        else
        {
            // Read on all processors (optional transforms)
            triSurface s2(io,actualFile,dict);
            triSurface::transfer(s2);

            if (debug)
            {
                Pout<< "triSurfaceMesh(const IOobject& io) :"
                    << " loaded triangles:" << triSurface::size() << endl;
            }
        }
    }

    const pointField& pts = triSurface::points();
    bounds() = boundBox(pts, false);

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
    clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
    outsideVolType_ = volumeType::UNKNOWN;
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::coordinates() const
{
    tmp<pointField> tPts(new pointField(8));
    pointField& pt = tPts.ref();

    // Use copy to calculate face centres so they don't get stored
    pt = PrimitivePatch<SubList<triSurface::FaceType>, const pointField&>
    (
        SubList<triSurface::FaceType>(*this, triSurface::size()),
        triSurface::points()
    ).faceCentres();

    return tPts;
}


void Foam::triSurfaceMesh::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres = coordinates();
    radiusSqr.setSize(size());
    radiusSqr = 0.0;

    const pointField& pts = triSurface::points();

    forAll(*this, facei)
    {
        const labelledTri& f = triSurface::operator[](facei);
        const point& fc = centres[facei];
        forAll(f, fp)
        {
            const point& pt = pts[f[fp]];
            radiusSqr[facei] = max(radiusSqr[facei], Foam::magSqr(fc-pt));
        }
    }

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::points() const
{
    return triSurface::points();
}


bool Foam::triSurfaceMesh::overlaps(const boundBox& bb) const
{
     const indexedOctree<treeDataTriSurface>& octree = tree();

     labelList indices = octree.findBox(treeBoundBox(bb));

     return !indices.empty();
}


void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    outsideVolType_ = volumeType::UNKNOWN;

    // Update local information (instance, event number)
    searchableSurface::instance() = objectRegistry::time().timeName();
    objectRegistry::instance() = searchableSurface::instance();

    label event = getEvent();
    searchableSurface::eventNo() = event;
    objectRegistry::eventNo() = searchableSurface::eventNo();

    // Clear additional addressing
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);

    bounds() = boundBox(triSurface::points(), false);
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::triSurfaceMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        // Boundary edges
        labelList bEdges
        (
            identity
            (
                nEdges()
               -nInternalEdges()
            )
          + nInternalEdges()
        );

        treeBoundBox bb(Zero, Zero);

        if (bEdges.size())
        {
            label nPoints;
            PatchTools::calcBounds
            (
                *this,
                bb,
                nPoints
            );

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.

            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        scalar oldTol = indexedOctree<treeDataEdge>::perturbTol();
        indexedOctree<treeDataEdge>::perturbTol() = tolerance();

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    localPoints(),  // points
                    bEdges          // selected edges
                ),
                bb,                 // bb
                maxTreeDepth(),     // maxLevel
                10,                 // leafsize
                3.0                 // duplicity
            )
        );

        indexedOctree<treeDataEdge>::perturbTol() = oldTol;
    }
    return edgeTree_();
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(patches().size());
        forAll(regions_, regionI)
        {
            regions_[regionI] = patches()[regionI].name();
        }
    }
    return regions_;
}


bool Foam::triSurfaceMesh::hasVolumeType() const
{
    if (surfaceClosed_ == -1)
    {
        if (isSurfaceClosed())
        {
            surfaceClosed_ = 1;
        }
        else
        {
            surfaceClosed_ = 0;
        }
    }

    return surfaceClosed_ == 1;
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    triSurfaceSearch::findNearest(samples, nearestDistSqr, info, threaded);
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    triSurfaceRegionSearch::findNearest
    (
        samples,
        nearestDistSqr,
        regionIndices,
        info
    );
}


void Foam::triSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    triSurfaceSearch::findLine(start, end, info);
}


void Foam::triSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    triSurfaceSearch::findLineAny(start, end, info, threaded);
}


void Foam::triSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    triSurfaceSearch::findLineAll(start, end, info);
}


void Foam::triSurfaceMesh::findSphere
(
    const pointField& centres,
    const scalarField& radiusSqr,
    List<labelList>& hitIndexes,
    labelList& nOverlapChecks
) const
{
    triSurfaceSearch::findSphere(centres,radiusSqr,hitIndexes,nOverlapChecks);
}


void Foam::triSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region,
    const bool threaded /*= false*/
) const
{
    region.setSize(info.size());

    if (threaded && info.size() > 2048)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                info.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, info.size(), grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    if (info[i].hit())
                    {
                        region[i] = triSurface::operator[](info[i].index()).region();
                    }
                    else
                    {
                        region[i] = -1;
                    }
                }
            },
            tbb::simple_partitioner()
        );
#endif
    }
    else
    {
        forAll(info, i)
        {
            if (info[i].hit())
            {
                region[i] = triSurface::operator[](info[i].index()).region();
            }
            else
            {
                region[i] = -1;
            }
        }
    }
}


void Foam::triSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const triSurface& s = *this;
    const pointField& pts = s.points();

    normal.setSize(info.size());

    if (minQuality_ >= 0)
    {
        // Make sure we don't use triangles with low quality since
        // normal is not reliable.

        const labelListList& faceFaces = s.faceFaces();

        forAll(info, i)
        {
            if (info[i].hit())
            {
                label facei = info[i].index();
                normal[i] = s[facei].areaNormal(pts);

                scalar qual = s[facei].tri(pts).quality();

                if (qual < minQuality_)
                {
                    // Search neighbouring triangles
                    const labelList& fFaces = faceFaces[facei];

                    forAll(fFaces, j)
                    {
                        label nbrI = fFaces[j];
                        scalar nbrQual = s[nbrI].tri(pts).quality();
                        if (nbrQual > qual)
                        {
                            qual = nbrQual;
                            normal[i] = s[nbrI].areaNormal(pts);
                        }
                    }
                }

                normal[i] /= mag(normal[i]) + VSMALL;
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
    else
    {
        forAll(info, i)
        {
            if (info[i].hit())
            {
                label facei = info[i].index();
                // Cached:
                //normal[i] = faceNormals()[facei];

                // Uncached
                normal[i] = s[facei].unitNormal(pts);
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
}


void Foam::triSurfaceMesh::setField(const labelList& values)
{
    if (foundObject<triSurfaceLabelField>("values"))
    {
        triSurfaceLabelField& fld = const_cast<triSurfaceLabelField&>
        (
            lookupObject<triSurfaceLabelField>
            (
                "values"
            )
        );
        fld.field() = values;
    }
    else
    {
        autoPtr<triSurfaceLabelField> fldPtr
        (
            new triSurfaceLabelField
            (
                IOobject
                (
                    "values",
                    objectRegistry::time().timeName(),  // instance
                    meshSubDir,                         // local
                    *this,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                *this,
                dimless,
                labelField(values)
            )
        );

        // Store field on triMesh
        fldPtr.ptr()->store();
    }
}


void Foam::triSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values,
    const bool threaded /*= false*/
) const
{
    if (foundObject<triSurfaceLabelField>("values"))
    {
        values.setSize(info.size());

        const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
        (
            "values"
        );

        if (threaded && info.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    info.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, info.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        if (info[i].hit())
                        {
                            values[i] = fld[info[i].index()];
                        }
                    }
                },
                tbb::simple_partitioner()
            );
#endif
        }
        else
        {
            forAll(info, i)
            {
                if (info[i].hit())
                {
                    values[i] = fld[info[i].index()];
                }
            }
        }
    }
}


void Foam::triSurfaceMesh::clearCurvaturePointField()
{
    if (foundObject<triSurfacePointScalarField>("curvatureValues"))
    {
        const_cast<triSurfacePointScalarField&>
        (
            lookupObject<triSurfacePointScalarField>
            (
                "curvatureValues"
            )
        ).checkOut();
    }
}


void Foam::triSurfaceMesh::setCurvaturePointField(const scalarField& curvature)
{
    if (foundObject<triSurfacePointScalarField>("curvatureValues"))
    {
        triSurfacePointScalarField& fld = const_cast<triSurfacePointScalarField&>
        (
            lookupObject<triSurfacePointScalarField>
            (
                "curvatureValues"
            )
        );
        fld.field() = curvature;
    }
    else
    {
        autoPtr<triSurfacePointScalarField> fldPtr
        (
            new triSurfacePointScalarField
            (
                IOobject
                (
                    "curvatureValues",
                    objectRegistry::time().timeName(),  // instance
                    meshSubDir,                         // local
                    *this,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                *this,
                dimless,
                scalarField(curvature)
            )
        );

        // Store field on triMesh
        fldPtr.ptr()->store();
    }
}


void Foam::triSurfaceMesh::getCurvatureField
(
    const List<pointIndexHit>& info,
    scalarField& curvature
) const
{
    if (foundObject<triSurfacePointScalarField>("curvatureValues"))
    {
        curvature.setSize(info.size());

        const triSurfacePointScalarField& fld =
            lookupObject<triSurfacePointScalarField>
            (
                "curvatureValues"
            );

        const triSurface& s = *this;
        const pointField& pts = s.points();
        forAll(info, i)
        {
            if (info[i].hit())
            {
                const labelledTri& f =  s[info[i].index()];
                point hPt = info[i].hitPoint();
                scalar fArea = mag((pts[f[1]]-pts[f[0]])^(pts[f[2]]-pts[f[0]]));
                scalar weight0 = mag((pts[f[1]]-hPt)^(pts[f[2]]-hPt))/fArea;
                scalar weight1 = mag((pts[f[0]]-hPt)^(pts[f[2]]-hPt))/fArea;
                scalar weight2 = mag((pts[f[0]]-hPt)^(pts[f[1]]-hPt))/fArea;

                curvature[i] = weight0*fld[f[0]] + weight1*fld[f[1]]
                    +weight2*fld[f[2]];
            }
            else
            {
                curvature[i] = -1;
            }
        }
    }
}


void Foam::triSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType,
    const bool threaded /*= false*/
) const
{
    // Precalculate and cache value for outside points
    if (outsideVolType_ == volumeType::UNKNOWN && hasVolumeType())
    {
        point midBbPt = 0.5*(bounds().max() + bounds().min());
        scalar dr = bounds().mag();
        pointField outsidePoints(3);
        outsidePoints[0] = midBbPt + point(dr, 0., 0.);
        outsidePoints[1] = midBbPt + point(0., dr, 0.);
        outsidePoints[2] = midBbPt + point(0., 0., dr);

        label nInside = 0;
        label nOutside = 0;
        forAll(outsidePoints, pointi)
        {
            point pt = outsidePoints[pointi];
            if (!tree().bb().contains(pt))
            {
                volumeType vtype = tree().shapes().getVolumeType
                (
                    tree(),
                    pt
                );
                if(vtype == volumeType::INSIDE)
                {
                    nInside++;
                }
                else if(vtype == volumeType::OUTSIDE)
                {
                    nOutside++;
                }
            }
        }
        if (nInside > nOutside)
        {
            outsideVolType_ = volumeType::INSIDE;
        }
        else if(nInside < nOutside)
        {
            outsideVolType_ = volumeType::OUTSIDE;
        }
    }

    volType.setSize(points.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    if (threaded && points.size() > 2048)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                points.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, points.size(), grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t pointi = r.begin(); pointi < r.end(); ++pointi)
                {
                    const point& pt = points[pointi];

                    if (!tree().bb().contains(pt))
                    {
                        if (hasVolumeType())
                        {
                            // Precalculate and cache value for this outside point
                            if (outsideVolType_ == volumeType::UNKNOWN)
                            {
                                outsideVolType_ = tree().shapes().getVolumeType(tree(), pt);
                            }
                            volType[pointi] = outsideVolType_;
                        }
                        else
                        {
                            // Have to calculate directly as outside the octree
                            volType[pointi] = tree().shapes().getVolumeType(tree(), pt);
                        }
                    }
                    else
                    {
                        // - use cached volume type per each tree node
                        volType[pointi] = tree().getVolumeType(pt);
                    }
                }
            },
            tbb::simple_partitioner()
        );
#endif
    }
    else
    {
        forAll(points, pointi)
        {
            const point& pt = points[pointi];

            if (!tree().bb().contains(pt))
            {
                if (hasVolumeType())
                {
                    // Precalculate and cache value for this outside point
                    if (outsideVolType_ == volumeType::UNKNOWN)
                    {
                        outsideVolType_ = tree().shapes().getVolumeType(tree(), pt);
                    }
                    volType[pointi] = outsideVolType_;
                }
                else
                {
                    // Have to calculate directly as outside the octree
                    volType[pointi] = tree().shapes().getVolumeType(tree(), pt);
                }
            }
            else
            {
                // - use cached volume type per each tree node
                volType[pointi] = tree().getVolumeType(pt);
            }
        }
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


bool Foam::triSurfaceMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    const Time& runTime = searchableSurface::time();
    const fileName& instance = searchableSurface::instance();

    if
    (
        instance != runTime.timeName()
     && instance != runTime.system()
     && instance != runTime.caseSystem()
     && instance != runTime.constant()
     && instance != runTime.caseConstant()
    )
    {
        const_cast<triSurfaceMesh&>(*this).searchableSurface::instance() =
            runTime.timeName();
        const_cast<triSurfaceMesh&>(*this).objectRegistry::instance() =
            runTime.timeName();
    }

    fileName fullPath;
    if (fName_.size())
    {
        // Override file name

        fullPath = fName_;

        fullPath.expand();
        if (!fullPath.isAbsolute())
        {
            // Add directory from regIOobject
            fullPath = searchableSurface::objectPath().path()/fullPath;
        }
    }
    else
    {
        fullPath = searchableSurface::objectPath();
    }

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
