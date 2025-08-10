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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cad/CADReader.H"
#include "distributedCADSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "hierarchGeomDecomp/hierarchGeomDecomp.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedCADSurface, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        distributedCADSurface,
        dict
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distributedCADSurface::decomposeBoundBox
(
    const boundBox& bounds
)
{
    // Get bb of all domains.
    procBb_.setSize(Pstream::nProcs());

    hierarchGeomDecomp decomp
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                searchableSurface::time().system(),
                searchableSurface::time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    const Vector<label>& n = decomp.n();

    if (n.x()*n.y()*n.z() != Pstream::nProcs())
    {
        FatalErrorInFunction
            << "Number of divisions " << n
            << " does not match number of processors: "
            << Pstream::nProcs()
            << abort(FatalError);
    }

    scalar deltax = bounds.span().x() / scalar(n.x());
    scalar deltay = bounds.span().y() / scalar(n.y());
    scalar deltaz = bounds.span().z() / scalar(n.z());

    label proci = 0;

    const scalar& bbXmin = bounds.min().x();
    const scalar& bbYmin = bounds.min().y();
    const scalar& bbZmin = bounds.min().z();

    for (label xi = 0; xi < n.x(); xi++)
    {
        scalar xmin = bbXmin + xi*deltax;
        scalar xmax = xmin + deltax;

        for (label yi = 0; yi < n.y(); yi++)
        {
            scalar ymin = bbYmin + yi*deltay;
            scalar ymax = ymin + deltay;

            for (label zi = 0; zi < n.z(); zi++)
            {
                scalar zmin = bbZmin + zi*deltaz;
                scalar zmax = zmin + deltaz;

                boundBox bb
                (
                    point(xmin, ymin, zmin),
                    point(xmax, ymax, zmax)
                );

                bb.inflate(0.001);

                if (debug)
                {
                    Info<<"\nprocessor" << proci <<endl;
                    Info<<"bounds " << bb <<endl;
                    Info<<"span " << bb.span() <<endl;
                    Info<<"midpoint " << bb.midpoint() <<endl;
                }

                procBb_[proci] =
                    treeBoundBoxList(1, treeBoundBox(bb));

                proci++;
            }
        }
    }
}

Foam::triSurface Foam::distributedCADSurface::triangulate
(
    CADReader& reader,
    const scalar& linDeflection,
    const scalar& angDeflection,
    const geometricSurfacePatchList& globalPatches
)
{
    pointField pointLst;
    List<labelledTri> faceLst;

    HashTable<label> patchIDs;

    // fill patchIDs in case globalPatches is not empty
    forAll(globalPatches, patchi)
    {
        patchIDs.insert
        (
            globalPatches[patchi].name(),
            globalPatches[patchi].index()
        );
    }

    reader.triangulate
    (
        linDeflection,
        angDeflection,
        pointLst,
        faceLst,
        patchIDs
    );

    geometricSurfacePatchList patches(patchIDs.size());

    forAllConstIters(patchIDs, iter)
    {
        const label patchIdx = iter.object();

        patches[patchIdx] =
            geometricSurfacePatch
            (
               iter.key(),
               patchIdx
            );
    }

    triSurface s
    (
        faceLst,
        patches,
        pointLst
    );

    return s;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedCADSurface::distributedCADSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    // if parallel run:
    //   each processor computes very coarse triangulation of CAD,
    //   so it knows the overall boundBox and patch names
    // if serial run triangulates with wanted accuracy
    triSurfaceMesh
    (
        io,
        (
            Pstream::parRun()
          ? triSurface()
          : triSurface
            (
                io.objectPath(),
                dict.lookupOrDefault<scalar>("linDeflection", 0.01),
                dict.lookupOrDefault<scalar>("angDeflection", 0.5),
                dict.lookupOrDefault<bool>("forceTriangulation", false),
                dict.lookupOrDefault<word>("nameRegionsBy", "solid")
            )
        )
    )
{
    if (!Pstream::parRun())
    {
        writeStats(Info);
        return;
    }

    fileName inputFile = io.globalFilePath(io.objectPath());

    const scalar linDeflection = dict.lookupOrDefault<scalar>("linDeflection", 0.01);
    const scalar angDeflection = dict.lookupOrDefault<scalar>("angDeflection", 0.5);
    const bool forceTriangulation = dict.lookupOrDefault<bool>("forceTriangulation", false);

    fileName triFile =
        inputFile.lessExt()
      + "_OCCT_" + Foam::name(linDeflection)
      + "_" + Foam::name(angDeflection)
      + "_proc"
      + Foam::name(Pstream::myProcNo())
      + ".obj";

    if (isFile(triFile) && !forceTriangulation)
    {
        Info<< "Skipping triangulation."
            << nl << "Found " << triFile
            << endl;

        triSurface s(triFile);
        triSurface::transfer(s);
        return;
    }

    //TODO currently all procs read the CAD file.
    // consider master only reads and send to other procs
    CADReader reader
    (
        inputFile,
        dict.lookupOrDefault<word>("nameRegionsBy", "solid"),
        true //distributed
    );

    // coarse triangulation to get all regions
    triSurface coarseSurf =
        triangulate
        (
            reader,
            10, //linDeflection,
            1   //angDeflection,
        );

    const word decompMethod = dict.lookupOrDefault<word>("decompositionMethod", "hierarchical");

    if (decompMethod == "hierarchical")
    {
        Info<<"\nCAD decomposition method: hierarchical" << endl;

        decomposeBoundBox
        (
            boundBox(coarseSurf.points(), false)
        );

        // OpenCASCADE triangulate with wanted accuracy only inside processor boundBox
        reader.distributeShapes
        (
            procBb_
        );
    }
    else if (decompMethod == "balanceArea")
    {
        Info<<"\nCAD decomposition method: balanceArea" << endl;

        // TODO redistribute with refinement level weights
        reader.distributeShapes(/*balanceArea*/true);
    }
    else if (decompMethod == "balanceTriangles")
    {
        Info<<"\nCAD decomposition method: balanceTriangles" << endl;

        // Each processor gets roughly same number of coarse triangles
        reader.distributeShapes();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown CAD decompositionMethod "
            << decompMethod
            << ". Available methods are hierarchical, balanceArea, balanceTriangles."
            << abort(FatalError);
    }

    triSurface s =
        triangulate
        (
            reader,
            linDeflection,
            angDeflection,
            coarseSurf.patches()
        );

    if (debug && Pstream::master())
    {
        word outNameCoarse = inputFile.nameLessExt() + std::string("_coarse.obj");

        fileName outFileCoarse{inputFile.path(), outNameCoarse};

        Info<< "\nWriting " << outFileCoarse <<endl;

        coarseSurf.write(outFileCoarse, /*sortByRegion*/ true, /*keepEmpty*/ true);
    }

    Pout<< "Writing distributed OCCT triangulation to "
        << triFile
        << endl;

    s.write(triFile, /*sortByRegion*/ true, /*keepEmpty*/ true);
    triSurface::transfer(s);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributedCADSurface::~distributedCADSurface()
{}


// ************************************************************************* //
