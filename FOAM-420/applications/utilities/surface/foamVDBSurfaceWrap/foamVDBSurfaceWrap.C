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

//OpenVDB
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridTransformer.h> //for Interpolation / GridSampler
#include <openvdb/tree/LeafManager.h> // used by sharpenFeatures
#include <openvdb/tools/Morphology.h> // dilateVoxels, erodeVoxels

//OpenFOAM
#include "foamVDB.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/argList/argList.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "searchableSurfaces/CADSurface/CADSurface.H"
#include "memInfo/memInfo.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "triSurface/triSurfaceTools/triSurfaceTools.H"

#include "geometryUtil.H"

using namespace Foam;

// * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();

    argList::noCheckProcessorDirectories();

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    memInfo mem;
    Info<< "Initial memory " << mem.update().size() << " kB" << endl;

    const word dictName("foamVDBDict");

    #include "include/setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;
    IOdictionary meshDict(dictIO);

    // all surface geometry
    dictionary& geometryDict = meshDict.subDict("geometry");

    bool debug = meshDict.lookupOrDefault<bool>("debug", false);

    const label memoryLimit = meshDict.lookupOrDefault<label>("memoryLimit", 1000000);

    const bool sharpenFeatures = meshDict.lookupOrDefault<bool>("sharpenFeatures", false);
    const scalar adaptivity= meshDict.lookupOrDefault<scalar>("adaptivity", 1e-3);

    // featureAngle: if ((faceNormal & neighbourFaceNormal) < edgeTolerance) -> is edgeVoxel
    const scalar edgetolerance = meshDict.lookupOrDefault<scalar>("edgetolerance", 0.5);

    const fileName outputFile = meshDict.lookup("outputFile");

    const scalar iso = meshDict.lookupOrDefault<scalar>("isoValue", 0);


    // Read geometry
    // ~~~~~~~~~~~~~
    autoPtr<searchableSurfaces> allGeometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                runTime.time().constant(),  // directory
                "triSurface",               // instance
                runTime.time(),             // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            geometryDict
        )
    );

    Info<< "Read geometry in = "
        << runTime.time().cpuTimeIncrement() << " s" << nl << endl;

    scalar minEdge = readScalar(meshDict.lookup("minEdge"));

    foamVDB ovdb(minEdge, /*maxLevel*/0);

    List<triSurface> surfList(allGeometryPtr().size());
    wordList allRegNames;
    label nRegions = 0;

    forAll(allGeometryPtr(), geomi)
    {
        const wordList& regNames = allGeometryPtr().regionNames()[geomi];

        nRegions += regNames.size();

        forAll(regNames, i)
        {
            allRegNames.append(regNames[i]);
        }

        const searchableSurface& s = allGeometryPtr()[geomi];

        if (isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(s);
            triSurface& inputSurf = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
            );

            surfList[geomi] = inputSurf;
        }
    }

    triSurface globalSurf = ovdb.addSurfaces(surfList);

    if (debug && isA<CADSurface>(allGeometryPtr()[0]))
    {
        //if input is stp or iges
        word outNameOcct= outputFile.nameLessExt() + std::string("_occt.obj");
        fileName outFileOcct{outputFile.path(), outNameOcct};
        Info<< "\nWriting OCCT triangulation to " << outFileOcct << endl;
        globalSurf.write(outFileOcct, /*sortByRegion*/true);
    }

    FloatGrid::Ptr surfLevelSet;

    {
        Timer timer("Calculate level-set");

        scalar halfWidth = max(3, (iso/minEdge) + 1);

        /*FloatGrid::Ptr*/ surfLevelSet =
            ovdb.meshToLevelSet
            (
                globalSurf,
                /*cellLevel*/0,
                halfWidth
            );

        Info<<"Active voxels: " << surfLevelSet->activeVoxelCount() <<endl;

        if (ovdb.isOpenSurface(surfLevelSet))
        {
            Info<< "WARNING: surface is NOT closed"
                << endl;
        }
    }

    #include "sharpenFeatures.H"

    Info<< "\nWriting " << outputFile << endl;
    wrapSurf.write(outputFile, /*sortByRegion*/true);
    wrapSurf.writeStats(Info);

    Info<< "\nelapsedCpuTime: " << runTime.elapsedCpuTime() << "s\n"
        << "elapsedClockTime: " << runTime.elapsedClockTime() << "s, ";

    scalar memoryPeak = mem.update().peak();

    if (memoryPeak > 1e6)
    {
        Info<< memoryPeak * 1e-6 << " Gb (peak)" << nl << endl;
    }
    else if (memoryPeak > 1e3)
    {
        Info<< memoryPeak * 1e-3 << " Mb (peak)" << nl << endl;
    }
    else
    {
        Info<< memoryPeak << " kB (peak)" << nl << endl;
    }

return 0;
}

// ************************************************************************* //
