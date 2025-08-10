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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Application
    surfacePatch

Group
    grpSurfaceUtilities

Description
    Patches (regionises) a surface using a user-selectable method.

    The current methods are either based on a geometric feature angle
    (a replacement for the surfaceAutoPatch functionality) or on intersecting
    a set of searchableSurfaces - this will re-triangulate the intersecting
    triangles and regonise the resulting triangles according to being
    inside/outside the surface.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "searchableSurfaceModifier/searchableSurfaceModifier.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::noParallel();
    #include "include/addDictOption.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    const word dictName("surfacePatchDict");
    #include "include/setSystemRunTimeDictionaryIO.H"
    const IOdictionary meshDict(dictIO);


    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            runTime.constant(),         // instance
            triSurfaceMesh::meshSubDir, // local
            runTime,                    // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshDict.subDict("geometry"),
        meshDict.lookupOrDefault("singleRegionName", true)
    );


    const dictionary& surfacesDict = meshDict.subDict("surfaces");

    forAllConstIter(dictionary, surfacesDict, surfacesIter)
    {
        if (surfacesIter().isDict())
        {
            const word& surfName = surfacesIter().keyword();
            const dictionary& surfDict = surfacesIter().dict();

            // Look up surface
            searchableSurface& surf = allGeometry[surfName];


            bool changed = false;

            // Check for modifier on global level
            if (surfDict.found("type"))
            {
                autoPtr<searchableSurfaceModifier> modifier
                (
                    searchableSurfaceModifier::New
                    (
                        surfDict.lookup("type"),
                        allGeometry,
                        surfDict
                    )
                );

                if (modifier().modify(identity(surf.regions().size()), surf))
                {
                    changed = true;
                }
            }

            // Check for modifier on a per-region level
            if (surfDict.found("regions"))
            {
                const dictionary& regionsDict = surfDict.subDict("regions");
                forAllConstIter(dictionary, regionsDict, regionsIter)
                {
                    const dictionary& regionDict = regionsIter().dict();
                    const keyType& regionName = regionsIter().keyword();

                    autoPtr<searchableSurfaceModifier> modifier
                    (
                        searchableSurfaceModifier::New
                        (
                            regionDict.lookup("type"),
                            allGeometry,
                            regionDict
                        )
                    );

                    labelList regionIDs =
                        findStrings(regionName, surf.regions());
                    if (modifier().modify(regionIDs, surf))
                    {
                        changed = true;
                    }
                }
            }


            if (changed)
            {
                const fileName name(surf.name());
                surf.rename(name.lessExt() + "_patched." + name.ext());

                Info<< "Writing repatched surface " << name << " to "
                    << surf.name() << nl << endl;

                surf.write();
            }
            else
            {
                Info<< "Surface " << surf.name() << " unchanged" << nl << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
