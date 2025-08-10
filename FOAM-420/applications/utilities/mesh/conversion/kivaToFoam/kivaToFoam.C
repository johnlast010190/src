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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Application
    kivaToFoam

Group
    grpMeshConversionUtilities

Description
    Converts a KIVA3v grid to OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/meshShapes/cellShape/cellShape.H"
#include "meshes/meshShapes/cellModeller/cellModeller.H"
#include "meshes/preservePatchTypes/preservePatchTypes.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/symmetry/symmetryPolyPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.H"
#include "mergedCyclic/mergedCyclicPolyPatch.H"
#include "mergedCyclic/polyMeshUnMergeCyclics.H"
#include "global/unitConversion/unitConversion.H"

using namespace Foam;

//- Supported KIVA versions
enum kivaVersions
{
    kiva3,
    kiva3v
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption
    (
        "file",
        "name",
        "specify alternative input file name - default is otape17"
    );
    argList::addOption
    (
        "version",
        "version",
        "specify kiva version [kiva3|kiva3v] - default is '3v'"
    );
    argList::addOption
    (
        "zHeadMin",
        "scalar",
        "minimum z-height for transferring liner faces to cylinder-head"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    const fileName kivaFileName =
        args.optionLookupOrDefault<fileName>("file", "otape17");

    kivaVersions kivaVersion = kiva3v;
    if (args.optionFound("version"))
    {
        const word versionName = args["version"];

        if (versionName == "kiva3")
        {
            kivaVersion = kiva3;
        }
        else if (versionName == "kiva3v")
        {
            kivaVersion = kiva3v;
        }
        else
        {
            FatalErrorInFunction
                << "KIVA file version " << versionName << " not supported"
                << exit(FatalError);

            args.printUsage();
            FatalError.exit(1);
        }
    }

    scalar zHeadMin = -GREAT;
    args.optionReadIfPresent("zHeadMin", zHeadMin);

    #include "readKivaGrid.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
