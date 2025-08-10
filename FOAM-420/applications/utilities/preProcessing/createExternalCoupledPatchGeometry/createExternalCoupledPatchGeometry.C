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
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2015 OpenFOAM Foundation

Application
    createExternalCoupledPatchGeometry.

Group
    grpPreProcessingUtilities

Description
    Application to generate the patch geometry (points and faces) for use
    with the externalCoupled functionObject.

Usage
    \verbatim
    createExternalCoupledPatchGeometry \<patchGroup\> [OPTION]
    \endverbatim

    \param -commsDir \<commsDir\> \n
    Specify an alternative communications directory (default is comms
    in the case directory)

    \param -region \<name\> \n
    Specify an alternative mesh region.

    \param -regions (\<name1\> \<name2\> .. \<namen\>) \n
    Specify alternative mesh regions. The region names will be sorted
    alphabetically and a single composite name will be created
        \<nameX\>_\<nameY\>.._\<nameZ\>

    On execution, the combined patch geometry (points and faces) are output
    to the communications directory.

Note:
    The addressing is patch-local, i.e. point indices for each patch point
    used for face addressing starts at index 0.

See also
    functionObjects::externalCoupled

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "externalCoupled/externalCoupled.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/addRegionOption.H"
    #include "addRegionsOption.H"
    argList::validArgs.append("patchGroup");
    argList::addOption
    (
        "commsDir",
        "dir",
        "specify alternate communications directory. default is 'comms'"
    );
    #include "include/setRootCase.H"
    #include "include/createTime.H"

    wordList regionNames(1, fvMesh::defaultRegion);
    if (!args.optionReadIfPresent("region", regionNames[0]))
    {
        args.optionReadIfPresent("regions", regionNames);
    }

    const wordRe patchGroup(args.argRead<wordRe>(1));

    fileName commsDir(runTime.path()/"comms");
    args.optionReadIfPresent("commsDir", commsDir);


    // Make sure region names are in canonical order
    stableSort(regionNames);


    PtrList<const fvMesh> meshes(regionNames.size());
    forAll(regionNames, i)
    {
        Info<< "Create mesh " << regionNames[i] << " for time = "
            << runTime.timeName() << nl << endl;

        meshes.set
        (
            i,
            new fvMesh
            (
                Foam::IOobject
                (
                    regionNames[i],
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }


    functionObjects::externalCoupled::writeGeometry
    (
        UPtrList<const fvMesh>(meshes),
        commsDir,
        patchGroup
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
