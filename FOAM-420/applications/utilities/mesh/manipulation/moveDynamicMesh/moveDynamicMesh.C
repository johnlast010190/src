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
    (c) 2016 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Application
    moveDynamicMesh

Group
    grpMeshManipulationUtilities

Description
    Mesh motion and topological mesh changes utility.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "sampledSurface/writers/vtk/vtkSurfaceWriter.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dump patch + weights to vtk file
void writeWeights
(
    const polyMesh& mesh,
    const scalarField& wghtSum,
    const primitivePatch& patch,
    const fileName& directory,
    const fileName& prefix,
    const word& timeName
)
{
    vtkSurfaceWriter writer;

    // Collect geometry
    labelList pointToGlobal;
    labelList uniqueMeshPointLabels;
    autoPtr<globalIndex> globalPoints;
    autoPtr<globalIndex> globalFaces;
    faceList mergedFaces;
    pointField mergedPoints;
    Foam::PatchTools::gatherAndMerge
    (
        mesh,
        patch.localFaces(),
        patch.meshPoints(),
        patch.meshPointMap(),

        pointToGlobal,
        uniqueMeshPointLabels,
        globalPoints,
        globalFaces,

        mergedFaces,
        mergedPoints
    );
    // Collect field
    scalarField mergedWeights;
    globalFaces().gather
    (
        UPstream::worldComm,
        labelList(UPstream::procID(UPstream::worldComm)),
        wghtSum,
        mergedWeights
    );

    if (Pstream::master())
    {
        writer.write
        (
            directory,
            prefix + "_" + timeName,
            meshedSurfRef
            (
                mergedPoints,
                mergedFaces
            ),
            "weightsSum",
            mergedWeights,
            false
        );
    }
}


void writeWeights(const polyMesh& mesh)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const word tmName(mesh.time().timeName());

    forAll(pbm, patchi)
    {
        if (isA<cyclicAMIPolyPatch>(pbm[patchi]))
        {
            const cyclicAMIPolyPatch& cpp =
                refCast<const cyclicAMIPolyPatch>(pbm[patchi]);

            if (cpp.owner())
            {
                Info<< "Calculating AMI weights between owner patch: "
                    << cpp.name() << " and neighbour patch: "
                    << cpp.nbrPatch().name() << endl;

                const AMIPatchToPatchInterpolation& ami =
                    cpp.AMI();

                writeWeights
                (
                    mesh,
                    ami.tgtWeightsSum(),
                    cpp.nbrPatch(),
                    "postProcessing",
                    "tgt",
                    tmName
                );
                writeWeights
                (
                    mesh,
                    ami.srcWeightsSum(),
                    cpp,
                    "postProcessing",
                    "src",
                    tmName
                );
            }
        }
    }
}



int main(int argc, char *argv[])
{
    #include "include/addRegionOption.H"
    argList::addBoolOption
    (
        "checkAMI",
        "check AMI weights"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createDynamicFvMesh.H"

    const bool checkAMI  = args.optionFound("checkAMI");

    if (checkAMI)
    {
        Info<< "Writing VTK files with weights of AMI patches." << nl << endl;
    }

    pimpleControl pimple(mesh);

    bool moveMeshOuterCorrectors
    (
        pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
    );

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
            }
        }

        mesh.checkMesh(true);

        if (checkAMI)
        {
            writeWeights(mesh);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
