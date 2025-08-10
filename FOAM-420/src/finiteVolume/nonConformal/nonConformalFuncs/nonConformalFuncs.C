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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/nonConformalFuncs/nonConformalFuncs.H"
#include "fvMesh/fvMeshStitchers/stationary/fvMeshStitchersStationary.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::nonConformalFuncs::generateCouples
(
    fvMesh& mesh,
    const List<Pair<word>>& patchNames,
    const wordList& cyclicNames,
    const List<cyclicTransform>& transforms
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> newPatches;

    // Find the first processor patch and face
    label firstProcPatchi = patches.size(), firstProcFacei = mesh.nFaces();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<processorPolyPatch>(pp) && firstProcPatchi == patches.size())
        {
            firstProcPatchi = patchi;
            firstProcFacei = pp.start();
        }

        if (!isA<processorPolyPatch>(pp) && firstProcPatchi != patches.size())
        {
            FatalErrorInFunction
                << "Processor patches do not follow boundary patches"
                << exit(FatalError);
        }
    }

    // Clone the non-processor patches
    for (label patchi = 0; patchi < firstProcPatchi; ++patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            newPatches.append
            (
                dpp.clone(patches, patchi, dpp.size(), dpp.start()).ptr()
            );
        }
    }

    // Convenience function to generate patch names for the owner or neighbour
    auto nccPatchNames = [&](const label i)
    {
        return
            Pair<word>
            (
                cyclicNames[i] + "_on_" + patchNames[i][0],
                cyclicNames[i] + "_on_" + patchNames[i][1]
            );
    };

    // Add the cyclic patches
    forAll(patchNames, i)
    {
        Info<< indent << "Adding "
            << nonConformalCyclicPolyPatch::typeName
            << " interfaces between patches: " << incrIndent << nl
            << indent << patchNames[i] << decrIndent << nl
            << indent << "Named:" << incrIndent << nl
            << indent << Pair<word>(nccPatchNames(i)) << decrIndent << nl
            << indent << "With transform: " << incrIndent << nl;
        transforms[i].write(Info);
        Info<< decrIndent << nl;

        newPatches.append
        (
            new nonConformalCyclicPolyPatch
            (
                nccPatchNames(i)[0],
                0,
                firstProcFacei,
                newPatches.size(),
                patches,
                nonConformalCyclicPolyPatch::typeName,
                nccPatchNames(i)[1],
                patchNames[i][0],
                transforms[i]
            )
        );
        newPatches.append
        (
            new nonConformalCyclicPolyPatch
            (
                nccPatchNames(i)[1],
                0,
                firstProcFacei,
                newPatches.size(),
                patches,
                nonConformalCyclicPolyPatch::typeName,
                nccPatchNames(i)[0],
                patchNames[i][1],
                inv(transforms[i])
            )
        );
    }

    // Add the error patches. Note there is only one for each source patch,
    // regardless of how many interfaces are attached to that patch.
    auto appendErrorPatches = [&](const bool owner)
    {
        wordHashSet patchANames;
        forAll(patchNames, i)
        {
            patchANames.insert(patchNames[i][!owner]);
        }
        forAllConstIter(wordHashSet, patchANames, iter)
        {
            newPatches.append
            (
                new nonConformalErrorPolyPatch
                (
                    nonConformalErrorPolyPatch::typeName + "_on_" + iter.key(),
                    0,
                    firstProcFacei,
                    newPatches.size(),
                    patches,
                    nonConformalErrorPolyPatch::typeName,
                    iter.key()
                )
            );
        }
    };
    appendErrorPatches(true);
    appendErrorPatches(false);

    // Clone the processor patches
    for (label patchi = firstProcPatchi; patchi < patches.size(); ++patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            newPatches.append
            (
                dpp.clone
                (
                    patches,
                    newPatches.size(),
                    dpp.size(),
                    dpp.start()
                ).ptr()
            );
        }
    }

    // Add the processor cyclic patches
    if (Pstream::parRun())
    {
        forAll(patchNames, i)
        {
            const polyPatch& patch1 = patches[patchNames[i][0]];
            const polyPatch& patch2 = patches[patchNames[i][1]];

            boolList procHasPatch1(Pstream::nProcs(), false);
            procHasPatch1[Pstream::myProcNo()] = !patch1.empty();
            Pstream::allGatherList(procHasPatch1);

            boolList procHasPatch2(Pstream::nProcs(), false);
            procHasPatch2[Pstream::myProcNo()] = !patch2.empty();
            Pstream::allGatherList(procHasPatch2);

            // Multiple cyclic interfaces must be ordered in a specific way for
            // processor communication to function correctly.
            //
            // A communication that is sent from the cyclic owner is received
            // on the cyclic neighbour and vice versa. Therefore, in a coupled
            // pair of processors if one sends the owner first the other must
            // receive the neighbour first.
            //
            // We ensure the above by ordering the patches so that for the
            // lower indexed processor the owner interface comes first, and for
            // the higher indexed processor the neighbour comes first.

            auto appendProcPatches = [&](const bool owner, const bool first)
            {
                const boolList& procHasPatchA =
                    owner ? procHasPatch1 : procHasPatch2;
                const boolList& procHasPatchB =
                    owner ? procHasPatch2 : procHasPatch1;

                if (procHasPatchA[Pstream::myProcNo()])
                {
                    forAll(procHasPatchB, proci)
                    {
                        if
                        (
                            (
                                (first && proci > Pstream::myProcNo())
                             || (!first && proci < Pstream::myProcNo())
                            )
                         && procHasPatchB[proci]
                        )
                        {
                            newPatches.append
                            (
                                new nonConformalProcessorCyclicPolyPatch
                                (
                                    0,
                                    mesh.nFaces(),
                                    newPatches.size(),
                                    patches,
                                    Pstream::myProcNo(),
                                    proci,
                                    nccPatchNames(i)[!owner],
                                    patchNames[i][!owner]
                                )
                            );
                        }
                    }
                }
            };

            appendProcPatches(true, true);
            appendProcPatches(false, true);
            appendProcPatches(false, false);
            appendProcPatches(true, false);
        }
    }

    // Re-patch the mesh
    mesh.removeFvBoundary();
    mesh.addFvPatches(newPatches);

    // Calculate the geometry for the patches (e.g. the transform), if needed
    mesh.recalculatePatches();
}


void Foam::nonConformalFuncs::createNonConformalCouples(fvMesh& mesh)
{
    Info<< "Generate non-conformal coupling" << endl;

    // Patch names between which to create couples, the associated cyclic name
    // prefix and transformation (if any)
    List<Pair<word>> patchNames;
    wordList cyclicNames;
    List<cyclicTransform> transforms;

    // Get patch names and transformation types from dictionary
    const word dictName("nonConformalCouplesDict");
    Info<< "Reading " << dictName << nl << endl;
    IOdictionary dict
    (
        IOobject
        (
            dictName,
            mesh.time().system(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict()) continue;

        const Pair<word> origPatches =
            iter().dict().lookup<Pair<word>>("patches");

        const word couplingType =
            iter().dict().lookupOrDefault<word>("couplingType", "standard");

        if (couplingType == "standard")
        {
            patchNames.append(origPatches);
            cyclicNames.append(iter().dict().dictName());
            transforms.append(cyclicTransform(iter().dict(), true));
        }
        else if (couplingType == "periodic")
        {
            const word transformType =
                iter().dict().lookup<word>("transformType");

            if (transformType != "rotational")
            {
                FatalErrorInFunction
                    << "Invalid transform type for periodic non-conformal"
                    << " coupling: only 'rotational' transform is supported."
                    << exit(FatalError);
            }

            const word baseName = iter().dict().dictName();

            // First, create a couple with no transform
            patchNames.append(origPatches);
            cyclicNames.append(baseName + "0");
            transforms.append(cyclicTransform(true));

            // Then, append rotational couples up to a full turn.
            dictionary coupleDict;
            coupleDict.add("patches", origPatches);
            coupleDict.add("transformType", transformType);
            coupleDict.add
            (
                "rotationAxis",
                iter().dict().lookup<vector>("rotationAxis")
            );
            coupleDict.add
            (
                "rotationCentre",
                iter().dict().lookup<point>("rotationCentre")
            );

            const scalar secAngle = iter().dict().lookup<scalar>("sectorAngle");

            scalar angle = secAngle;
            while (Foam::mag(angle) < 360)
            {
                coupleDict.add("rotationAngle", angle, true);

                patchNames.append(origPatches);
                cyclicNames.append
                (
                    word(baseName + name(std::round(Foam::mag(angle))))
                );
                transforms.append(cyclicTransform(coupleDict, true));

                angle += secAngle;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Invalid coupling type for NCC. Valid types are 'standard'"
                << " or 'periodic'." << exit(FatalError);
        }
    }

    // Save mesh instance for writing
    const word meshInstance = mesh.pointsInstance();

    // Make sure the mesh is not connected before couples are added
    fvMeshStitchers::stationary stitcher(mesh);
    stitcher.disconnect(false, false);

    // Generate all non-conformal couples
    generateCouples(mesh, patchNames, cyclicNames, transforms);

    // Connect the mesh so that the new stitching topology gets written out
    stitcher.connect(false, true, false);

    // Write resulting mesh
    mesh.setInstance(meshInstance);
    Info<< "Writing non-conformal mesh to "
        << mesh.pointsInstance() << nl << endl;

    mesh.write();
}


// ************************************************************************* //
