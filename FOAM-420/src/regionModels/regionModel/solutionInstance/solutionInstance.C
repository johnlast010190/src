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
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solutionInstance/solutionInstance.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/UnbracketedTuple2/UnbracketedTuple2.H"
#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionInstance::solutionInstance
(
    const Time& runTime,
    const List<List<word>>& instanceRegions,
    const dictionary& dict,
    const word& name,
    const bool backwardCompatibility
)
:
    solutionInstanceRegistry(runTime, instanceRegions, dict, name),
    solutionRegistryPtrList_(regionNames_.size()),
    meshPtrList_(meshNames_.size()),
    backwardCompatibility_(backwardCompatibility)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionInstance::~solutionInstance() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvSolutionRegistry &Foam::solutionInstance::getSolutionRegistry
(
    const Time& runTime
)
{
    const word& meshName = meshNames_[0];
    const word& regionName = solutionRegions_[0][0];
    fvMesh* meshPtr = nullptr;
    if (runTime.foundObject<fvMesh>(meshName))
    {
        Info<< "Using Existing Mesh " << meshName << nl<<endl;
        meshPtr = runTime.lookupObjectRefPtr<fvMesh>(meshName);
    }
    else
    {
        IOobject dynamicMeshDictObj
        (
            "dynamicMeshDict",
            runTime.constant(),
            meshName == fvMesh::defaultRegion ? "" : meshName,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!dynamicMeshDictObj.typeHeaderOk<IOdictionary>(true))
        {
            Info<<"Constructing Static Mesh "<<meshName<<nl<<endl;
            meshPtr =
                new fvMesh
                (
                    IOobject
                    (
                        meshName,
                        runTime.timeName(),
                        runTime,
                        IOobject::MUST_READ
                    )
                );

            meshPtr->objectRegistry::store();
        }
        else
        {
            Info<<"Constructing Dynamic Mesh "<<meshName<<nl<<endl;
            meshPtr =
                 dynamicFvMesh::New
                 (
                    IOobject
                    (
                        meshName,
                        runTime.timeName(),
                        runTime,
                        IOobject::MUST_READ
                    )
                ).ptr();

            meshPtr->objectRegistry::store();
        }
    }
    meshPtrList_.set(0, meshPtr);

    if (!solutionRegistryPtrList_.set(0))
    {
        solutionRegistryPtrList_.set
        (
            0,
            new fvSolutionRegistry
            (
                runTime, *meshPtr, regionName, backwardCompatibility_
            )
        );
        Info<< "Constructing Registry: ";
        Info<< solutionRegistryPtrList_[0].regionName()<<nl<<endl;
    }
    active_ = true;
    return solutionRegistryPtrList_[0];
}


const Foam::PtrList<Foam::fvSolutionRegistry>&
Foam::solutionInstance::getSolutionRegistryList(const Time& runTime)
{
    if (backwardCompatibility_)
    {
        return getSolutionRegistryListOldFormat(runTime);
    }

    label registryCounter = 0;
    forAll(meshNames_, meshi)
    {
        forAll(solutionRegions_[meshi], regioni)
        {
            const word& meshName = meshNames_[meshi];
            const word& regionName = solutionRegions_[meshi][regioni];

            if (!solutionRegistryPtrList_.set(registryCounter))
            {
                solutionRegistryPtrList_.set
                (
                    registryCounter,
                    new fvSolutionRegistry
                    (
                        runTime,
                        *lookupOrCreateMesh(runTime, meshName),
                        regionName,
                        backwardCompatibility_
                    )
                );
                Info<< "Constructing Registry: "
                 << solutionRegistryPtrList_[registryCounter].regionName()<<nl<<endl;
                registryCounter++;
            }
        }
        active_ = true;
    }

    return solutionRegistryPtrList_;
}

const Foam::PtrList<Foam::fvSolutionRegistry>&
Foam::solutionInstance::getSolutionRegistryListOldFormat(const Time& runTime)
{
    forAll(meshNames_, meshi)
    {
        const word& meshName = meshNames_[meshi];
        const word& regionName = regionNames_[meshi];
        if (!solutionRegistryPtrList_.set(meshi))
        {
            solutionRegistryPtrList_.set
            (
                meshi,
                new fvSolutionRegistry
                (
                runTime,
                *lookupOrCreateMesh(runTime, meshName),
                regionName,
                backwardCompatibility_
                )
            );
            Info<< "Constructing Registry: ";
            Info<< solutionRegistryPtrList_[meshi].regionName()<<nl<<endl;
        }
        active_ = true;
    }
    return solutionRegistryPtrList_;
}


Foam::fvMesh* Foam::solutionInstance::lookupOrCreateMesh
(
    const Time& runTime,
    const word& meshName
)
{
    fvMesh* meshPtr = nullptr;

    if (runTime.foundObject<fvMesh>(meshName))
    {
        Info<< "Using existing mesh " << meshName << nl << endl;
        meshPtr = runTime.lookupObjectRefPtr<fvMesh>(meshName);
    }
    else
    {
        const bool meshWithChangers =
            runTime.controlDict().lookupOrDefault<bool>
            (
                "meshChangers",
                false
            );

        IOobject dynamicMeshDictObj
        (
            "dynamicMeshDict",
            runTime.constant(),
            meshName == fvMesh::defaultRegion ? "" : meshName,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!dynamicMeshDictObj.typeHeaderOk<IOdictionary>(true))
        {
            if (meshWithChangers)
            {
                // Create static mesh with dummy changers. Make the solver
                // compatible with static non-conformal couplings.
                Info<< "Create static mesh with changers for time = "
                    << runTime.timeName() << nl << endl;
            }
            else
            {
                Info<< "Constructing static mesh " << meshName << nl << endl;
            }

            meshPtr =
                new fvMesh
                (
                    IOobject
                    (
                        meshName,
                        runTime.timeName(),
                        runTime,
                        IOobject::MUST_READ
                    ),
                    meshWithChangers
                );

            meshPtr->objectRegistry::store();
            Info<< endl;
        }
        else
        {
            if (meshWithChangers)
            {
                // Create mesh with changers
                // (mesh movers and/or topo changers)
                Info<< "Create mesh with changers for time = "
                    << runTime.timeName() << nl << endl;

                meshPtr = new fvMesh
                (
                    IOobject
                    (
                        meshName,
                        runTime.timeName(),
                        runTime,
                        IOobject::MUST_READ
                    ),
                    true
                );

                meshPtr->objectRegistry::store();
                Info<< endl;
            }
            else
            {
                Info<< "Constructing dynamic mesh " << meshName << nl << endl;

                meshPtr =
                    dynamicFvMesh::New
                    (
                        IOobject
                        (
                            meshName,
                            runTime.timeName(),
                            runTime,
                            IOobject::MUST_READ
                        )
                    ).ptr();

                meshPtr->objectRegistry::store();
            }
        }

        // Test for the presence of non-conformal patches in a mesh
        // without changers
        if (!meshWithChangers)
        {
            forAll(meshPtr->boundaryMesh(), patchi)
            {
                const polyPatch& pp = meshPtr->boundaryMesh()[patchi];

                if (isA<Foam::nonConformalPolyPatch>(pp))
                {
                    FatalErrorInFunction
                        << "Non-conformal patches are not compatible with old "
                        << "dynamic mesh structure, without changers." << nl
                        << "Please set the 'meshChangers' flag to 'true' "
                        << "in the controlDict dictionary or remove any "
                        << "non-conformal patch."
                        << exit(FatalError);
                }
            }
        }
    }

    auto meshi = std::find(meshNames_.begin(), meshNames_.end(), meshName);
    if (meshi != meshNames_.end())
    {
        meshPtrList_.set(meshi - meshNames_.begin(), meshPtr);
    }

    return meshPtr;
}


void Foam::solutionInstance::releaseUnusedMeshes
(
    const Time& runTime,
    const word &meshName
)
{
    forAll(meshNames_, meshI)
    {
        if (meshNames_[meshI] == meshName)
        {
            runTime.objectRegistry::checkOut(meshNames_[meshI]);
            // Pointer was deleted as owned by the registry
            meshPtrList_.set(meshI, nullptr);
        }
    }
}

void Foam::solutionInstance::deactivate()
{
    solutionRegistryPtrList_.clear();
    active_ = false;
}


// ************************************************************************* //
