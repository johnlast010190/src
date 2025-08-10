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

#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/UnbracketedTuple2/UnbracketedTuple2.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::label Foam::solutionInstanceRegistry::getTotalNumberOfRegions()
{
    label totalNumberOfRegions = 0;
    forAll(solutionRegions_, ri)
    {
        totalNumberOfRegions += solutionRegions_[ri].size();
    }
    return totalNumberOfRegions;
}


void Foam::solutionInstanceRegistry::generateMeshList()
{
    forAll(regionNames_, rN)
    {
        word meshName = solutionMeshDict_.lookup(regionNames_[rN]);
        if (!meshNames_.found(meshName))
        {
            meshNames_.append(meshName);
        }
    }
}


Foam::label Foam::solutionInstanceRegistry::whichRegionIndex
(
    const word& regionName
) const
{
    label regionIndex = -1;
    forAll(regionNames_, rN)
    {
        const word& regionNameIn = regionNames_[rN];
        if (regionNameIn==regionName)
        {
            regionIndex = rN;
        }
    }
    return regionIndex;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionInstanceRegistry::solutionInstanceRegistry
(
    const objectRegistry& obr,
    const List<List<word>>& instanceRegions,
    const dictionary& dict,
    const word& name
)
:
    objectRegistry
    (
        IOobject
        (
            name,
            obr.path(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    solutionMeshDict_(dict),
    meshNames_(0),
    solutionRegions_(instanceRegions),
    regionNames_( getTotalNumberOfRegions()),
    active_(false),
    selected_(true)
{
    label regionCounter = 0;
    forAll(solutionRegions_, ri)
    {
        forAll(solutionRegions_[ri], rii)
        {
            regionNames_[regionCounter++] = solutionRegions_[ri][rii];
        }
    }
    generateMeshList();
    Info<<"Building instance: "<<IOobject::name()<<nl<<endl;
}

Foam::solutionInstanceRegistry::~solutionInstanceRegistry() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::solutionInstanceRegistry::whichMesh(const word& regionName)
const
{
    label regionIndex = whichRegionIndex(regionName);

    if (regionIndex == -1)
    {
        FatalErrorInFunction
            << "Requested mesh name for region " << regionName
            << ", does not exist!"
            << exit(FatalError);
    }

    word meshName = solutionMeshDict_.lookup(regionNames_[regionIndex]);

    return meshName;
}
