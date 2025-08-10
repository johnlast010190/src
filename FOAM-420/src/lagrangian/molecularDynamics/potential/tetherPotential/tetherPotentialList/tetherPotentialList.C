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

\*---------------------------------------------------------------------------*/

#include "tetherPotential/tetherPotentialList/tetherPotentialList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetherPotentialList::readTetherPotentialDict
(
    const List<word>& siteIdList,
    const dictionary& tetherPotentialDict,
    const List<word>& tetherSiteIdList
)
{

    Info<< nl << "Building tether potentials." << endl;

    idMap_ = List<label>(siteIdList.size(), -1);

    label tetherMapIndex = 0;

    forAll(tetherSiteIdList, t)
    {
        word tetherPotentialName = tetherSiteIdList[t];

        label tetherId = findIndex(siteIdList, tetherPotentialName);

        if (tetherId == -1)
        {
            FatalErrorInFunction
                << nl
                << "No matching entry found in siteIdList for tether name "
                << tetherPotentialName
                << abort(FatalError);
        }
        else if (!tetherPotentialDict.found(tetherPotentialName))
        {
            FatalErrorInFunction
                << nl << "tether potential specification subDict "
                << tetherPotentialName << " not found"
                << abort(FatalError);
        }
        else
        {
            this->set
            (
                tetherMapIndex,
                tetherPotential::New
                (
                    tetherPotentialName,
                    tetherPotentialDict.subDict(tetherPotentialName)
                )
            );
        }

        idMap_[tetherId] = tetherMapIndex;

        tetherMapIndex++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetherPotentialList::tetherPotentialList()
:
    PtrList<tetherPotential>(),
    idMap_()
{}


Foam::tetherPotentialList::tetherPotentialList
(
    const List<word>& siteIdList,
    const dictionary& tetherPotentialDict,
    const List<word>& tetherSiteIdList
)
:
    PtrList<tetherPotential>(),
    idMap_()
{
    buildPotentials(siteIdList, tetherPotentialDict, tetherSiteIdList);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetherPotentialList::~tetherPotentialList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetherPotentialList::buildPotentials
(
    const List<word>& siteIdList,
    const dictionary& tetherPotentialDict,
    const List<word>& tetherSiteIdList
)
{
    setSize(tetherSiteIdList.size());

    readTetherPotentialDict(siteIdList, tetherPotentialDict, tetherSiteIdList);
}


const Foam::tetherPotential& Foam::tetherPotentialList::tetherPotentialFunction
(
    const label a
) const
{
    return (*this)[tetherPotentialIndex(a)];
}


Foam::vector Foam::tetherPotentialList::force
(
    const label a,
    const vector rIT
) const
{
    return (*this)[tetherPotentialIndex(a)].force(rIT);
}


Foam::scalar Foam::tetherPotentialList::energy
(
    const label a,
    const vector rIT
) const
{
    return (*this)[tetherPotentialIndex(a)].energy(rIT);
}


// ************************************************************************* //
