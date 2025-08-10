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
    (c) 2009-2009 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedProfile/pointProfileDataList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointProfileDataList::pointProfileDataList()
:
    PtrList<pointProfileData>()
{}


Foam::pointProfileDataList::pointProfileDataList(ISstream& is)
:
    PtrList<pointProfileData>()
{
    read(is);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointProfileDataList::~pointProfileDataList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarField* Foam::pointProfileDataList::mapData
(
    const pointField& points,
    word fieldName,
    word regionName
)
{

    PtrList<pointProfileData>& ppdl(*this);

    List<string> availableRegions(ppdl.size());

    forAll(ppdl, rI)
    {
        if (regionName == ppdl[rI].regionName())
        {
            return ppdl[rI].mapData(points, fieldName);
        }
        else
        {
            availableRegions[rI] = ppdl[rI].regionName();
        }
    }

    WarningInFunction
        << "Could not find boundary region matching descriptor "
        <<  regionName
        << nl << "Available regions are: " << availableRegions
        << endl;


    return nullptr;
}


void Foam::pointProfileDataList::read(Foam::ISstream& is)
{

    PtrList<pointProfileData>& ppdl(*this);
    label nRegions = 0;
    while (!is.eof())
    {
        //read leading bracket
        token currToken(is);
        if (is.eof()) break;

        pointProfileData::checkIO(currToken, token::BEGIN_LIST);

        nRegions++;

        ppdl.setSize(nRegions);

        ppdl.set(nRegions-1, new pointProfileData(is));

    }

}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// ************************************************************************* //
