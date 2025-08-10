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

\*---------------------------------------------------------------------------*/

#include "slidingInterface/enrichedPatch/enrichedPatch.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcPointPoints() const
{
    // Calculate point-point addressing
    if (pointPointsPtr_)
    {
        FatalErrorInFunction
            << "Point-point addressing already calculated."
            << abort(FatalError);
    }

    // Algorithm:
    // Go through all faces and add the previous and next point as the
    // neighbour for each point. While inserting points, reject the
    // duplicates (as every internal edge will be visited twice).
    List<DynamicList<label, primitiveMesh::edgesPerPoint_>>
        pp(meshPoints().size());

    const faceList& lf = localFaces();

    bool found = false;

    forAll(lf, facei)
    {
        const face& curFace = lf[facei];

        forAll(curFace, pointi)
        {
            DynamicList<label, primitiveMesh::edgesPerPoint_>&
                curPp = pp[curFace[pointi]];

            // Do next label
            label next = curFace.nextLabel(pointi);

            found = false;

            forAll(curPp, i)
            {
                if (curPp[i] == next)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                curPp.append(next);
            }

            // Do previous label
            label prev = curFace.prevLabel(pointi);
            found = false;

            forAll(curPp, i)
            {
                if (curPp[i] == prev)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                curPp.append(prev);
            }
        }
    }

    // Re-pack the list
    pointPointsPtr_ = new labelListList(pp.size());
    labelListList& ppAddr = *pointPointsPtr_;

    forAll(pp, pointi)
    {
        ppAddr[pointi].transfer(pp[pointi]);
    }
}


// ************************************************************************* //
