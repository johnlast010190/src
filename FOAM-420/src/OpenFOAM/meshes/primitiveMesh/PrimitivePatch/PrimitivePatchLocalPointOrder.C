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
    (c) 2022 Esi Ltd.

Description
     Orders the local points on the patch for most efficient search

\*---------------------------------------------------------------------------*/

#include "containers/LinkedLists/user/SLList.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcLocalPointOrder() const
{
    // Note: Cannot use bandCompressing as point-point addressing does
    // not exist and is not considered generally useful.

    if (debug)
    {
        PoutInFunction << "Calculating local point order" << endl;
    }

    if (localPointOrderPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "local point order already calculated"
            << abort(FatalError);
    }

    const List<FaceType>& lf = localFaces();

    const labelListList& ff = faceFaces();

    boolList visitedFace(lf.size(), false);

    localPointOrderPtr_ = new labelList(meshPoints().size(), -1);

    labelList& pointOrder = *localPointOrderPtr_;

    boolList visitedPoint(pointOrder.size(), false);

    label nPoints = 0;

    forAll(lf, facei)
    {
        if (!visitedFace[facei])
        {
            SLList<label> faceOrder(facei);

            do
            {
                const label curFace = faceOrder.first();

                faceOrder.removeHead();

                if (!visitedFace[curFace])
                {
                    visitedFace[curFace] = true;

                    const labelList& curPoints = lf[curFace];

                    // mark points
                    forAll(curPoints, pointi)
                    {
                        if (!visitedPoint[curPoints[pointi]])
                        {
                            visitedPoint[curPoints[pointi]] = true;

                            pointOrder[nPoints] = curPoints[pointi];

                            nPoints++;
                        }
                    }

                    // add face neighbours to the list
                    const labelList& nbrs = ff[curFace];

                    forAll(nbrs, nbrI)
                    {
                        if (!visitedFace[nbrs[nbrI]])
                        {
                            faceOrder.append(nbrs[nbrI]);
                        }
                    }
                }
            } while (faceOrder.size());
        }
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating local point order" << endl;
    }
}


// ************************************************************************* //
