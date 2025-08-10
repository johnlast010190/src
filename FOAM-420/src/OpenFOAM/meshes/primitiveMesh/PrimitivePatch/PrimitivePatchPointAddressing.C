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
     Point addressing on the patch: pointEdges and pointFaces.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
#include "containers/LinkedLists/user/SLList.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcPointEdges() const
{
    if (debug)
    {
        InfoInFunction << "Calculating pointEdges" << endl;
    }

    if (pointEdgesPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "pointEdges already calculated"
            << abort(FatalError);
    }

    pointEdgesPtr_ = new labelListList(meshPoints().size());

    labelListList& pe = *pointEdgesPtr_;

    invertManyToMany(pe.size(), edges(), pe);

    if (debug)
    {
        Info<< "    Finished." << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcPointFaces() const
{
    if (debug)
    {
        InfoInFunction << "Calculating pointFaces" << endl;
    }

    if (pointFacesPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "pointFaces already calculated"
            << abort(FatalError);
    }

    const List<FaceType>& f = localFaces();

    // set up storage for pointFaces
    List<SLList<label>> pointFcs(meshPoints().size());

    forAll(f, facei)
    {
        const FaceType& curPoints = f[facei];

        forAll(curPoints, pointi)
        {
            pointFcs[curPoints[pointi]].append(facei);
        }
    }

    // sort out the list
    pointFacesPtr_ = new labelListList(pointFcs.size());

    labelListList& pf = *pointFacesPtr_;

    forAll(pointFcs, pointi)
    {
        pf[pointi].setSize(pointFcs[pointi].size());

        label i = 0;
        forAllIter(SLList<label>, pointFcs[pointi], curFacesIter)
        {
            pf[pointi][i++] = curFacesIter();
        }
    }

    if (debug)
    {
        Info<< "    Finished." << endl;
    }
}


// ************************************************************************* //
