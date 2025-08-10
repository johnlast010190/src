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

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcBdryPoints() const
{
    if (debug)
    {
        InfoInFunction << "Calculating boundary points" << endl;
    }

    if (boundaryPointsPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "edge types already calculated"
            << abort(FatalError);
    }

    const edgeList& e = edges();

    labelHashSet bp(2*e.size());

    for (label edgeI = nInternalEdges_; edgeI < e.size(); edgeI++)
    {
        const edge& curEdge = e[edgeI];

        bp.insert(curEdge.start());
        bp.insert(curEdge.end());
    }

    boundaryPointsPtr_ = new labelList(bp.toc());
    sort(*boundaryPointsPtr_);

    if (debug)
    {
        Info<< "    Finished." << endl;
    }
}


// ************************************************************************* //
