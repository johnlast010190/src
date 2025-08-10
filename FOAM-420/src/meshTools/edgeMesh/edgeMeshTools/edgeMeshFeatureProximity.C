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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "edgeMesh/edgeMeshTools/edgeMeshTools.H"

#include "edgeMesh/extendedEdgeMesh/extendedEdgeMesh.H"
#include "triSurface/triSurface.H"
#include "triSurface/fields/triSurfaceFields.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "MeshedSurface/MeshedSurface.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

static scalar calcProximityOfFeaturePoints
(
    const List<pointIndexHit>& hitList,
    const scalar defaultCellSize
)
{
    scalar minDist = defaultCellSize;

    for
    (
        label hI1 = 0;
        hI1 < hitList.size() - 1;
        ++hI1
    )
    {
        const pointIndexHit& pHit1 = hitList[hI1];

        if (pHit1.hit())
        {
            for
            (
                label hI2 = hI1 + 1;
                hI2 < hitList.size();
                ++hI2
            )
            {
                const pointIndexHit& pHit2 = hitList[hI2];

                if (pHit2.hit())
                {
                    scalar curDist = mag(pHit1.hitPoint() - pHit2.hitPoint());

                    minDist = min(curDist, minDist);
                }
            }
        }
    }

    return minDist;
}


scalar calcProximityOfFeatureEdges
(
    const edgeMesh& emesh,
    const List<pointIndexHit>& hitList,
    const scalar defaultCellSize
)
{
    scalar minDist = defaultCellSize;

    for
    (
        label hI1 = 0;
        hI1 < hitList.size() - 1;
        ++hI1
    )
    {
        const pointIndexHit& pHit1 = hitList[hI1];

        if (pHit1.hit())
        {
            const edge& e1 = emesh.edges()[pHit1.index()];

            for
            (
                label hI2 = hI1 + 1;
                hI2 < hitList.size();
                ++hI2
            )
            {
                const pointIndexHit& pHit2 = hitList[hI2];

                if (pHit2.hit())
                {
                    const edge& e2 = emesh.edges()[pHit2.index()];

                    // Don't refine if the edges are connected to each other
                    if (!e1.connects(e2))
                    {
                        scalar curDist =
                            mag(pHit1.hitPoint() - pHit2.hitPoint());

                        minDist = min(curDist, minDist);
                    }
                }
            }
        }
    }

    return minDist;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::edgeMeshTools::featureProximity
(
    const extendedEdgeMesh& emesh,
    const triSurface& surf,
    const scalar searchDistance
)
{
    tmp<scalarField> tfld(new scalarField(surf.size(), searchDistance));
    scalarField& featureProximity = tfld.ref();

    Info<< "Extracting proximity of close feature points and "
        << "edges to the surface" << endl;

    forAll(surf, fI)
    {
        const triPointRef& tri = surf[fI].tri(surf.points());
        const point& triCentre = tri.circumCentre();

        const scalar radiusSqr = min
        (
            sqr(4*tri.circumRadius()),
            sqr(searchDistance)
        );

        List<pointIndexHit> hitList;

        emesh.allNearestFeatureEdges(triCentre, radiusSqr, hitList);

        featureProximity[fI] =
            calcProximityOfFeatureEdges
            (
                emesh,
                hitList,
                featureProximity[fI]
            );

        emesh.allNearestFeaturePoints(triCentre, radiusSqr, hitList);

        featureProximity[fI] =
            calcProximityOfFeaturePoints
            (
                hitList,
                featureProximity[fI]
            );
    }

    return tfld;
}


Foam::tmp<Foam::scalarField> Foam::edgeMeshTools::writeFeatureProximity
(
    const Time& runTime,
    const word& basename,
    const extendedEdgeMesh& emesh,
    const triSurface& surf,
    const scalar searchDistance
)
{
    Info<< nl << "Extracting curvature of surface at the points."
        << endl;


    tmp<scalarField> tfld =
        edgeMeshTools::featureProximity(emesh, surf, searchDistance);
    scalarField& featureProximity = tfld.ref();

    triSurfaceScalarField outputField
    (
        IOobject
        (
            basename + ".featureProximity",
            runTime.constant(),
            "triSurface",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        scalarField()
    );

    outputField.swap(featureProximity);
    outputField.write();
    outputField.swap(featureProximity);

    return tfld;
}

// ************************************************************************* //
