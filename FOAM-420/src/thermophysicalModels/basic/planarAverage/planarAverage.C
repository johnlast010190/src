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
    (c) 2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "planarAverage.H"
#include "fields/Fields/primitiveFields.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "fvMesh/fvMesh.H"
#include "faMesh/faMesh.H"
#include "faEdgeMesh/faEdgeMesh.H"
#include "fields/areaFields/areaFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/slicedVolFields.H"

namespace Foam
{
namespace planarAverage
{
int debug(Foam::debug::debugSwitch("planarAverage", 0));
}
registerDebugSwitchWithName(planarAverage, planarAverage, "planarAverage");
}

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

void Foam::planarAverage::addToAverage
(
    const vector& ownerPoint,
    const vector& nbrPoint,
    const scalar& hMin,
    const scalar& hMax,
    const label nPlanes,
    const vector& n,
    const scalar& ownerf,
    const scalar& nbrf,
    scalarField& avg,
    labelField& nFound
)
{
    // Start/end segment node of the edge connecting cell centres
    scalar nodea = ((ownerPoint & n) - hMin)/(hMax-hMin)*(nPlanes-1);
    scalar nodeb = ((nbrPoint & n) - hMin)/(hMax-hMin)*(nPlanes-1);
    if (nodea > nodeb)
    {
        scalar temp = nodea;
        nodea = nodeb;
        nodeb = temp;
    }
    // Make sure nodes are selected that are only just overlapping the mesh
    nodea -= ROOTSMALL;
    nodeb += ROOTSMALL;
    label nodej = ceil(nodea);
    label nodek = floor(nodeb);
    if (debug)
    {
        if (nodej < 0 || nodej >= avg.size())
        {
            Pout<< "Unexpected out-of-range: " << nodej << " " << nodea << " "
                << ownerPoint << " " << n << " " << hMin << " " << hMax << " "
                << nPlanes << endl;
        }
        if (nodek < 0 || nodek >= avg.size())
        {
            Pout<< "Unexpected out-of-range: " << nodek << " " << nodeb << " "
                << nbrPoint << " " << n << " " << hMin << " " << hMax << " "
                << nPlanes << endl;
        }
    }
    nodej = min(max(nodej, label(0)), avg.size()-1);
    nodek = min(max(nodek, label(0)), avg.size()-1);

    for (label nodel = nodej; nodel <= nodek; nodel++)
    {
        scalar height = nodel*(hMax-hMin)/(nPlanes-1)+hMin;
        scalar dh = (nbrPoint & n) - (ownerPoint & n);
        scalar r = (mag(dh) < VSMALL ? 0.5 : (height - (ownerPoint & n))/dh);
        r = max(min(r, scalar(1.01)), scalar(-0.01));
        avg[nodel] += r*(nbrf-ownerf)+ownerf;
        ++nFound[nodel];
    }
}


template <class GeoScalarField, class GeoVectorField>
Foam::tmp<Foam::scalarField> Foam::planarAverage::createHeightProfile
(
    const GeoScalarField& f,
    const label nPlanes,
    const scalar& hMin,
    const scalar& hMax,
    const GeoVectorField& points,
    const labelList& ownerList,
    const labelList& neighbourList,
    const vector& n
)
{
    tmp<scalarField> tavg(new scalarField(nPlanes, 0.0));
    scalarField& avg = tavg.ref();
    labelField nFound(nPlanes, label(0));
    // Loop over all edges (i.e. connections between points) and assign
    // interpolated rho to nodes on the height profile that fall within the edge
    forAll(neighbourList, edgei)
    {
        label own = ownerList[edgei];
        label nei = neighbourList[edgei];

        addToAverage
        (
            points[own],
            points[nei],
            hMin,
            hMax,
            nPlanes,
            n,
            f[own],
            f[nei],
            avg,
            nFound
        );
    }

    const typename GeoScalarField::Boundary& bf(f.boundaryField());
    forAll(bf, patchi)
    {
        tmp<vectorField> ownerPoints =
            points.boundaryField()[patchi].patchInternalField();
        tmp<scalarField> ownerf =
            bf[patchi].patchInternalField();
        if (bf[patchi].patch().coupled() && bf[patchi].coupled())
        {
            tmp<vectorField> nbrPoints =
                points.boundaryField()[patchi].patchNeighbourField();
            tmp<scalarField> nbrf =
                bf[patchi].patchNeighbourField();
            forAll(bf[patchi], facei)
            {
                addToAverage
                (
                    ownerPoints()[facei],
                    nbrPoints()[facei],
                    hMin,
                    hMax,
                    nPlanes,
                    n,
                    ownerf()[facei],
                    nbrf()[facei],
                    avg,
                    nFound
                );
            }
        }
        else
        {
            const vectorField& bPoints = points.boundaryField()[patchi];
            const scalarField& pf = bf[patchi];
            forAll(pf, facei)
            {
                addToAverage
                (
                    ownerPoints()[facei],
                    bPoints[facei],
                    hMin,
                    hMax,
                    nPlanes,
                    n,
                    ownerf()[facei],
                    pf[facei],
                    avg,
                    nFound
                );
            }
        }
    }

    listReduce(avg, sumOp<scalar>());
    listReduce(nFound, sumOp<label>());

    forAll(avg, nodei)
    {
        // Nodes may not have any contributions if bbox too big
        if (nFound[nodei])
        {
            avg[nodei] /= nFound[nodei];
        }
    }

    // Interpolate/extrapolate to fill in any possible 'holes'
    label lastValid = -1;
    forAll(avg, nodei)
    {
        if (nFound[nodei])
        {
            for (label nodej = lastValid+1; nodej < nodei; nodej++)
            {
                if (lastValid == -1)
                {
                    avg[nodej] = avg[nodei];
                }
                else
                {
                    avg[nodej] =
                        avg[lastValid]
                      + (avg[nodei]-avg[lastValid])
                       *(nodej-lastValid)/(nodei-lastValid);
                }
            }
            lastValid = nodei;
        }
        if (debug)
        {
            if (!nFound[nodei])
            {
                Pout<< "Not found at " << nodei << endl;
            }
        }
    }

    if (avg.size() && !nFound[avg.size()-1])
    {
        if (lastValid != -1)
        {
            for (label nodej = lastValid+1; nodej < avg.size(); nodej++)
            {
                avg[nodej] = avg[lastValid];
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unable to create height profile."
                << exit(FatalError);
        }
    }

    return tavg;
}


void Foam::planarAverage::interpolateHeightProfileToPoints
(
    const scalarField& profileNodes,
    const scalar& hMin,
    const scalar& hMax,
    const vectorField& points,
    const vector& n,
    scalarField& f
)
{
    const label nPlanes = profileNodes.size();
    forAll(f, pi)
    {
        scalar node = ((points[pi] & n) - hMin)/(hMax-hMin)*(nPlanes-1);
        label node0 = floor(node);
        node0 = max(min(node0, nPlanes-2), 0);
        scalar r = node-node0;
        f[pi] =
            r*(profileNodes[node0+1]-profileNodes[node0])+profileNodes[node0];
    }
}


Foam::tmp<Foam::scalarField>
Foam::planarAverage::planarAverage
(
    const faMesh& aMesh,
    const scalarField& f,
    const vector& n,
    label nPlanes
)
{
    if (debug)
    {
        Pout<< "Enter planarAverage(faMesh& ...)" << endl;
    }

    scalar hMin = gMin(aMesh.points() & n);
    scalar hMax = gMax(aMesh.points() & n);

    if (nPlanes < 0)
    {
        scalar minLe = gMin(aMesh.magLe().primitiveField());
        nPlanes = 2*round(mag(hMax-hMin)/minLe)+1;
    }

    areaScalarField af
    (
        "af",
        aMesh.thisDb(),
        aMesh,
        dimless,
        "zeroGradient"
    );
    af.primitiveFieldRef() = f;
    af.correctBoundaryConditions();

    tmp<scalarField> profile =
        createHeightProfile
        (
            af,
            nPlanes,
            hMin,
            hMax,
            aMesh.areaCentres(),
            aMesh.edgeOwner(),
            aMesh.edgeNeighbour(),
            n
        );

    tmp<scalarField> fAvg(new scalarField(f));

    // Interpolate from the height segments to cell centres
    interpolateHeightProfileToPoints
    (
        profile, hMin, hMax, aMesh.areaCentres(), n, fAvg.ref()
    );

    if (debug)
    {
        Pout<< "Done planarAverage(faMesh& ...)" << endl;
    }

    return fAvg;
}


Foam::tmp<Foam::volScalarField>
Foam::planarAverage::planarAverage
(
    const volScalarField& f,
    const vector& n,
    label nPlanes
)
{
    if (debug)
    {
        Pout<< "Enter planarAverage(volScalarField& ...)" << endl;
    }

    const fvMesh& mesh = f.mesh();

    boundBox globalBounds = mesh.bounds();
    globalBounds.reduce();
    scalar hMin = globalBounds.min() & n;
    scalar hMax = globalBounds.max() & n;

    if (nPlanes < 0)
    {
        // Default to 1000
        nPlanes = 1000;
        // NOTE: We could default to the smallest edge length as below,
        // but this is likely to be overkill
        //nPlanes =
        //    mag(hMax-hMin)/gMin(mag(mesh.delta())().primitiveField())+1;
    }

    // Unfortunately mesh.C() does not treat coupled non-processor patches,
    // so we create our own which does
    /*
    slicedVolVectorField C
    (
        IOobject
        (
            "C",
            mesh.time().timeName(),
            f.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.C(),
        true                //preserveCouples
    );
    */
    // The above might have been responsible for issues, so we play it safe
    // and create a copy in a volScalarField instead
    volVectorField C
    (
        "C",
        f.db(),
        mesh,
        mesh.C().dimensions()
    );
    C.primitiveFieldRef() = mesh.C().primitiveField();
    volVectorField::Boundary& bf = C.boundaryFieldRef();
    forAll(bf, patchi)
    {
        bf[patchi].forceAssign(mesh.C().boundaryField()[patchi]);
    }
    C.correctBoundaryConditions();

    if (debug)
    {
        Pout<< "Call createHeightProfile" << endl;
    }

    tmp<scalarField> profile =
        createHeightProfile
        (
            f,
            nPlanes,
            hMin,
            hMax,
            C,
            mesh.owner(),
            mesh.neighbour(),
            n
        );

    if (debug)
    {
        Pout<< "Done createHeightProfile" << endl;
    }

    tmp<volScalarField> fAvg(new volScalarField(f));

    if (debug)
    {
        Pout<< "Call interpolateHeightProfileToPoints for int field" << endl;
    }

    // Interpolate from the height segments to cell centres
    interpolateHeightProfileToPoints
    (
        profile(), hMin, hMax, C, n, fAvg->primitiveFieldRef()
    );

    forAll(fAvg().boundaryField(), patchi)
    {
        if (debug)
        {
            Pout<< "Call interpolateHeightProfileToPoints for boundary field "
                << patchi << " " << mesh.boundary()[patchi].name() << " "
                << mesh.boundary()[patchi].type() << " "
                << fAvg->boundaryField()[patchi].type() << endl;
        }

        interpolateHeightProfileToPoints
        (
            profile,
            hMin,
            hMax,
            mesh.boundary()[patchi].Cf(),
            n,
            fAvg->boundaryFieldRef()[patchi]
        );
    }

    if (debug)
    {
        Pout<< "Sync boundary values" << endl;
    }

    // Sync processor and coupled values
    fAvg->correctBoundaryConditions();

    if (debug)
    {
        Pout<< "Exiting planarAverage(volScalarField& ...)" << endl;
    }

    return fAvg;
}

// ************************************************************************* //
