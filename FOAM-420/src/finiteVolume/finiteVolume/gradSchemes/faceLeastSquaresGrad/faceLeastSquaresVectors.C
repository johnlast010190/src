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
    (c) 2011 OpenFOAM Foundation
    (c) 2014 CSIR
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/faceLeastSquaresGrad/faceLeastSquaresVectors.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(faceLeastSquaresVectors, 0);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::faceLeastSquaresVectors::faceLeastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, MoveableMeshObject, faceLeastSquaresVectors>(mesh),
    pVectorsPtr_(nullptr),
    nVectorsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::faceLeastSquaresVectors::~faceLeastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "faceLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField& Sf = mesh.Sf();

    const surfaceScalarField& weight = mesh.magSf();
    tmp<surfaceScalarField> tmWeight = fvc::applyFaceMask(weight);
    const surfaceScalarField& mWeight = tmWeight();

    // Set up temporary storage for mean of edge centres and calculate

    vectorField meanPoints(mesh_.nCells(), vector::zero);
    scalarField weightSum = fvc::surfaceSum(weight)().internalField();
    surfaceVectorField Ce
    (
        IOobject
        (
            "Ce",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimLength
    );

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        Ce[facei] = w[facei]*C[own] + (1-w[facei])*C[nei];
        meanPoints[own] += mWeight[facei]*Ce[facei];
        meanPoints[nei] += mWeight[facei]*Ce[facei];
    }

    forAll(w.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        fvsPatchVectorField& pCe = Ce.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& pmWeight = mWeight.boundaryField()[patchi];
        const labelUList& faceCells = pw.patch().faceCells();

        tmp<vectorField> tpd = pw.patch().delta();
        vectorField& pd = tpd.ref();

        if (pw.coupled())
        {
            forAll(pw, patchFacei)
            {
                pCe[patchFacei] = C[faceCells[patchFacei]] + (1-pw[patchFacei])*pd[patchFacei];
                meanPoints[faceCells[patchFacei]] += pmWeight[patchFacei]*pCe[patchFacei];
            }
        }
        else
        {
            forAll(pw, patchFacei)
            {
                //pCe[patchFacei] = C[faceCells[patchFacei]] + pd[patchFacei]; // Should coincide with Cf on boundary
                // Boundary values are assumed specified at a point projected orthogonally to the boundary from the internal point
                // Should use deltaCoeff since this might be stabilised and not end up quite on the wall - for consistency
                // with snGrad
                pCe[patchFacei] = C[faceCells[patchFacei]] +
                                  ( Sf.boundaryField()[patchi][patchFacei] /
                                    (magSf.boundaryField()[patchi][patchFacei]*pw.patch().deltaCoeffs()[patchFacei]) );

                meanPoints[faceCells[patchFacei]] += pmWeight[patchFacei]*pCe[patchFacei];
            }
        }
    }
    meanPoints /= weightSum;

    // Set up temporary storage for the dd tensor (before inversion)

    tensorField dd(mesh_.nCells(), symmTensor::zero);

    // Division by norm factor here and below to prevent numerical issues with
    // inversion of tiny tensor
    scalarField norm =
        fvc::surfaceSum(weight*magSqr(mesh.delta()))->primitiveField();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        dd[own] += mWeight[facei]/norm[own]*Ce[facei]*(Ce[facei]-meanPoints[own]);
        dd[nei] += mWeight[facei]/norm[nei]*Ce[facei]*(Ce[facei]-meanPoints[nei]);
    }

    surfaceVectorField::Boundary& blsP = lsP.boundaryFieldRef();

    forAll(blsP, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchVectorField& pCe = Ce.boundaryField()[patchi];
        const fvsPatchScalarField& pmWeight = mWeight.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.patch().faceCells();

        forAll(pw, patchFacei)
        {
            dd[faceCells[patchFacei]] +=
                pmWeight[patchFacei]/norm[faceCells[patchFacei]]*pCe[patchFacei]*(pCe[patchFacei]-meanPoints[faceCells[patchFacei]]);
        }
    }


    // Invert the dd tensor
    // First stabilise in any singular directions
    vector diagStab = (Vector<label>::one-mesh_.geometricD())*scalar(0.5);
    dd += tensor(diagStab.x(), 0, 0, 0, diagStab.y(), 0, 0, 0, diagStab.z());
    const tensorField invDd(inv(dd));


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        lsP[facei] = mWeight[facei]/norm[own]*(invDd[own] & Ce[facei]);
        lsN[facei] = mWeight[facei]/norm[nei]*(invDd[nei] & Ce[facei]);
    }

    forAll(blsP, patchi)
    {
        fvsPatchVectorField& patchLsP = blsP[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pmWeight = mWeight.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        forAll(pw, patchFacei)
        {
            patchLsP[patchFacei] =
               pmWeight[patchFacei]/norm[faceCells[patchFacei]]*(invDd[faceCells[patchFacei]] & Ce.boundaryField()[patchi][patchFacei]);
        }
    }

    if (debug)
    {
        Info<< "faceLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField& Foam::faceLeastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::faceLeastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::faceLeastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}


// ************************************************************************* //
