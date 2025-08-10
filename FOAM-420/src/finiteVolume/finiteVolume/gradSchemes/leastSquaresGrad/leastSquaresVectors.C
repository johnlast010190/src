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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/leastSquaresGrad/leastSquaresVectors.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::leastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresVectors>(mesh),
    pVectorsPtr_(nullptr),
    nVectorsPtr_(nullptr)
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresVectors::calcLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresVectors::calcLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

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
        dimensionedVector("zero", dimless/dimLength, Zero)
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
        dimensionedVector("zero", dimless/dimLength, Zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const surfaceScalarField& w = mesh_.weights();
//     const surfaceScalarField& magSf = mesh_.magSf();


    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), Zero);

    tmp<surfaceVectorField> td = mesh_.delta();
    const surfaceVectorField& d = td();
    tmp<surfaceVectorField> tmd = fvc::applyFaceMask(d);
    const surfaceVectorField& md = tmd();
    surfaceSymmTensorField mdSqr(sqr(d));
    fvc::applyFaceMaskTo(mdSqr);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        symmTensor wdd = (1.0/magSqr(d[facei]))*mdSqr[facei];

        dd[own] += wdd;
        dd[nei] += wdd;
    }


    surfaceVectorField::Boundary& pVectorsBf =
        lsP.boundaryFieldRef();

    forAll(pVectorsBf, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong because of incorrect weighting.  HJ, 23/Oct/2008
//         const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        const vectorField& pd = d.boundaryField()[patchi];
        const symmTensorField& pmdSqr = mdSqr.boundaryField()[patchi];

        forAll(pd, patchFacei)
        {
            dd[faceCells[patchFacei]] +=
                (1.0/magSqr(pd[patchFacei]))*pmdSqr[patchFacei];
        }
    }


    // Invert the dd tensor
    const symmTensorField invDd(inv(dd));
    // Fix: householder inverse.  HJ, 3/Nov/2009
    //symmTensorField invDd = hinv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar magSfByMagSqrd = 1.0/magSqr(d[facei]);

        lsP[facei] = magSfByMagSqrd*(invDd[own] & md[facei]);
        lsN[facei] = -magSfByMagSqrd*(invDd[nei] & md[facei]);
    }

    forAll(pVectorsBf, patchi)
    {
        fvsPatchVectorField& patchLsP = pVectorsBf[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        // Note: least squares in 1.4.1 and other OpenCFD versions
        // are wrong because of incorrect weighting.  HJ, 23/Oct/2008
//         const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        const vectorField& pd = d.boundaryField()[patchi];
        const vectorField& pmd = md.boundaryField()[patchi];

        forAll(pd, patchFacei)
        {
            patchLsP[patchFacei] =
                (1.0/magSqr(pd[patchFacei]))
                *(invDd[faceCells[patchFacei]] & pmd[patchFacei]);
        }
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing least square gradient vectors" << endl;
    }
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        calcLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        calcLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::leastSquaresVectors::movePoints()
{
    if (debug)
    {
        InfoIn("bool leastSquaresVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}



// ************************************************************************* //
