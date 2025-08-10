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
    (c) 2016-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/taylorGauss/taylorGaussData.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(taylorGaussData, 0);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::taylorGaussData::taylorGaussData(const fvMesh& mesh)
:
    MeshObject<fvMesh, MoveableMeshObject, taylorGaussData>(mesh),
    invPseudoVolPtr_(nullptr),
    fieldInvPseudoVolPtrs_(new HashPtrTable<Field<tensor>, word, string::hash>())
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::taylorGaussData::~taylorGaussData()
{
    deleteDemandDrivenData(invPseudoVolPtr_);
    deleteDemandDrivenData(fieldInvPseudoVolPtrs_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::tensorField> Foam::taylorGaussData::pseudoVolumes() const
{
    if (debug)
    {
        Info<< "taylorGaussNormals::makeInvPseudoVolumes() :"
            << "Constructing inverse pseudo-volume tensors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceVectorField& Sf = mesh.Sf();

    // Set up temporary storage for edge centres and calculate
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

    tmp<tensorField> pseudoVolPtr
    (
        new tensorField(mesh_.nCells(), tensor::zero)
    );

    tensorField& pseudoVol(pseudoVolPtr.ref());

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(Sf);
    const surfaceVectorField& mSf = tmSf();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        Ce[facei] = w[facei]*C[own] + (1-w[facei])*C[nei];
        pseudoVol[own] += mSf[facei]*(Ce[facei]-C[own]);
        pseudoVol[nei] -= mSf[facei]*(Ce[facei]-C[nei]);
    }

    forAll(w.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        fvsPatchVectorField& pCe = Ce.boundaryFieldRef()[patchi];
        const fvsPatchVectorField& pmSf = mSf.boundaryField()[patchi];
        const labelUList& faceCells = pw.patch().faceCells();

        tmp<vectorField> tpd = pw.patch().delta();
        const vectorField& pd = tpd();

        if (pw.coupled())
        {
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];
                pCe[patchFacei] = C[own]+(1-pw[patchFacei])*pd[patchFacei];
                pseudoVol[own] += pmSf[patchFacei]*(pCe[patchFacei]-C[own]);
            }
        }
        else
        {
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];

                // Boundary values are assumed specified at a point
                // projected orthogonally to the boundary from the internal
                // point
                pCe[patchFacei] = C[own] + pd[patchFacei];
                pseudoVol[own] += pmSf[patchFacei]*(pCe[patchFacei]-C[own]);
            }
        }
    }

    mesh.stabiliseEmptyDirections(pseudoVol);

    if (debug)
    {
        Info<< "taylorGaussNormals::makeCoeffs() :"
            << "Finished constructing skew Gauss data"
            << endl;
    }

    return pseudoVolPtr;
}


template<class Type>
Foam::tmp<Foam::tensorField> Foam::taylorGaussData::pseudoVolumes
(
    GeometricField<Type, fvPatchField, volMesh>  vf
) const
{
    const fvMesh& mesh(vf.mesh());

    // treatment for extrapolated boundaries
    tmp<tensorField> pseudoVolPtr(pseudoVolumes());
    tensorField& pseudoVol(pseudoVolPtr.ref());

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh_.Sf());
    const surfaceVectorField& mSf = tmSf();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatch& patch(mesh.boundary()[patchi]);

        const labelUList& pFaceCells = patch.faceCells();

        const vectorField& pmSf = mSf.boundaryField()[patchi];

        if (vf.boundaryField()[patchi].extrapolated())
        {
            tmp<vectorField> delta(patch.delta());

            forAll(patch, facei)
            {
                pseudoVol[pFaceCells[facei]]
                    -= (pmSf[facei]*delta->operator[](facei));
            }
        }
    }

    mesh.stabiliseEmptyDirections(pseudoVol);

    return pseudoVolPtr;
}


const Foam::tensorField& Foam::taylorGaussData::invPseudoVolumes() const
{
    if (!invPseudoVolPtr_)
    {
        // Store inverse of tensor
        invPseudoVolPtr_ = inv(pseudoVolumes()).ptr();
    }

    return *invPseudoVolPtr_;
}

const Foam::tensorField& Foam::taylorGaussData::invPseudoVolumes
(
    const volScalarField& vsf
) const
{
    // check if field specific modification as required
    bool extrapolated = false;
    forAll(vsf.boundaryField(), patchi)
    {
        if (vsf.boundaryField()[patchi].extrapolated())
        {
            extrapolated = true;
            break;
        }
    }

    if (extrapolated)
    {
        if (!fieldInvPseudoVolPtrs_->found(vsf.name()))
        {
            fieldInvPseudoVolPtrs_->insert
            (
                vsf.name(),
                (inv(pseudoVolumes(vsf))).ptr()
            );
        }

        return *(fieldInvPseudoVolPtrs_->operator[](vsf.name()));
    }

    return invPseudoVolumes();
}

const Foam::tensorField& Foam::taylorGaussData::invPseudoVolumes
(
    const volVectorField& vvf
) const
{
    // check if field specific modification as required
    bool extrapolated = false;
    forAll(vvf.boundaryField(), patchi)
    {
        if (vvf.boundaryField()[patchi].extrapolated())
        {
            extrapolated = true;
            break;
        }
    }

    if (extrapolated)
    {
        if (!fieldInvPseudoVolPtrs_->found(vvf.name()))
        {
            fieldInvPseudoVolPtrs_->insert
            (
                vvf.name(),
                (inv(pseudoVolumes(vvf))).ptr()
            );
        }

        return *(fieldInvPseudoVolPtrs_->operator[](vvf.name()));
    }

    return invPseudoVolumes();
}


bool Foam::taylorGaussData::movePoints()
{
    deleteDemandDrivenData(invPseudoVolPtr_);
    deleteDemandDrivenData(fieldInvPseudoVolPtrs_);

    return true;
}


// ************************************************************************* //
