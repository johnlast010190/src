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
    (c) 2018-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "GRF.H"
#include "fvMatrices/fvMatrices.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/faceSet.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::GRF::init()
{
    coorFramePtr_ = coordinateFrame::lookupNew(mesh_, GRFDict_);
    coorFrame().registry().attachObject
    (
        this->name(),
        frameSourceFaces::typeName
    );

    MRF_ = true;
    fieldNames_.setSize(0);
    applied_.setSize(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::GRF::GRF
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    frameSources(name, modelType, dict, mesh),
    GRFDict_(dict)
{
    init();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::GRF::applyRotationalFlux(surfaceScalarField& phiOmega) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();
    const vectorField& Cfi = Cf;
    const vectorField& Sfi = Sf;
    const labelList& internalFaces = frameSourceFaces_.internalFaces();
    forAll(internalFaces, fI)
    {
        label fL = internalFaces[fI];
        face f  = mesh_.faces()[fL];
        label nPoints = f.size();
        if (nPoints > 3 && getFrameSwitches().conservative())
        {
            phiOmega[fL] += getTriangularFlux(f);
        }
        else
        {
            phiOmega[fL] +=
                coorFrame().frameVelocity(Cfi[fL], false)&Sfi[fL];
        }
    }

    const labelListList& includedFaces = frameSourceFaces_.includedFaces();
    forAll(includedFaces, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        forAll(includedFaces[patchi], i)
        {
            label patchFacei = includedFaces[patchi][i];
            face f = mesh_.faces()[patchFacei + pp.start()];
            label nPoints = f.size();
            if (nPoints > 3 && getFrameSwitches().conservative())
            {
                phiOmega.boundaryFieldRef()[patchi][patchFacei] +=
                    getTriangularFlux(f);
            }
            else
            {
                phiOmega.boundaryFieldRef()[patchi][patchFacei] +=
                    (
                        coorFrame().frameVelocity
                        (
                            Cf.boundaryField()[patchi][patchFacei],
                            false
                        )
                    )
                  & Sf.boundaryField()[patchi][patchFacei];
            }
        }
    }
    const labelListList& excludedFaces = frameSourceFaces_.excludedFaces();
    forAll(excludedFaces, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        forAll(excludedFaces[patchi], i)
        {
            label patchFacei = excludedFaces[patchi][i];
            face f = mesh_.faces()[patchFacei + pp.start()];
            label nPoints = f.size();
            if (nPoints > 3 && getFrameSwitches().conservative())
            {
                phiOmega.boundaryFieldRef()[patchi][patchFacei] +=
                    getTriangularFlux(f);
            }
            else
            {
                phiOmega.boundaryFieldRef()[patchi][patchFacei] +=
                    (
                        coorFrame().frameVelocity
                        (
                            Cf.boundaryField()[patchi][patchFacei],
                            false
                        )
                    )
                  & Sf.boundaryField()[patchi][patchFacei];
            }
        }
    }
}


Foam::scalar Foam::fv::GRF::getTriangularFlux(const face& f) const
{
    const pointField& p = mesh_.points();
    scalar totalRotFlux = 0;
    point fCentre = p[f[0]];
    label nPoints = f.size();
    for (label pi = 1; pi < nPoints; pi++)
    {
        fCentre += p[f[pi]];
    }
    fCentre /= nPoints;
    for (label pi = 0; pi < nPoints; pi++)
    {
        const point& nextPoint = p[f[(pi + 1) % nPoints]];
        vector c = p[f[pi]] + nextPoint + fCentre;
        c /= 3.;
        vector sf = 0.5*(nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
        totalRotFlux += coorFrame().frameVelocity(c, false) & sf;
    }
    return totalRotFlux;
}


void Foam::fv::GRF::relaxBoundary
(
    fvBlockMatrix<vector>& UEqn,
    const scalar relax
) const
{
    CoeffField<vector>& diagCoeff = UEqn.diag();
    label nCoeffs = pTraits<vector>::nComponents;

    tensorField& bDiag = diagCoeff.asSquare();
    const vectorField& psiInt = UEqn.psi();

    vectorField& source = UEqn.source();

    const scalar oneMrelax = 1 - relax;
    const scalar twoMrelax = 2 - relax;

    boolList markCells(psiInt.size(), false);

    const labelListList& includedFaces = frameSourceFaces_.includedFaces();

    forAll(includedFaces, patchi)
    {
        const fvPatch& fvp =  mesh_.boundary()[patchi];
        const labelUList& faceCells = fvp.faceCells();
        forAll(includedFaces[patchi], i)
        {
            label patchFacei = includedFaces[patchi][i];
            const label cI = faceCells[patchFacei];

            if (!markCells[cI])
            {
                for (label iC=0; iC<nCoeffs; iC++)
                {
                    source[cI][iC] +=
                        bDiag[cI][iC*nCoeffs+iC]*psiInt[cI][iC]*oneMrelax;
                }
                markCells[cI] = true;
            }
        }
    }

    markCells = false;
    forAll(includedFaces, patchi)
    {
        const fvPatch& fvp =  mesh_.boundary()[patchi];
        const labelUList& faceCells = fvp.faceCells();
        forAll(includedFaces[patchi], i)
        {
            label patchFacei = includedFaces[patchi][i];
            const label cI = faceCells[patchFacei];

            if (!markCells[cI])
            {
                for (label iC=0; iC<nCoeffs; iC++)
                {
                    bDiag[cI][iC*nCoeffs+iC] *= twoMrelax;
                }
                markCells[cI] = true;
            }
        }
    }
}


void Foam::fv::GRF::makeRelative(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        label celli = cells_[i];
        U[celli] -= coorFrame().frameVelocity(C[celli], false);
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();
    const labelListList& includedFaces = frameSourceFaces_.includedFaces();
    forAll(includedFaces, patchi)
    {
        forAll(includedFaces[patchi], i)
        {
            label patchFacei = includedFaces[patchi][i];
            Ubf[patchi][patchFacei] -=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }

    const labelListList& excludedFaces = frameSourceFaces_.excludedFaces();
    forAll(excludedFaces, patchi)
    {
        forAll(excludedFaces[patchi], i)
        {
            label patchFacei = excludedFaces[patchi][i];
            Ubf[patchi][patchFacei] -=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }
}


void Foam::fv::GRF::zero
(
    surfaceScalarField& phi
) const
{
    if (!active_ || getFrameSwitches().ddtPhiCorr())
    {
        return;
    }

    Field<scalar>& phii = phi.primitiveFieldRef();

    forAll(frameSourceFaces_.internalFaces(), i)
    {
        phii[frameSourceFaces_.internalFaces()[i]] = Zero;
    }

    typename GeometricField<scalar, fvsPatchField, surfaceMesh>::Boundary&
    phibf = phi.boundaryFieldRef();

    const labelListList& includedFaces = frameSourceFaces_.includedFaces();
    forAll(includedFaces, patchi)
    {
        forAll(includedFaces[patchi], i)
        {
            phibf[patchi][includedFaces[patchi][i]] = Zero;
        }
    }
    const labelListList& excludedFaces = frameSourceFaces_.excludedFaces();
    forAll(excludedFaces, patchi)
    {
        forAll(excludedFaces[patchi], i)
        {
            phibf[patchi][excludedFaces[patchi][i]] = Zero;
        }
    }
}


bool Foam::fv::GRF::read()
{
    // Use this instead of on-construction initialisation
    if (cellSetOption::read(GRFDict_))
    {
        init();

        return true;
    }
    else
    {
        return false;
    }
    return false;
}


// ************************************************************************* //
