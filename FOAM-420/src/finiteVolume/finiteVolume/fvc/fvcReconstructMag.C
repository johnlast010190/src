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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcReconstruct.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> reconstructMag(const surfaceScalarField& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(Sf);
    const surfaceVectorField& mSf = tmSf();

    tmp<volScalarField> treconField
    (
        new volScalarField
        (
            IOobject
            (
                "reconstruct("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "0",
                ssf.dimensions()/dimArea,
                scalar(0)
            ),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField& rf = treconField.ref();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        rf[own] += (mSf[facei] & (Cf[facei] - C[own]))*ssf[facei]/magSf[facei];
        rf[nei] -= (mSf[facei] & (Cf[facei] - C[nei]))*ssf[facei]/magSf[facei];
    }

    const surfaceScalarField::Boundary& bsf = ssf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvsPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];
        const vectorField& pmSf = mSf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];
            rf[own] +=
                (pmSf[pFacei] & (pCf[pFacei] - C[own]))
               *psf[pFacei]/pMagSf[pFacei];
        }
    }

    rf /= mesh.V();

    treconField.ref().correctBoundaryConditions();

    return treconField;
}


tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>& tssf)
{
    tmp<volScalarField> tvf
    (
        fvc::reconstructMag(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
