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

#include "finiteVolume/gradSchemes/faceLeastSquaresGrad/faceLeastSquaresGrad.H"
#include "finiteVolume/gradSchemes/faceLeastSquaresGrad/faceLeastSquaresVectors.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fields/Fields/zeroField/zeroField.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchField.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::initWithZero
(
    Field<Type>& v
)
{
    forAll(v, i)
    {
        v[i] = Type::zero;
    }
}

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::faceLeastSquaresGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tlsGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vsf.instance(),
                vsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad.ref();

    // Get reference to least square vectors
    const faceLeastSquaresVectors& lsv = faceLeastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
      = tinterpScheme_().interpolate(vsf);
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf.ref();

    const surfaceScalarField& weight = mesh.magSf();
    tmp<surfaceScalarField> tmWeight = fvc::applyFaceMask(weight);
    const surfaceScalarField& mWeight = tmWeight();

    // Calculate face area weighted average of edge values

    Field<Type> mean(mesh.nCells());
    initWithZero(mean);

    forAll(own, facei)
    {
        label ownFacei = own[facei];
        label neiFacei = nei[facei];

        mean[ownFacei] += mWeight[facei]*ssf[facei];
        mean[neiFacei] += mWeight[facei]*ssf[facei];
    }

    forAll(ssf.boundaryField(), patchi)
    {
        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
        const labelUList& faceCells = pssf.patch().faceCells();

        const fvsPatchScalarField& pmWeight = mWeight.boundaryField()[patchi];

        forAll(pssf, patchFacei)
        {
            mean[faceCells[patchFacei]] += pmWeight[patchFacei]*pssf[patchFacei];
        }
    }
    mean /= fvc::surfaceSum(weight);

    // Calculate gradient

    forAll(own, facei)
    {
        label ownFaceI = own[facei];
        label neiFaceI = nei[facei];

        lsGrad[ownFaceI] += ownLs[facei]*(ssf[facei]-mean[ownFaceI]);
        lsGrad[neiFaceI] += neiLs[facei]*(ssf[facei]-mean[neiFaceI]);
    }

    // Boundary faces
    forAll(ssf.boundaryField(), patchi)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];
        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        const labelUList& faceCells =
            lsGrad.boundaryField()[patchi].patch().faceCells();

        forAll(pssf, patchFaceI)
        {
            lsGrad[faceCells[patchFaceI]] +=
                patchOwnLs[patchFaceI]
               *(
                    pssf[patchFaceI]-mean[faceCells[patchFaceI]]
                );
        }

    }

    lsGrad.correctBoundaryConditions();

    return tlsGrad;
}


// ************************************************************************* //
