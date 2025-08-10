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

#include "finiteVolume/gradSchemes/twistCorrectedGauss/twistCorrectedGaussGrad.H"
#include "finiteVolume/gradSchemes/twistCorrectedGauss/twistCorrectedGaussData.H"
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
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::twistCorrectedGaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        tinterpScheme_().interpolate(vsf)
    );
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf();

    //modification to add boundary influence to gradient calculation
    // allows difference between explicit terms and implicit terms
    forAll(tssf->boundaryField(), pi)
    {
        if (!tssf->boundaryField()[pi].coupled())
        {
            tssf->boundaryFieldRef()[pi]
                += vsf.boundaryField()[pi].gradientBoundaryValue()
                - Field<Type>(vsf.boundaryField()[pi]);
        }
    }

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    // Make the initial correction gradient
    tmp<GradFieldType> tcorrGrad;
    tcorrGrad =
        tgradCorrectionScheme_->grad(vsf, "grad("+vsf.name()+")Corrector");

    tmp<GradFieldType> tgGrad;

    for (label iter = 0; iter < nIterations_; iter++)
    {
        const GradFieldType& corrGrad = tcorrGrad();

        tgGrad = tmp<GradFieldType>
        (
            new GradFieldType
            (
                IOobject
                (
                    name,
                    ssf.instance(),
                    vsf.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<GradType>
                (
                    "0",
                    ssf.dimensions()/dimLength,
                    pTraits<GradType>::zero
                ),
                zeroGradientFvPatchField<GradType>::typeName,
                ssf.boundaryField().patchTypes()
            )
        );
        GradFieldType& gGrad = tgGrad.ref();

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        // Get reference to skew face vectors
        const twistCorrectedGaussData& data = twistCorrectedGaussData::New(mesh);
        const surfaceTensorField& skewCorrTensors = data.skewCorrTensors();

        const surfaceVectorField& Sf = mesh.Sf();
        tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(Sf);
        const surfaceVectorField& mSf = tmSf();

        const surfaceScalarField& w = mesh.weights();

        Field<GradType>& igGrad = gGrad;
        const Field<Type>& issf = ssf;

        forAll(owner, facei)
        {
            GradType Sfssf = mSf[facei]*issf[facei] +
                  w[facei]*(skewCorrTensors[facei] & corrGrad[owner[facei]]) +
                  (1-w[facei])*(skewCorrTensors[facei] & corrGrad[neighbour[facei]]);
            igGrad[owner[facei]] += Sfssf;
            igGrad[neighbour[facei]] -= Sfssf;
        }

        forAll(mesh.boundary(), patchi)
        {
            const labelUList& pFaceCells =
                mesh.boundary()[patchi].faceCells();

            const vectorField& pmSf = mSf.boundaryField()[patchi];
            const tensorField& pSkewCorrTensors = skewCorrTensors.boundaryField()[patchi];
            const scalarField& pw = w.boundaryField()[patchi];

            const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

            if (mesh.boundary()[patchi].coupled())
            {
                const Field<GradType> neiCorrGrad
                (
                    corrGrad.boundaryField()[patchi].patchNeighbourField()
                );
                forAll(mesh.boundary()[patchi], facei)
                {
                    igGrad[pFaceCells[facei]] += pmSf[facei]*pssf[facei] +
                          pw[facei]*(pSkewCorrTensors[facei] & corrGrad[pFaceCells[facei]]) +
                          (1-pw[facei])*(pSkewCorrTensors[facei] & neiCorrGrad[facei]);
                }
            }
            else
            {
                forAll(mesh.boundary()[patchi], facei)
                {
                    igGrad[pFaceCells[facei]] += pmSf[facei]*pssf[facei] +
                          (pSkewCorrTensors[facei] & corrGrad[pFaceCells[facei]]);
                }
            }
        }

        igGrad /= data.consistentVols();

        gGrad.correctBoundaryConditions();

        // Replace boundary surface normal component with snGrad
        gaussGrad<Type>::correctBoundaryConditions(vsf, gGrad);

        tcorrGrad = tgGrad; //Transfer to corrector for subsequent iteration
    }

    return tcorrGrad;
}


// ************************************************************************* //
