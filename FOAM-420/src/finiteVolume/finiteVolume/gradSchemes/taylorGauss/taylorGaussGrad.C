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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/taylorGauss/taylorGaussGrad.H"
#include "finiteVolume/gradSchemes/taylorGauss/taylorGaussData.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchField.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"

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
Foam::fv::taylorGaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gaussGrad<Type>::gradFinternal(ssf, name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    Field<GradType>& igGrad = gGrad;

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pmSf = mSf.boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
        fvsPatchField<Type> pssftmp = pssf;
        pssftmp = pssf.gradientBoundaryValue();

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pmSf[facei]*pssftmp[facei];
        }
    }

    // Get reference to skew face vectors
    const taylorGaussData& data = taylorGaussData::New(mesh);
    const tensorField& invPseudoVol = data.invPseudoVolumes();
    igGrad = invPseudoVol & igGrad;

    gGrad.correctBoundaryConditions();

    return tgGrad;
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
Foam::fv::taylorGaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        tinterpScheme_().interpolate(vsf)
    );

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gaussGrad<Type>::gradFinternal(tssf(), name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    Field<GradType>& igGrad = gGrad;

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    //apply boundary contributions
    forAll(vsf.boundaryField(), patchi)
    {
        const fvPatch& patch(vsf.mesh().boundary()[patchi]);

        const labelUList& pFaceCells = patch.faceCells();

        const vectorField& pmSf = mSf.boundaryField()[patchi];

        if (vsf.boundaryField()[patchi].extrapolated())
        {
            const tmp<vectorField> delta(patch.delta());

            forAll(patch, facei)
            {
                label celli = pFaceCells[facei];

                //zero gradient contributions - correction in tensor
                igGrad[celli] += (pmSf[facei]*vsf[celli]);
            }
        }
        else if (vsf.boundaryField()[patchi].coupled())
        {
            const fvsPatchField<Type>& pssf = tssf->boundaryField()[patchi];

            forAll(patch, facei)
            {
                igGrad[pFaceCells[facei]] += pmSf[facei]*pssf[facei];
            }
        }
        else
        {
            tmp<Field<Type>> gbv;
            if (isA<indirectPolyPatch>(patch.patch()))
            {
                const fvsPatchField<Type>& pssf = tssf->boundaryField()[patchi];
                gbv = pssf.gradientBoundaryValue().ptr();
            }
            else
            {
                gbv = vsf.boundaryField()[patchi].gradientBoundaryValue().ptr();
            }

            forAll(patch, facei)
            {
                igGrad[pFaceCells[facei]]
                    += pmSf[facei]*gbv->operator[](facei);
            }
        }
    }

    // Get reference to skew face vectors
    const taylorGaussData& data = taylorGaussData::New(mesh);
    const tensorField& invPseudoVol = data.invPseudoVolumes(vsf);

    igGrad = (invPseudoVol & igGrad);

    // Replace boundary surface normal component with snGrad
    gGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
Foam::tmp
<
    Foam::BlockLduSystem
    <
        Foam::vector, typename Foam::outerProduct<Foam::vector, Type>::type
    >
> Foam::fv::taylorGaussGrad<Type>::fvmGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorInFunction
        << "Implicit gradient operator defined only for scalar."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType>> tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
}


// ************************************************************************* //
