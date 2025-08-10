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
    (c) 2010-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchField.H"
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
Foam::fv::gaussGrad<Type>::gradFinternal
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                ssf.instance(),
                ssf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimLength,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName,
            ssf.boundaryField().patchTypes()
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(Sf);
    const vectorField& mSf = tmSf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, facei)
    {
        GradType Sfssf = mSf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

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
Foam::fv::gaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gradFinternal(ssf, name)
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

    igGrad /= mesh.V();

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
Foam::fv::gaussGrad<Type>::calcGrad
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
        gradFinternal(tssf(), name)
    );

    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    Field<GradType>& igGrad = gGrad;

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();


    //apply boundary contributions

    // treatment for extrapolated boundaries
    autoPtr<tensorField> gradScalePtr(nullptr);
    const scalarField& Vol(mesh.V());


    forAll(vsf.boundaryField(), patchi)
    {
        const fvPatch& patch(mesh.boundary()[patchi]);

        const labelUList& pFaceCells = patch.faceCells();

        const vectorField& pmSf = mSf.boundaryField()[patchi];

        if (vsf.boundaryField()[patchi].extrapolated())
        {
            if (!gradScalePtr.valid())
            {
                gradScalePtr.reset(new tensorField(vsf.size(), tensor::I));
            }

            tmp<vectorField> delta(patch.Cf() - patch.Cn());

            forAll(patch, facei)
            {
                label celli = pFaceCells[facei];

                //scaling contribution
                gradScalePtr->operator[](celli)
                    -= (pmSf[facei]*delta->operator[](facei))/Vol[celli];

                //zero gradient contributions
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
            tmp<Field<Type>> gbv =
                vsf.boundaryField()[patchi].gradientBoundaryValue();

            forAll(patch, facei)
            {
                igGrad[pFaceCells[facei]]
                    += pmSf[facei]*gbv().operator[](facei);
            }
        }
    }

    //since multiple extrapolation boundaries can contribute to the same cell
    //the multiplication has to be done after the summation is complete
    if (gradScalePtr.valid())
    {
        forAll(igGrad, celli)
        {
            if (gradScalePtr->operator[](celli) != tensor::I)
            {
                igGrad[celli]
                    = (inv(gradScalePtr->operator[](celli)) & igGrad[celli]);
            }
        }
    }

    igGrad /= Vol;

    gGrad.correctBoundaryConditions();
    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
Foam::tmp
<
    Foam::BlockLduSystem
    <
        Foam::vector, typename Foam::outerProduct<Foam::vector, Type>::type
    >
> Foam::fv::gaussGrad<Type>::fvmGrad
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


template<class Type>
void Foam::fv::gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    typename GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >::Boundary& gGradbf = gGrad.boundaryFieldRef();

    forAll(vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
     }
}


// ************************************************************************* //
