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

#include "finiteVolume/gradSchemes/phasicGaussGrad/phasicGaussGrad.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchField.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::phasicGaussGrad<Type>::phasicFaceWeights
(
    const scalar alpha_own,
    const scalar alpha_nei,
    const scalar alphaCutoff,
    scalar& w_cov,
    scalar& w_own,
    scalar& w_nei
)
{
    //for every face there are 3 relevant fraction
    // - covered:
    scalar alpha_cov = min(alpha_own, alpha_nei);
    // - uncovered:
    scalar alpha_unc = max(alpha_own, alpha_nei) - alpha_cov;
    // - empty: (undefined contribution, cancelled by reducing volume)
    //scalar alpha_emp = 1 - max(alpha_own, alpha_nei)


    //face gradient contributions in phasic space are asymmetric in normal
    //space, since we ignore out-of-phase contributions.
    //Every face value contribution to every cell will thus be made up
    //of a combination of linear interpolation for the covered fraction
    //and extrapolation for the uncovered section of the larger cell
    //
    //Cells below the cutoff will be treated as if they have zero alpha
    //Dont forget to scale the cell volume with the phase fraction at the
    //end or the integration will be wrong!!

    //both cell values need to be above cutoff for symmetric contribution
    //to be valid
    w_cov =
        (
            (alpha_cov > alphaCutoff)
            ? scalar(1) //condition true
            : scalar(0) //condition false
        )
      * alpha_cov;

    //remember to replace the below-cutoff covered face portion with
    //zero gradient to ensure consistency (inline below)

    w_own =
        //face area scaling
        (
            (alpha_own > alpha_nei)
            ?
            (
                alpha_unc //add covered area if covered is small (see above)
                + ((alpha_cov <= alphaCutoff) ? alpha_cov : scalar(0))
            )
            : scalar(0)
        )
        //make 100% sure out-of-phase gradients are identically zero
        *
        (
            (alpha_own > alphaCutoff)
            ? scalar(1)
            : scalar(0)
        );


    w_nei =
        //face area scaling
        (
            (alpha_nei > alpha_own)
            ?
            (
                alpha_unc
                //add covered area if covered is small (see above)
                + ((alpha_cov <= alphaCutoff) ? alpha_cov : scalar(0))
            )
            : scalar(0)
        )
        //make 100% sure out-of-phase gradients are identically zero
        *
        (
            (alpha_nei > alphaCutoff)
            ? scalar(1)
            : scalar(0)
        );
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
Foam::fv::phasicGaussGrad<Type>::gradFinternal
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name,
    const scalar alphaCutoff,
    const word alphaName
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
                vsf.db(),
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
    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const vectorField& mSf = tmSf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    const volScalarField& alpha
    (
        vsf.db().objectRegistry::template lookupObject<volScalarField>
        (
            alphaName
        )
    );

    //phasic gradient implementation
    //internal faces
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        scalar alpha_own = alpha[own];
        scalar alpha_nei = alpha[nei];

        //define weights
        scalar w_cov(0); //overlapping
        scalar w_own(0); //non-overlapping owner
        scalar w_nei(0); //non-overlapping neighbour

        //calculate phase weights
        phasicFaceWeights
            (alpha_own, alpha_nei, alphaCutoff, w_cov, w_own, w_nei);

        //gradient contributions
        GradType Sfssf = w_cov * mSf[facei]*issf[facei];
        GradType Sfssf_o = w_own * mSf[facei] * vsf[own];
        GradType Sfssf_n = w_nei * mSf[facei] * vsf[nei];

        // add contributions to cell gradients
        igGrad[own] += (Sfssf + Sfssf_o);
        igGrad[nei] -= (Sfssf + Sfssf_n);
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
Foam::fv::phasicGaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
)
{
    FatalErrorInFunction
        << " does not support gradf function"
        << exit(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gradFinternal(ssf, vsf, name, 1e-5, "alpha.foam")
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    Field<GradType>& igGrad = gGrad;

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch(mesh.boundary()[patchi]);
        if (isA<indirectPolyPatch>(patch.patch()))
        {
            FatalErrorInFunction
                << " does not support indirectPolyPatch (GIB)"
                << exit(FatalError);
        }

        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        fvsPatchField<Type> pssftmp = pssf;
        pssftmp = pssf.gradientBoundaryValue();

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssftmp[facei];
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
Foam::fv::phasicGaussGrad<Type>::calcGrad
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
        gradFinternal(tssf(), vsf, name, alphaCutoff_, alphaName_)
    );

    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    Field<GradType>& igGrad = gGrad;


    //apply boundary contributions
    const volScalarField& alpha
    (
        vsf.db().objectRegistry::template lookupObject<volScalarField>(alphaName_)
    );

    // treatment for extrapolated boundaries
    //autoPtr<tensorField> gradScalePtr(nullptr);

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    forAll(vsf.boundaryField(), patchi)
    {
        const fvPatch& patch(mesh.boundary()[patchi]);

        const labelUList& pFaceCells = patch.faceCells();
        const vectorField& mSfp = mSf.boundaryField()[patchi];

        const fvPatchField<Type>& vsfp(vsf.boundaryField()[patchi]);
        tmp<Field<Type>> vsfp_own = vsfp.patchInternalField();

        const fvPatchField<scalar>& alphap(alpha.boundaryField()[patchi]);
        tmp<scalarField> alpha_own = alphap.patchInternalField();

        if (vsfp.extrapolated())
        {
            FatalErrorInFunction << this->typeName
                << " does not support extrapolated"
                << exit(FatalError);
        }
        else if (vsfp.coupled())
        {
            tmp<scalarField> alpha_nei = alphap.patchNeighbourField();

            const fvsPatchField<Type>& pssf = tssf->boundaryField()[patchi];

            forAll(patch, facei)
            {
                //define weights
                scalar w_cov(0); //overlapping
                scalar w_own(0); //non-overlapping owner
                scalar w_nei(0); //non-overlapping neighbour

                //calculate phase weights
                phasicFaceWeights
                (
                    alpha_own().operator[](facei),
                    alpha_nei().operator[](facei),
                    alphaCutoff_,
                    w_cov,
                    w_own,
                    w_nei
                );

                //gradient contributions
                GradType Sfssf = w_cov * mSfp[facei] * pssf[facei];
                GradType Sfssf_o =
                    w_own * mSfp[facei] * vsfp_own().operator[](facei);

                igGrad[pFaceCells[facei]] += (Sfssf + Sfssf_o);
            }
        }
        else
        {
            tmp<Field<Type>> gbvTmp = vsfp.gradientBoundaryValue();

            forAll(patch, facei)
            {
                //define weights
                scalar w_cov(0); //overlapping
                scalar w_own(0); //non-overlapping owner
                scalar w_nei(0); //non-overlapping neighbour

                //calculate phase weights
                phasicFaceWeights
                (
                    alpha_own().operator[](facei),
                    alphap[facei],
                    alphaCutoff_,
                    w_cov,
                    w_own,
                    w_nei
                );

                //gradient contributions
                GradType Sfssf = w_cov * mSfp[facei] * gbvTmp()[facei];
                GradType Sfssf_o =
                    w_own * mSfp[facei] * vsfp_own().operator[](facei);

                igGrad[pFaceCells[facei]] += (Sfssf + Sfssf_o);
            }
        }
    }

    //since multiple extrapolation boundaries can contribute to the same cell
    //the multiplication has top be done after the summation is complete
    /* this is for extrapolated boundaries - disable for now
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
    */

    //phasic 4-volume is proportional to phase fraction
    const scalarField& Vol(mesh.V());
    igGrad /= (Vol*max(alphaCutoff_, alpha));

    gGrad.correctBoundaryConditions();
    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void Foam::fv::phasicGaussGrad<Type>::correctBoundaryConditions
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
