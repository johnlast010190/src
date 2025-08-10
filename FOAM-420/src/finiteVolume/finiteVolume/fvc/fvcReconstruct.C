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
    (c) 2017-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcReconstruct.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
reconstruct
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());

    tmp<volTensorField> ttf = surfaceSum(SfHat*mesh.Sf());
    mesh.stabiliseEmptyDirections(ttf.ref());

    tmp<GeometricField<GradType, fvPatchField, volMesh>> treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(ttf)&surfaceSum(SfHat*ssf),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );

    treconField.ref().correctBoundaryConditions();

    return treconField;

}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
reconstruct
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh>> tvf
    (
        fvc::reconstruct(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
reconstructV
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tmssf =
        fvc::applyFaceMask(ssf);
    const GeometricField<Type, fvsPatchField, surfaceMesh>& mssf = tmssf();
    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0", ssf.dimensions(), pTraits<GradType>::zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );

    GeometricField<GradType, fvPatchField, volMesh>& rf = treconField();

    Field<tensor> faceVol(rf.size(), tensor::zero);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vector fDeltaO(mesh.Cf()[facei] - mesh.C()[owner[facei]]);
        vector fDeltaN(mesh.C()[neighbour[facei]] - mesh.Cf()[facei]);

        faceVol[owner[facei]] += fDeltaO * mSf[facei];
        faceVol[neighbour[facei]] += fDeltaN * mSf[facei];

        rf[owner[facei]] += fDeltaO*mssf[facei];
        rf[neighbour[facei]] += fDeltaN*mssf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& faceCells(mesh.boundary()[patchi].faceCells());

        forAll(faceCells, facei)
        {
            label celli = faceCells[facei];
            vector fDeltaO
            (
                mesh.boundary()[patchi].Cf()[facei] - mesh.C()[celli]
            );

            faceVol[celli] += fDeltaO*mSf.boundaryField()[patchi][facei];
            rf[celli] += fDeltaO*mssf.boundaryField()[patchi][facei];
        }
    }

    rf.internalField() = (inv(faceVol) & rf.internalField());

    rf.correctBoundaryConditions();

    return treconField;

}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
reconstructV
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh>> tvf
    (
        fvc::reconstructV(tssf())
    );
    tssf.clear();
    return tvf;
}

template<class Type>
tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<surfaceVectorField> tmSf = fvc::applyFaceMask(mesh.Sf());
    const surfaceVectorField& mSf = tmSf();

    tmp<GeometricField<Type, fvPatchField, volMesh>> trsf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0", vsf.dimensions(), pTraits<Type>::zero
            ),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );

    GeometricField<Type, fvPatchField, volMesh>& rf = trsf.ref();

    Field<scalar> faceVol(rf.size(), 0.0);

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    scalar oneThird(1.0/3.0);
    scalar oneTwelfth(1.0/12.0);

    forAll(owner, facei)
    {
        scalar baseHeightO
        (
            (mesh.Cf()[facei] - mesh.C()[owner[facei]])& mSf[facei]
        );
        scalar baseHeightN
        (
            (mesh.C()[neighbour[facei]] - mesh.Cf()[facei])& mSf[facei]
        );

        faceVol[owner[facei]] += oneThird * baseHeightO;
        faceVol[neighbour[facei]] += oneThird * baseHeightN;

        rf[owner[facei]]
            += 0.25*baseHeightO*ssf[facei]
            + oneTwelfth*baseHeightO*vsf[owner[facei]];
        rf[neighbour[facei]]
            += 0.25*baseHeightN*ssf[facei]
            + oneTwelfth*baseHeightN*vsf[neighbour[facei]];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& faceCells(mesh.boundary()[patchi].faceCells());

        forAll(faceCells, facei)
        {
            label celli = faceCells[facei];
            scalar baseHeightO
            (
                (mesh.boundary()[patchi].Cf()[facei] - mesh.C()[celli])
                & mSf.boundaryField()[patchi][facei]
            );

            faceVol[celli] += oneThird * baseHeightO;
            rf[celli]
                += 0.25*baseHeightO*ssf.boundaryField()[patchi][facei]
                + oneTwelfth*baseHeightO*vsf[celli];
        }
    }

    rf.primitiveFieldRef() = (rf.primitiveField() / faceVol);

    rf.correctBoundaryConditions();

    return trsf;

}

template<class Type>
tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
(
    const GeometricField<Type, fvPatchField, volMesh> & vf,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tsf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> trvf
    (
        reconstructSmooth(vf, tsf())
    );
    tsf.clear();
    return trvf;
}


template<class Type>
tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return reconstructSmooth
    (
        vf,
        fvc::interpolate(vf, "interpolate(" + vf.name() + "_reconstruct)")
    );
}


template<class Type>
tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> trvf
    (
        fvc::reconstructSmooth(tvf())
    );
    tvf.clear();
    return trvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
