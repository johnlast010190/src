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
    (c) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    (c) 2014-2018 Johan Heyns - CSIR, South Africa
    (c) 2021-2022 Esi Ltd.

Description
    Handy low-level functions to avoid writing face loops

\*---------------------------------------------------------------------------*/

#ifndef fvcCellFaceFunctions_H
#define fvcCellFaceFunctions_H

#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fields/surfaceFields/slicedSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{

//- Return the face-interpolate of the given cell field
//-  with the given owner and neighbour multipliers

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
ownerNeighbourSum
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tsfP,
    const tmp<surfaceScalarField>& tsfN
)
{
    const surfaceScalarField& sfP = tsfP();
    const surfaceScalarField& sfN = tsfN();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& sfPi = sfP.internalField();
    const scalarField& sfNi = sfN.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                vf.name()+'*'+sfP.name()+'+'+vf.name()+'*'+sfN.name(),
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()*sfP.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    if (sfP.dimensions() != sfN.dimensions())
    {
        FatalErrorInFunction << nl
                     << "Dimensions of sfP and sfN differ: " << nl
                     << sfP.dimensions() << " != " << sfN.dimensions() << nl
                     << exit(FatalError);
    }

    Field<Type>& sfi = sf.internalField();

    forAll(P, fi)
    {
        sfi[fi] = sfPi[fi]*vfi[P[fi]] + sfNi[fi]*vfi[N[fi]];
    }

    forAll(sfP.boundaryField(), pi)
    {
        const fvsPatchScalarField& psfP = sfP.boundaryField()[pi];
        const fvsPatchScalarField& psfN = sfN.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            sf.boundaryField()[pi] =
                psfP*vf.boundaryField()[pi].patchInternalField()
              + psfN*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = (psfP+psfN)*vf.boundaryField()[pi];
        }
    }

    tsfP.clear();
    tsfN.clear();

    return tsf;
}


template<class Type>
tmp<Field<Type>> ownerField
(
    const Field<Type>& vfi, const polyMesh& mesh
)
{
    const labelUList& P = mesh.faceOwner();
    tmp<Field<Type>> tsf(new Field<Type>(P.size()));
    Field<Type>& sf = tsf.ref();

    forAll(P, fi)
    {
        sf[fi] = vfi[P[fi]];
    }

    return tsf;
}


//- Return the owner cell values at each face
//  Boundaries contain internal cell values

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
ownerField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "owner("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    forAll(P, fi)
    {
        sf[fi] = vfi[P[fi]];
    }

    forAll(vf.boundaryField(), pi)
    {
        sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi].patchInternalField();
    }

    return tsf;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
ownerField
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        ownerField(tvf());
    tvf.clear();
    return tsf;
}


template<class Type>
tmp<Field<Type>> neighbourField
(
    const Field<Type>& vfi, const polyMesh& mesh
)
{
    const labelUList& N = mesh.faceNeighbour();
    tmp<Field<Type>> tsf(new Field<Type>(N.size()));
    Field<Type>& sf = tsf.ref();

    forAll(N, fi)
    {
        sf[fi] = vfi[N[fi]];
    }

    return tsf;
}


//- Return the neighbour cell values at each face
//  Boundaries contain boundary face values

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
neighbourField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "neighbour("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    forAll(N, fi)
    {
        sfi[fi] = vfi[N[fi]];
    }

    forAll(vf.boundaryField(), pi)
    {
        if (vf.boundaryField()[pi].coupled())
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    return tsf;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
neighbourField
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        neighbourField(tvf());
    tvf.clear();
    return tsf;
}


//- Return the surface sum taking into account different owner and
//- neighbour fields

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfO,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfN
)
{
    const fvMesh& mesh = sfO.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+sfO.name()+","+sfN.name()+')',
                sfO.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", sfO.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    if (sfO.dimensions() != sfN.dimensions())
    {
        FatalErrorInFunction << nl
                     << "Dimensions of sfO and sfN differ: " << nl
                     << sfO.dimensions() << " != " << sfN.dimensions() << nl
                     << exit(FatalError);
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += sfO[facei];
        vf[neighbour[facei]] += sfN[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& psfO = sfO.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += psfO[facei];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceSum
(
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfO,
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfN
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf =
        surfaceSum(tsfO(), tsfN());
    tsfO.clear();
    tsfN.clear();
    return tvf;
}


//- Return the surface sum taking into account different owner and
//- neighbour fields, doing only coupled boundary faces and ignoring
//- physical boundaries

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
intSurfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfO,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfN
)
{
    const fvMesh& mesh = sfO.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+sfO.name()+","+sfN.name()+')',
                sfO.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", sfO.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    if (sfO.dimensions() != sfN.dimensions())
    {
        FatalErrorInFunction << nl
                     << "Dimensions of sfO and sfN differ: " << nl
                     << sfO.dimensions() << " != " << sfN.dimensions() << nl
                     << exit(FatalError);
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += sfO[facei];
        vf[neighbour[facei]] += sfN[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& psfO = sfO.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                vf[pFaceCells[facei]] += psfO[facei];
            }
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
intSurfaceSum
(
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfO,
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfN
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf = intSurfaceSum(tsfO(), tsfN());
    tsfO.clear();
    tsfN.clear();
    return tvf;
}


//- Calculate surface sum including only coupled boundaries

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
intSurfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = sf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+sf.name()+')',
                sf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", sf.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += sf[facei];
        vf[neighbour[facei]] += sf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& psf = sf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                vf[pFaceCells[facei]] += psf[facei];
            }
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


// Select from one or the other surface field comparing to a threshold value
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
surfaceFieldSelect
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfL,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfR,
    const surfaceScalarField& select,
    const scalar thres
)
{
    const fvMesh& mesh = sfL.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "select("+sfL.name()+","+sfR.name()+')',
                sfL.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            sfL.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    if (sfL.dimensions() != sfR.dimensions())
    {
        FatalErrorInFunction << "Dimensions of sfL and sfR differ: " << nl
                     << sfL.dimensions() << " != " << sfL.dimensions() << nl
                     << exit(FatalError);
    }

    forAll(select, facei)
    {
        sf[facei] = (select[facei] >= thres ? sfL[facei] : sfR[facei]);
    }

    forAll(select.boundaryField(), patchi)
    {
        Field<Type>& psf = sf.boundaryFieldRef()[patchi];
        const Field<Type>& psfL = sfL.boundaryField()[patchi];
        const Field<Type>& psfR = sfR.boundaryField()[patchi];
        const scalarField& pSel = select.boundaryField()[patchi];
        forAll(pSel, facei)
        {
            psf[facei] = (pSel[facei] >= thres ? psfL[facei] : psfR[facei]);
        }
    }

    return tsf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
surfaceFieldSelect
(
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& sfL,
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& sfR,
    const surfaceScalarField& select,
    const scalar thres
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        surfaceFieldSelect(sfL(), sfR(), select, thres);
    sfL.clear();
    sfR.clear();
    return tsf;
}


// Call a function for each face of a field; return result as a field
template<class Type, class Lambda>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
calcForEachFace
(
    Lambda fn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = sf.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfout
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "calc("+sf.name()+')',
                sf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            sf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfout = tsfout.ref();

    forAll(sfout, facei)
    {
        sfout[facei] = fn(sf[facei]);
    }

    forAll(sfout.boundaryField(), patchi)
    {
        Field<Type>& psfout = sfout.boundaryFieldRef()[patchi];
        const Field<Type>& psf = sf.boundaryField()[patchi];
        forAll(psfout, facei)
        {
            psfout[facei] = fn(psf[facei]);
        }
    }

    return tsfout;
}


// Call a function for each face of a field; return result as a field
template<class Type, class Lambda>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
calcForEachFace
(
    Lambda fn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf1,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf2
)
{
    const fvMesh& mesh = sf1.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfout
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "calc("+sf1.name()+','+sf2.name()+')',
                sf1.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            sf1.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfout = tsfout.ref();

    forAll(sfout, facei)
    {
        sfout[facei] = fn(sf1[facei], sf2[facei]);
    }

    forAll(sfout.boundaryField(), patchi)
    {
        Field<Type>& psfout = sfout.boundaryFieldRef()[patchi];
        const Field<Type>& psf1 = sf1.boundaryField()[patchi];
        const Field<Type>& psf2 = sf2.boundaryField()[patchi];
        forAll(psfout, facei)
        {
            psfout[facei] = fn(psf1[facei], psf2[facei]);
        }
    }

    return tsfout;
}


// Sum the co-incident internal and GIB face contributions and sync both to the
// owner and neighbour fields
template<class Type>
void sumIndirectAndInternalFaces
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfOwn,
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfNei
)
{
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bsfOwn =
        sfOwn.boundaryFieldRef();
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bsfNei =
        sfNei.boundaryFieldRef();

    const fvMesh& mesh = sfOwn.mesh();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvP = mesh.boundary()[patchi];

        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            const indirectPolyPatch& ipp =
                refCast<const indirectPolyPatch>(fvP.patch());
            const labelList& addr = ipp.addressing();
            const boolList& fm = ipp.fm();

            fvsPatchField<Type>& psfOwn = bsfOwn[patchi];
            fvsPatchField<Type>& psfNei = bsfNei[patchi];

            if (fvP.isActive())
            {
                bool master = (ipp.indirectPolyPatchType() == "master");
                forAll(addr, pfacei)
                {
                    const label& fi = addr[pfacei];
                    if (mesh.isInternalFace(fi))
                    {
                        if (master != fm[pfacei])
                        {
                            Type sum = sfOwn[fi] + psfOwn[pfacei];
                            sfOwn[fi] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                        else
                        {
                            Type sum = sfNei[fi] + psfOwn[pfacei];
                            sfNei[fi] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                    }
                    else
                    {
                        label pI = mesh.boundaryMesh().whichPatch(fi);
                        label lpfI = mesh.boundaryMesh()[pI].whichFace(fi);
                        if (master != fm[fi])
                        {
                            Type sum = bsfOwn[pI][lpfI] + psfOwn[pfacei];
                            bsfOwn[pI][lpfI] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                        if (master != fm[fi])
                        {
                            Type sum = bsfNei[pI][lpfI] + psfOwn[pfacei];
                            bsfNei[pI][lpfI] = sum;
                            psfOwn[pfacei] = sum;
                            psfNei[pfacei] = sum;
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void applyFaceMaskTo
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = sf.mesh();
    Field<Type>& isf = sf.primitiveFieldRef();
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bsf =
        sf.boundaryFieldRef();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvP = mesh.boundary()[patchi];

        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            // This blending factor is the proportion of the internal (non indirect-
            // patch) face value that is included for a 'porous' GIB boundary.
            // Currently a placeholder - still to be made settable and generalised
            // to a per-patch/per-face basis
            scalar blendingFactor = 0;

            const indirectPolyPatch& ipp =
                refCast<const indirectPolyPatch>(fvP.patch());
            const labelList& addr = ipp.addressing();

            if (fvP.isActive())
            {
                bool master = (ipp.indirectPolyPatchType() == "master");
                forAll(addr, pfacei)
                {
                    const label& fi = addr[pfacei];
                    if (mesh.isInternalFace(fi))
                    {
                        if (master)
                        {
                            isf[fi] *= blendingFactor;
                        }
                    }
                    else
                    {
                        label pI = mesh.boundaryMesh().whichPatch(fi);
                        label lpfI = mesh.boundaryMesh()[pI].whichFace(fi);
                        bsf[pI][lpfI] *= blendingFactor;
                    }
                }

                bsf[patchi] *= scalar(1)-blendingFactor;
            }
            else
            {
                bsf[patchi] *= Zero;
            }
        }
    }
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
applyFaceMask
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = sf.mesh();
    bool haveIndirectPatches = false;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvP = mesh.boundary()[patchi];

        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            haveIndirectPatches = true;
            break;
        }
    }

    if (!haveIndirectPatches)
    {
        // Return a const ref; no copy
        return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>(sf);
    }

    typedef
        SlicedGeometricField
        <
            Type,
            fvsPatchField,
            slicedFvsPatchField,
            surfaceMesh
        >
        slicedFieldType;
    if (isA<slicedFieldType>(sf))
    {
        // For a sliced field, we need to force a deep copy, otherwise the
        // original referenced field will be altered
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tmsf
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                sf,
                sf.mesh(),
                sf.dimensions()
            )
        );
        tmsf.ref() = sf;
        applyFaceMaskTo(tmsf.ref());
        return tmsf;
    }
    else
    {
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tmsf
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>(sf)
        );
        applyFaceMaskTo(tmsf.ref());
        return tmsf;
    }
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
applyFaceMask
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tsf
)
{
    applyFaceMaskTo(tsf.ref());
    return tsf;
}


}  // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

