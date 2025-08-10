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
    (c) 2016 OpenCFD Ltd.
    (c) 2010-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#include "fields/fvPatchFields/basic/coupled/coupledFvPatchFields.H"
#include "containers/Lists/UIndirectList/UIndirectList.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "containers/Lists/CompactListList/CompactListList.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvm/fvmDdt.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<template<class> class ListType>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const ListType<Type>& values
)
{
    // Mark the GIB faces
    const boolList gibFaces = markGIBFaces();

    // Fix the values
    forAll(cellLabels, i)
    {
        setValue(cellLabels[i], values[i], gibFaces);
    }
}


template<class Type>
template<template<class> class ListType>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const ListType<Type>& values,
    const scalarList& fractions,
    const bool hasDdt
)
{
    scalarField ddtDiag(diag().size(), 0);

    bool hasFrac = false;
    forAll(fractions, i)
    {
        if (fractions[i] < 1 - ROOTSMALL)
        {
            hasFrac = true;
            break;
        }
    }

    if (hasFrac)
    {
        // Get the proportion of the diagonal associated with iterative update
        const scalar alpha = relaxationFactor();
        if (alpha > 0)
        {
            ddtDiag += (1 - alpha)*diag();
        }
        if (hasDdt)
        {
            const fvMatrix<Type> ddtEqn(fvm::ddt(psi_));
            if (ddtEqn.hasDiag())
            {
                ddtDiag += ddtEqn.diag();
            }
        }
    }

    // Mark the GIB faces
    const boolList gibFaces = markGIBFaces();

    forAll(cellLabels, i)
    {
        if (- ROOTSMALL < fractions[i] && fractions[i] < ROOTSMALL)
        {
            // Do nothing
        }
        else if (1 - ROOTSMALL < fractions[i] && fractions[i] < 1 + ROOTSMALL)
        {
            // Fix the values
            setValue(cellLabels[i], values[i], gibFaces);
        }
        else
        {
            // Fractionally fix the values
            setValue(cellLabels[i], values[i], fractions[i], ddtDiag);
        }
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm1,
    const fvMatrix<Type>& fvm2,
    const char* op
)
{
    if (&fvm1.psi() != &fvm2.psi())
    {
        FatalErrorInFunction
            << "incompatible fields for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << "] "
            << op
            << " [" << fvm2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fvm1.dimensions() != fvm2.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << fvm1.dimensions()/dimVolume << " ] "
            << op
            << " [" << fvm2.psi().name() << fvm2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm,
    const DimensionedField<Type, volMesh>& df,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != df.dimensions())
    {
        FatalErrorInFunction
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << df.name() << df.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != dt.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    fvMatrix<Type>& fvm,
    const dictionary& solverControls
)
{
    return fvm.solve(solverControls);
}

template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    const tmp<fvMatrix<Type>>& tfvm,
    const dictionary& solverControls
)
{
    SolverPerformance<Type> solverPerf =
        const_cast<fvMatrix<Type>&>(tfvm()).solve(solverControls);

    tfvm.clear();

    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve(fvMatrix<Type>& fvm)
{
    return fvm.solve();
}

template<class Type>
Foam::SolverPerformance<Type> Foam::solve(const tmp<fvMatrix<Type>>& tfvm)
{
    SolverPerformance<Type> solverPerf =
        const_cast<fvMatrix<Type>&>(tfvm()).solve();

    tfvm.clear();

    return solverPerf;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::evaluate
(
    const fvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& eqnVar
)
{
    tmp<Foam::fvMatrix<Type>> tAeval
    (
        new fvMatrix<Type>(eqnVar, A.dimensions())
    );
    tAeval.ref() += (A & A.psi());

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().schemes().fluxRequired(A.psi().name())
    )
    {
        tAeval.ref().faceFluxCorrectionPtr() = (A.flux()).ptr();
    }

    return tAeval;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::evaluate
(
    const tmp<fvMatrix<Type>>& tA,
    const GeometricField<Type, fvPatchField, volMesh>& eqnVar
)
{
    tmp<Foam::fvMatrix<Type>> tAeval = evaluate(tA(), eqnVar);
    tA.clear();
    return tAeval;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::correction
(
    const fvMatrix<Type>& A
)
{
    tmp<Foam::fvMatrix<Type>> tAcorr = A - (A & A.psi());

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().schemes().fluxRequired(A.psi().name())
    )
    {
        tAcorr.ref().faceFluxCorrectionPtr() = (-A.flux()).ptr();
        if (A.faceFluxCorrectionPtr())
        {
            *tAcorr().faceFluxCorrectionPtr() += *A.faceFluxCorrectionPtr();
        }
    }

    return tAcorr;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::correction
(
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<Foam::fvMatrix<Type>> tAcorr(tA.ptr());
    tAcorr.ref() -= (tAcorr() & tAcorr().psi());

    // Note the matrix coefficients are still that of matrix A
    const fvMatrix<Type>& A = tAcorr();

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().schemes().fluxRequired(A.psi().name())
    )
    {
        if (A.faceFluxCorrectionPtr())
        {
            *tAcorr().faceFluxCorrectionPtr() -= A.flux();
        }
        else
        {
            tAcorr.ref().faceFluxCorrectionPtr() = (-A.flux()).ptr();
        }
    }

    return tAcorr;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::changeVariable
(
    const tmp<fvMatrix<Type>>& tA,
    const GeometricField<scalar, fvPatchField, volMesh>& dOlddNew,
    const GeometricField<Type, fvPatchField, volMesh>& newVar
)
{
    const fvMesh& mesh = newVar.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxCorr =
        tA().flux();

    const GeometricField<Type, fvPatchField, volMesh>& oldVar =
        tA().psi();

    tmp<Foam::fvMatrix<Type>> tALin
    (
        new fvMatrix<Type>(newVar, tA)
    );
    fvMatrix<Type>& ALin = tALin.ref();

    if (ALin.hasLower() || ALin.hasUpper())
    {
        const labelUList& owner = mesh.owner();
        scalarField& lower = ALin.lower();
        forAll(owner, fi)
        {
            lower[fi] *= dOlddNew[owner[fi]];
        }

        const labelUList& neighbour = mesh.neighbour();
        scalarField& upper = ALin.upper();
        forAll(neighbour, fi)
        {
            upper[fi] *= dOlddNew[neighbour[fi]];
        }
    }
    if (ALin.hasDiag())
    {
        ALin.diag() *= dOlddNew;
    }

    FieldField<Field,Type>& bc = ALin.boundaryCoeffs();
    FieldField<Field,Type>& ic = ALin.internalCoeffs();
    const GeometricField<scalar, fvPatchField, volMesh>::Boundary& bdOlddNew =
        dOlddNew.boundaryField();
    forAll(bc, pi)
    {
        ic[pi] *= bdOlddNew[pi].patchInternalField();
        if (newVar.boundaryField()[pi].coupled())
        {
            bc[pi] *= bdOlddNew[pi].patchNeighbourField();
        }
        else if (newVar.boundaryField()[pi].regionCoupled())
        {
            bc[pi] *= bdOlddNew[pi];
        }
    }

    fluxCorr.ref() -= ALin.flux();
    DimensionedField<Type, volMesh> intFluxCorr
    (
        IOobject
        (
            "intFluxCorr",
            fluxCorr().instance(),
            fluxCorr().db()
        ),
        mesh,
        dimensioned<Type>("0", fluxCorr().dimensions()/dimVol, Zero)
    );
    // Using this form prevents an unnecessary boundary update:
    fvc::surfaceIntegrate(intFluxCorr, fluxCorr());
    ALin += intFluxCorr;

    if (ALin.faceFluxCorrectionPtr())
    {
        *ALin.faceFluxCorrectionPtr() += fluxCorr();
    }
    else if
    (
        (ALin.hasUpper() || ALin.hasLower())
     && (
            mesh.schemes().fluxRequired(oldVar.name())
         || mesh.schemes().fluxRequired(newVar.name())
        )
    )
    {
        ALin.faceFluxCorrectionPtr() = fluxCorr.ptr();
    }

    return tALin;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += A.psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC().psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const zero&
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const zero&
)
{
    return tA;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const volScalarField::Internal& dsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField::Internal>& tdsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const volScalarField::Internal& dsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField::Internal>& tdsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type>> tC(new fvMatrix<Type>(A));
    tC.ref() *= ds;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "M&" + psi.name(),
                psi.instance(),
                psi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi.mesh(),
            M.dimensions()/dimVol,
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Mphi = tMphi.ref();

    // Loop over field components
    if (M.hasDiag())
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            scalarField psiCmpt(psi.field().component(cmpt));
            scalarField boundaryDiagCmpt(M.diag());
            M.addBoundaryDiag(boundaryDiagCmpt, cmpt);
            Mphi.primitiveFieldRef().replace(cmpt, -boundaryDiagCmpt*psiCmpt);
        }
    }
    else
    {
        Mphi.primitiveFieldRef() = Zero;
    }

    Mphi.primitiveFieldRef() += M.lduMatrix::H(psi.field()) + M.source();
    M.addBoundarySource(Mphi.primitiveFieldRef(), M.psi());

    Mphi.primitiveFieldRef() /= -psi.mesh().V();
    Mphi.correctBoundaryConditions();

    return tMphi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<DimensionedField<Type, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const tmp<DimensionedField<Type, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvMatrix<Type>& fvm)
{
    os  << static_cast<const lduMatrix&>(fvm) << nl
        << fvm.dimensions_ << nl
        << fvm.source_ << nl
        << fvm.internalCoeffs_ << nl
        << fvm.boundaryCoeffs_ << endl;

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
