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
    (c) Vuko Vukcevic, FMENA Zagreb.
    (c) Update by Hrvoje Jasak
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvBlockMatrix<Type>::negate()
{
    BlockLduSystem<Type, Type>::negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvBlockMatrix<Type>::operator=
(
    const fvBlockMatrix<Type>& bxs
)
{
    if (this == &bxs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    BlockLduSystem<Type, Type>::operator=(bxs);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator=
(
    const tmp<fvBlockMatrix<Type>>& tbxs
)
{
    operator=(tbxs());
    tbxs.clear();
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const fvBlockMatrix<Type>& bxs
)
{
    BlockLduSystem<Type, Type>::operator+=(bxs);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const tmp<fvBlockMatrix<Type>>& tbxs
)
{
    operator+=(tbxs());
    tbxs.clear();
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator+=
(
    fvBlockMatrix<Type2>& bxs
)
{
    bIndex_.check();
    direction sI = bIndex_.findStartIndex(bxs.psi().name());
    this->insertEquation(sI, bxs);
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator+=
(
    tmp<fvBlockMatrix<Type2>>& tbxs
)
{
    operator+=(tbxs.ref());
    tbxs.clear();
}



template<class Type>
void Foam::fvBlockMatrix<Type>::operator-=
(
    const fvBlockMatrix<Type>& bxs
)
{
    BlockLduSystem<Type, Type>::operator-=(bxs);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator-=
(
    tmp<fvBlockMatrix<Type>>& tbxs
)
{
    operator-=(tbxs());
    tbxs.clear();
}

template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator-=
(
    const fvBlockMatrix<Type2>& bxs
)
{
    bIndex_.check();
    direction sI = bIndex_.findStartIndex(bxs.psi().name());
    tmp<fvBlockMatrix<Type2>> tmbxs = -bxs;
    this->insertEquation(sI, tmbxs.ref());
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator-=
(
    tmp<fvBlockMatrix<Type2>>& tbxs
)
{
    operator-=(tbxs.ref());
    tbxs.clear();
}




template<class Type>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const fvMatrix<Type>& cbxs
)
{
    fvMatrix<Type>& bxs
    (
        const_cast<fvMatrix<Type>&>(cbxs)
    );
    insertEquation(0, bxs);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const tmp<fvMatrix<Type>>& tbxs
)
{
    operator+=(tbxs());
    tbxs.clear();
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const fvMatrix<Type2>& cbxs
)
{
    bIndex_.check();
    fvMatrix<Type2>& bxs
    (
        const_cast<fvMatrix<Type2>&>(cbxs)
    );
    direction sI = bIndex_.findStartIndex(bxs.psi().name());

    this->insertEquation(sI, bxs);
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const tmp<fvMatrix<Type2>>& tbxs
)
{
    operator+=(tbxs());
    tbxs.clear();
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator+=
(
    const DimensionedField<Type, volMesh>& su
)
{
//    checkMethod(*this, su, "+=");
    this->source() -= su.mesh().V()*su.field();
}




template<class Type>
void Foam::fvBlockMatrix<Type>::operator-=
(
    const fvMatrix<Type>& cbxs
)
{
    fvMatrix<Type>& bxs
    (
        const_cast<fvMatrix<Type>&>(cbxs)
    );
    insertEquation(0, -bxs);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator-=
(
    const tmp<fvMatrix<Type>>& tbxs
)
{
    operator-=(tbxs());
    tbxs.clear();
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator-=
(
    const fvMatrix<Type2>& cbxs
)
{
    bIndex_.check();
    fvMatrix<Type2>& bxs
    (
        const_cast<fvMatrix<Type2>&>(cbxs)
    );
    direction sI = bIndex_.findStartIndex(bxs.psi().name());

    this->insertEquation(sI, -bxs);
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator-=
(
    const tmp<fvMatrix<Type2>>& tbxs
)
{
    operator-=(tbxs());
    tbxs.clear();
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator*=
(
    const scalarField& sf
)
{
    BlockLduSystem<Type, Type>::operator*=(sf);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator*=
(
    const tmp<scalarField>& tsf
)
{
    operator*=(tsf);
    tsf.clear();
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator*=
(
    const scalar s
)
{
    BlockLduSystem<Type, Type>::operator*=(s);
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::operator()
(
    const word eqName,
    fvBlockMatrix<Type2>& bxs
)
{
    bIndex_.check();
    direction sI = bIndex_.findStartIndex(eqName);
    this->insertEquation(sI, bxs);
}


template<class Type>
template<class Type2, class Type3>
void Foam::fvBlockMatrix<Type>::operator()
(
    const word eqName,
    const word fieldName,
    BlockLduSystem<Type2, Type3>& bs
)
{
    bIndex_.check();
    direction iI = bIndex_.findIndex(eqName);
    direction iJ = bIndex_.findIndex(fieldName);
    direction sI = bIndex_.findStartIndex(eqName);
    direction sJ = bIndex_.findStartIndex(fieldName);
    bool increment = bIndex_.increment(iI, iJ);
    //Info<< sI << tab << sJ << tab << increment <<endl;

    this->insertBlockCoupling(sI, sJ, bs, increment);
}


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const fvBlockMatrix<Type>& A,
    const fvBlockMatrix<Type>& B,
    const char* op
)
{
    if (&A.psi() != &B.psi())
    {
        FatalErrorInFunction
            << "incompatible fields for operation "
            << endl << "    "
            << "[" << A.psi().name() << "] "
            << op
            << " [" << B.psi().name() << "]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvBlockMatrix<Type>& A,
    const fvMatrix<Type>& B,
    const char* op
)
{
    if (&A.psi() != &B.psi())
    {
        FatalErrorInFunction
            << "incompatible fields for operation "
            << endl << "    "
            << "[" << A.psi().name() << "] "
            << op
            << " [" << B.psi().name() << "]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const fvBlockMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const fvBlockMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += A.psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC().psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const fvBlockMatrix<Type>& A,
    const zero&
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator==
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const zero&
)
{
    return tA;
}



template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvBlockMatrix<Type>& A,
    const fvBlockMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const fvBlockMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvBlockMatrix<Type>& A,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvBlockMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvBlockMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tB(), A(), "+");
    tmp<fvBlockMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tB(), tA(), "+");
    tmp<fvBlockMatrix<Type>> tC(tB.ptr());
    tC.ref() += tA;
    tA.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvBlockMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvBlockMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const fvBlockMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A,
    const fvBlockMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const fvBlockMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvBlockMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tB(), A(), "-");
    tmp<fvBlockMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvBlockMatrix<Type>>& tB
)
{
    checkMethod(tB(), tA(), "-");
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB;
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const fvBlockMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<fvBlockMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvBlockMatrix<Type>& A
)
{
    tmp<fvBlockMatrix<Type>> tC(new fvBlockMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvBlockMatrix<Type>>& tA
)
{
    tmp<fvBlockMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
