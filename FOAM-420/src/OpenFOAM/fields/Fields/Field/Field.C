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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/Field/Field.H"
#include "fields/Fields/Field/FieldM.H"
#include "db/dictionary/dictionary.H"
#include "primitives/contiguous/contiguous.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributeBase.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class Type>
const char* const Foam::Field<Type>::typeName("Field");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::UList<Type>& Foam::Field<Type>::copySelf
(
    const UList<Type>& mapF,
    tmp<Field<Type>>& tmapF
) const
{
    if (static_cast<const UList<Type>*>(this) == &mapF)
    {
        tmapF = clone();
    }

    return tmapF.valid() ? tmapF() : mapF;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>::Field()
:
    List<Type>()
{}


template<class Type>
Foam::Field<Type>::Field(const label size)
:
    List<Type>(size)
{}


template<class Type>
Foam::Field<Type>::Field(const label size, const Type& t)
:
    List<Type>(size, t)
{}


template<class Type>
Foam::Field<Type>::Field(const label size, const zero)
:
    List<Type>(size, Zero)
{}


template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const labelUList& mapAddressing
)
:
    List<Type>(mapAddressing.size())
{
    map(mapF, mapAddressing);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing
)
:
    List<Type>(mapAddressing.size())
{
    map(tmapF, mapAddressing);
}


template<class Type>
Foam::Field<Type>::Field
(
    const UList<Type>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
:
    List<Type>(mapAddressing.size())
{
    map(mapF, mapAddressing, mapWeights);
}


template<class Type>
Foam::Field<Type>::Field
(
    const tmp<Field<Type>>& tmapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
:
    List<Type>(mapAddressing.size())
{
    map(tmapF, mapAddressing, mapWeights);
}


template<class Type>
Foam::Field<Type>::Field(const Field<Type>& f)
:
    tmp<Field<Type>>::refCount(),
    List<Type>(f)
{}


template<class Type>
Foam::Field<Type>::Field(Field<Type>&& f)
:
    tmp<Field<Type>>::refCount(),
    List<Type>(std::move(f))
{}


template<class Type>
Foam::Field<Type>::Field(Field<Type>& f, bool reuse)
:
    List<Type>(f, reuse)
{}


template<class Type>
Foam::Field<Type>::Field(const Xfer<List<Type>>& f)
:
    List<Type>(f)
{}


template<class Type>
Foam::Field<Type>::Field(const Xfer<Field<Type>>& f)
:
    List<Type>()
{
    List<Type>::transfer(f());
}


template<class Type>
Foam::Field<Type>::Field(const UList<Type>& list)
:
    List<Type>(list)
{}


template<class Type>
Foam::Field<Type>::Field(const UIndirectList<Type>& list)
:
    List<Type>(list)
{}


#ifndef NoConstructFromTmp
template<class Type>
Foam::Field<Type>::Field(const tmp<Field<Type>>& tf)
:
    List<Type>(const_cast<Field<Type>&>(tf()), tf.isTmp())
{
    tf.clear();
}
#endif


template<class Type>
Foam::Field<Type>::Field(Istream& is)
:
    List<Type>(is)
{}


template<class Type>
Foam::Field<Type>::Field
(
    const word& keyword,
    const dictionary& dict,
    const label s
)
{
    if (s)
    {
        ITstream& is = dict.lookup(keyword);

        // Read first token
        token firstToken(is);

        if (firstToken.isWord())
        {
            if (firstToken.wordToken() == "uniform")
            {
                this->setSize(s);
                operator=(pTraits<Type>(is));
            }
            else if (firstToken.wordToken() == "nonuniform")
            {
                is >> static_cast<List<Type>&>(*this);
                if (this->size() != s)
                {
                    FatalIOErrorInFunction
                    (
                        dict
                    )   << "size " << this->size()
                        << " is not equal to the given value of " << s
                        << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', found "
                    << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "expected keyword 'uniform' or 'nonuniform', found "
                << firstToken.info()
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Field<Type>::clone() const
{
    return tmp<Field<Type>>(new Field<Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Field<Type>::map
(
    const UList<Type>& mapF0,
    const labelUList& mapAddressing
)
{
    Field<Type>& f = *this;

    tmp<Field<Type>> tmapF;
    const UList<Type>& mapF = copySelf(mapF0, tmapF);

    if (f.size() != mapAddressing.size())
    {
        f.setSize(mapAddressing.size(), Type(Zero));
    }

    if (mapF.size() > 0)
    {
        forAll(f, i)
        {
            const label mapI = mapAddressing[i];

            if (mapI >= 0)
            {
                f[i] = mapF[mapI];
            }
            else
            {
                f[i] = Type(Zero);
            }
        }
    }
    else
    {
        f = Type(Zero);
    }
}


template<class Type>
void Foam::Field<Type>::map
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing
)
{
    map(tmapF(), mapAddressing);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::map
(
    const UList<Type>& mapF0,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
{
    if (mapWeights.size() != mapAddressing.size())
    {
        FatalErrorInFunction
            << mapWeights.size() << " map size: " << mapAddressing.size()
            << abort(FatalError);
    }

    Field<Type>& f = *this;

    tmp<Field<Type>> tmapF;
    const UList<Type>& mapF = copySelf(mapF0, tmapF);

    if (this->size() != mapAddressing.size())
    {
        this->setSize(mapAddressing.size());
    }

    forAll(f, i)
    {
        const labelList& localAddrs = mapAddressing[i];
        const scalarList& localWeights = mapWeights[i];

        f[i] = Zero;

        forAll(localAddrs, j)
        {
            f[i] += localWeights[j]*mapF[localAddrs[j]];
        }
    }
}


template<class Type>
void Foam::Field<Type>::map
(
    const tmp<Field<Type>>& tmapF,
    const labelListList& mapAddressing,
    const scalarListList& mapWeights
)
{
    map(tmapF(), mapAddressing, mapWeights);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const UList<Type>& mapF0,
    const labelUList& mapAddressing
)
{
    Field<Type>& f = *this;

    tmp<Field<Type>> tmapF;
    const UList<Type>& mapF = copySelf(mapF0, tmapF);

    forAll(mapF, i)
    {
        const label mapI = mapAddressing[i];

        if (mapI >= 0)
        {
            f[mapI] = mapF[i];
        }
    }
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing
)
{
    rmap(tmapF(), mapAddressing);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const UList<Type>& mapF0,
    const labelUList& mapAddressing,
    const UList<scalar>& mapWeights
)
{
    Field<Type>& f = *this;

    tmp<Field<Type>> tmapF;
    const UList<Type>& mapF = copySelf(mapF0, tmapF);

    f = Zero;

    forAll(mapF, i)
    {
        f[mapAddressing[i]] += mapF[i]*mapWeights[i];
    }
}


template<class Type>
void Foam::Field<Type>::rmap
(
    const tmp<Field<Type>>& tmapF,
    const labelUList& mapAddressing,
    const UList<scalar>& mapWeights
)
{
    rmap(tmapF(), mapAddressing, mapWeights);
    tmapF.clear();
}


template<class Type>
void Foam::Field<Type>::negate()
{
    TFOR_ALL_F_OP_OP_F(Type, *this, =, -, Type, *this)
}


template<class Type>
Foam::tmp<Foam::Field<typename Foam::Field<Type>::cmptType>>
Foam::Field<Type>::component
(
    const direction d
) const
{
    tmp<Field<cmptType>> Component(new Field<cmptType>(this->size()));
    ::Foam::component(Component.ref(), *this, d);
    return Component;
}


template<class Type>
void Foam::Field<Type>::replace
(
    const direction d,
    const UList<cmptType>& sf
)
{
    TFOR_ALL_F_OP_FUNC_S_F(Type, *this, ., replace, const direction, d,
        cmptType, sf)
}


template<class Type>
void Foam::Field<Type>::replace
(
    const direction d,
    const tmp<Field<cmptType>>& tsf
)
{
    replace(d, tsf());
    tsf.clear();
}


template<class Type>
void Foam::Field<Type>::replace
(
    const direction d,
    const cmptType& c
)
{
    TFOR_ALL_F_OP_FUNC_S_S(Type, *this, ., replace, const direction, d,
        cmptType, c)
}


template<class Type>
template<class VSForm>
VSForm Foam::Field<Type>::block(const label start) const
{
    VSForm vs;
    for (direction i=0; i<VSForm::nComponents; i++)
    {
        vs[i] = this->operator[](start + i);
    }
    return vs;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Field<Type>::T() const
{
    tmp<Field<Type>> transpose(new Field<Type>(this->size()));
    ::Foam::T(transpose.ref(), *this);
    return transpose;
}


template<class Type>
void Foam::Field<Type>::writeEntry
(
    const word& keyword,
    Ostream& os,
    bool flushStream
) const
{
    os.writeKeyword(keyword);

    // Can the contents be considered 'uniform' (ie, identical)?
    bool uniform = (this->size() && contiguous<Type>());
    if (uniform)
    {
        forAll(*this, i)
        {
            if (this->operator[](i) != this->operator[](0))
            {
                uniform = false;
                break;
            }
        }
    }

    if (uniform)
    {
        os << "uniform " << this->operator[](0) << token::END_STATEMENT;
    }
    else
    {
        os << "nonuniform ";
        List<Type>::writeEntry(os);
        os << token::END_STATEMENT;
    }

    if (flushStream)
    {
        os << endl;
    }
    else
    {
        os << nl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Field<Type>::operator=(const Field<Type>& rhs)
{
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    List<Type>::operator=(rhs);
}


template<class Type>
void Foam::Field<Type>::operator=(Field<Type>&& rhs)
{
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    List<Type>::operator=(std::move(rhs));
}


template<class Type>
void Foam::Field<Type>::operator=(const SubField<Type>& rhs)
{
    List<Type>::operator=(rhs);
}


template<class Type>
void Foam::Field<Type>::operator=(const UList<Type>& rhs)
{
    List<Type>::operator=(rhs);
}


template<class Type>
void Foam::Field<Type>::operator=(List<Type>&& rhs)
{
    List<Type>::operator=(std::move(rhs));
}


template<class Type>
void Foam::Field<Type>::operator=(const tmp<Field>& rhs)
{
    if (this == &(rhs()))
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (rhs.isTmp())
    {
        List<Type>::transfer(rhs.ref());
    }
    else
    {
        List<Type>::operator=(rhs());
    }
    rhs.clear();
}


template<class Type>
void Foam::Field<Type>::operator=(const Type& t)
{
    List<Type>::operator=(t);
}


template<class Type>
void Foam::Field<Type>::operator=(const zero)
{
    List<Type>::operator=(Zero);
}


template<class Type>
template<class Form, class Cmpt, Foam::direction nCmpt>
void Foam::Field<Type>::operator=(const VectorSpace<Form,Cmpt,nCmpt>& vs)
{
    TFOR_ALL_F_OP_S(Type, *this, =, VSType, vs)
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type>                                                           \
void Foam::Field<Type>::operator op(const UList<TYPE>& f)                      \
{                                                                              \
    TFOR_ALL_F_OP_F(Type, *this, op, TYPE, f)                                  \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::Field<Type>::operator op(const tmp<Field<TYPE>>& tf)                \
{                                                                              \
    operator op(tf());                                                         \
    tf.clear();                                                                \
}                                                                              \
                                                                               \
template<class Type>                                                           \
void Foam::Field<Type>::operator op(const TYPE& t)                             \
{                                                                              \
    TFOR_ALL_F_OP_S(Type, *this, op, TYPE, t)                                  \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const Field<Type>& f)
{
    os  << static_cast<const List<Type>&>(f);
    return os;
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const tmp<Field<Type>>& tf)
{
    os  << tf();
    tf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/Fields/Field/FieldFunctions.C"

// ************************************************************************* //
