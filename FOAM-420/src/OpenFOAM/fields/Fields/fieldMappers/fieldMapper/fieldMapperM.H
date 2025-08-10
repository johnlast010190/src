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
    (c) 2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "VectorN/Fields/VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DECLARE_FIELD_MAPPER_OPERATOR(type, Type, Modifier, args...)           \
                                                                               \
    virtual void operator()                                                    \
    (                                                                          \
        Field<type>& f,                                                        \
        const Field<type>& mapF                                                \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<type>> operator()                                        \
    (                                                                          \
        const Field<type>& mapF                                                \
    ) const Modifier;


#define DEFINE_FIELD_MAPPER_OPERATOR(type, Type, mapperType, args...)          \
                                                                               \
    void Foam::mapperType::operator()                                          \
    (                                                                          \
        Field<type>& f,                                                        \
        const Field<type>& mapF                                                \
    ) const                                                                    \
    {                                                                          \
        map(f, mapF);                                                          \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::type>> Foam::mapperType::operator()            \
    (                                                                          \
        const Field<type>& mapF                                                \
    ) const                                                                    \
    {                                                                          \
        return map(mapF);                                                      \
    }


#define DECLARE_DISTRIBUTED_MAPPER_OPERATOR(type, Type, Modifier, args...)     \
                                                                               \
    virtual void operator()                                                    \
    (                                                                          \
        Field<type>& f,                                                        \
        const Field<type>& mapF                                                \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<type>> operator()                                        \
    (                                                                          \
        const Field<type>& mapF                                                \
    ) const Modifier;                                                          \
                                                                               \
    virtual void operator()                                                    \
    (                                                                          \
        Field<type>& f,                                                        \
        const Field<type>& mapF,                                               \
        const bool applyFlip                                                   \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<type>> operator()                                        \
    (                                                                          \
        const Field<type>& mapF,                                               \
        const bool applyFlip                                                   \
    ) const Modifier;


#define DEFINE_DISTRIBUTED_MAPPER_OPERATOR(type, Type, mapperType, args...)    \
                                                                               \
    void Foam::mapperType::operator()                                          \
    (                                                                          \
        Field<type>& f,                                                        \
        const Field<type>& mapF                                                \
    ) const                                                                    \
    {                                                                          \
        map(f, mapF, true);                                                    \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::type>> Foam::mapperType::operator()            \
    (                                                                          \
        const Field<type>& mapF                                                \
    ) const                                                                    \
    {                                                                          \
        return map(mapF, true);                                                \
    }                                                                          \
                                                                               \
    void Foam::mapperType::operator()                                          \
    (                                                                          \
        Field<type>& f,                                                        \
        const Field<type>& mapF,                                               \
        const bool applyFlip                                                   \
    ) const                                                                    \
    {                                                                          \
        map(f, mapF, applyFlip);                                               \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::type>> Foam::mapperType::operator()            \
    (                                                                          \
        const Field<type>& mapF,                                               \
        const bool applyFlip                                                   \
    ) const                                                                    \
    {                                                                          \
        return map(mapF, applyFlip);                                           \
    }


// ************************************************************************* //
