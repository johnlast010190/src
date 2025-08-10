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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
inline const Foam::DimensionedField<Type, GeoMesh>&
Foam::DimensionedField<Type, GeoMesh>::null()
{
    return NullObjectRef<DimensionedField<Type, GeoMesh>>();
}


template<class Type, class GeoMesh>
inline const typename GeoMesh::Mesh&
Foam::DimensionedField<Type, GeoMesh>::mesh() const
{
    return mesh_;
}


template<class Type, class GeoMesh>
inline const Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh>::dimensions() const
{
    return dimensions_;
}


template<class Type, class GeoMesh>
inline Foam::dimensionSet& Foam::DimensionedField<Type, GeoMesh>::dimensions()
{
    return dimensions_;
}


template<class Type, class GeoMesh>
inline const Foam::orientedType&
Foam::DimensionedField<Type, GeoMesh>::oriented() const
{
    return oriented_;
}


template<class Type, class GeoMesh>
inline Foam::orientedType& Foam::DimensionedField<Type, GeoMesh>::oriented()
{
    return oriented_;
}


template<class Type, class GeoMesh>
inline void Foam::DimensionedField<Type, GeoMesh>::setOriented
(
    const bool oriented
)
{
    oriented_.setOriented(oriented);
}


template<class Type, class GeoMesh>
inline const Foam::Field<Type>&
Foam::DimensionedField<Type, GeoMesh>::field() const
{
    return *this;
}


template<class Type, class GeoMesh>
inline Foam::Field<Type>&
Foam::DimensionedField<Type, GeoMesh>::field()
{
    return *this;
}


// ************************************************************************* //
