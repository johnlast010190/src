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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::pTraits<Type>::labelType Foam::fvMesh::validComponents() const
{
    return pow
    (
        this->solutionD(),
        pTraits
        <
            typename powProduct<Vector<label>,
            pTraits<Type>::rank>::type
        >::zero
    );
}


template<class Type>
typename Foam::pTraits<Type>::labelType Foam::fvMesh::validComponents2() const
{
    return pow
    (
        (this->solutionD() + Vector<label>::one)/2,
        pTraits
        <
            typename powProduct<Vector<label>,
            pTraits<Type>::rank>::type
        >::zero
    );
}


template<class Type>
void Foam::fvMesh::stabiliseEmptyDirections(Field<Type>& t) const
{
    const typename pTraits<Type>::labelType v = validComponents2<Type>();
    const Type ident = pTraits<Type>::I;
    for (label cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        if (!v[cmpt])
        {
            t.replace(cmpt, ident[cmpt]);
        }
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fvMesh::stabiliseEmptyDirections
(
    GeometricField<Type, PatchField, GeoMesh>& t
) const
{
    const typename pTraits<Type>::labelType v = validComponents2<Type>();
    for (label cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        if (!v[cmpt])
        {
            const dimensioned<typename pTraits<Type>::cmptType> dimCmpt
            (
                "I", t.dimensions(), pTraits<Type>::I[cmpt]
            );
            t.replace(cmpt, dimCmpt);
        }
    }
}


// ************************************************************************* //
