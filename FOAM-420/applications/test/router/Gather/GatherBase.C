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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "GatherBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type GatherBase::flatten(const List<Type> lst)
{
    label sum = 0;

    forAll(lst, lstI)
    {
        sum += lst[lstI].size();
    }

    Type result(sum);

    label index = 0;

    forAll(lst, lstI)
    {
        const Type& sub = lst[lstI];

        forAll(sub, subI)
        {
            result[index++] = sub[subI];
        }
    }

    return result;
}


template<class DataType, class IndexType, class AddOp>
IndexType GatherBase::offset
(
    const List<DataType>& values,
    const List<IndexType>& indices,
    AddOp aop
)
{
    if (values.size() != indices.size())
    {
        FatalErrorInFunction
            << "Input data and indices lists not equal size." << endl
            << "data size:" << values.size()
            << "  indices:" << indices.size()
            << abort(FatalError);
    }


    label sum = 0;

    forAll(indices, lstI)
    {
        sum += indices[lstI].size();
    }

    IndexType result(sum);

    label index = 0;

    label offset = 0;

    forAll(indices, lstI)
    {
        const IndexType& sub = indices[lstI];

        forAll(sub, subI)
        {
            result[index++] = aop(sub[subI], offset);
        }
        offset += values[lstI].size();
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
