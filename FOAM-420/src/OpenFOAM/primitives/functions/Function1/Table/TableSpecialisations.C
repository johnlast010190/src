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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Table/TableBase.H"


// * * * * * * * * * * * * Template specialisations  * * * * * * * * * * * * //

template<>
bool Foam::Function1Types::TableBase<bool>::value
(
    const scalar x
) const
{
    scalar xDash = x;

    if (checkMinBounds(x, xDash))
    {
        return table_.first().second();
    }

    if (checkMaxBounds(xDash, xDash))
    {
        return table_.last().second();
    }

    // Use interpolator
    interpolator().valueWeights(xDash, currentIndices_, currentWeights_);

    for (label i = 0; i < currentIndices_.size(); i++)
    {
        if (currentWeights_[i]>0)
        {
            return table_[currentIndices_[i]].second();
        }
    }

    WarningInFunction
        << "There are not positive weights: something is wrong"
        << endl;

    return false;
}

template<>
bool Foam::Function1Types::TableBase<bool>::integrate(const scalar x1, const scalar x2) const
{
    FatalErrorInFunction
        << "Table with bool type is not supported for integration"
        << exit(FatalError);
    return false;
}

template<>
bool Foam::Function1Types::TableBase<bool>::integrateYoverX(const scalar x1, const scalar x2) const
{
    FatalErrorInFunction
        << "Table with bool type is not supported for integration"
        << exit(FatalError);
    return false;
}


// ************************************************************************* //
