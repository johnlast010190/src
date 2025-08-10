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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/ranges/scalarRange/scalarRanges.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarRanges::scalarRanges()
:
    List<scalarRange>(0)
{}


Foam::scalarRanges::scalarRanges(Istream& is)
:
    List<scalarRange>(0)
{
    DynamicList<scalarRange> lst;

    while (is.good())
    {
        scalarRange sr(is);
        if (sr.valid())
        {
            lst.append(sr);
        }
    }

    transfer(lst);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::scalarRanges::selected(const scalar value) const
{
    forAll(*this, i)
    {
        if (operator[](i).selected(value))
        {
            return true;
        }
    }

    return false;
}


Foam::List<bool> Foam::scalarRanges::selected
(
    const UList<scalar>& values
) const
{
    List<bool> lst(values.size(), false);

    // check ranges
    forAll(values, i)
    {
        if (selected(values[i]))
        {
            lst[i] = true;
        }
    }

    // check specific values
    forAll(*this, rangeI)
    {
        if (operator[](rangeI).isExact())
        {
            scalar target = operator[](rangeI).value();

            int nearestIndex = -1;
            scalar nearestDiff = Foam::GREAT;

            forAll(values, timeIndex)
            {
                scalar diff = fabs(values[timeIndex] - target);
                if (diff < nearestDiff)
                {
                    nearestDiff = diff;
                    nearestIndex = timeIndex;
                }
            }

            if (nearestIndex >= 0)
            {
                lst[nearestIndex] = true;
            }
        }
    }

    return lst;
}


Foam::List<Foam::scalar> Foam::scalarRanges::select
(
    const List<scalar>& values
) const
{
    return subset(selected(values), values);
}


void Foam::scalarRanges::inplaceSelect
(
    List<scalar>& values
) const
{
    inplaceSubset(selected(values), values);
}


// ************************************************************************* //
