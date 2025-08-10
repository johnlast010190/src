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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "regionSizeDistribution/regionSizeDistribution.H"
#include "regionSplit/regionSplit.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Map<Type> Foam::functionObjects::regionSizeDistribution::regionSum
(
    const regionSplit& regions,
    const Field<Type>& fld
) const
{
    // Per region the sum of fld
    Map<Type> regionToSum(regions.nRegions()/Pstream::nProcs());

    forAll(fld, celli)
    {
        label regioni = regions[celli];

        typename Map<Type>::iterator fnd = regionToSum.find(regioni);
        if (fnd == regionToSum.end())
        {
            regionToSum.insert(regioni, fld[celli]);
        }
        else
        {
            fnd() += fld[celli];
        }
    }
    Pstream::mapCombineGather(regionToSum, plusEqOp<Type>());
    Pstream::mapCombineScatter(regionToSum);

    return regionToSum;
}


template<class Type>
Foam::List<Type> Foam::functionObjects::regionSizeDistribution::extractData
(
    const UList<label>& keys,
    const Map<Type>& regionData
) const
{
    List<Type> sortedData(keys.size());

    forAll(keys, i)
    {
        sortedData[i] = regionData[keys[i]];
    }
    return sortedData;
}


// ************************************************************************* //
