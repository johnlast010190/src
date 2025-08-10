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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "sampledSetWriters/csv/csvSetWriter.H"
#include "coordSet/coordSet.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::csvSetWriter<Type>::csvSetWriter()
:
    writer<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::csvSetWriter<Type>::~csvSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const char* Foam::csvSetWriter<Type>::ext() const {
    return ".csv";
}

template<class Type>
void Foam::csvSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    writeHeader(points,valueSetNames,os);

    // Collect sets into columns
    List<const List<Type>*> columns(valueSets.size());

    forAll(valueSets, i)
    {
        columns[i] = valueSets[i];
    }

    this->writeTable(points, columns, os);
}


template<class Type>
void Foam::csvSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& points,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    writeHeader(points[0],valueSetNames,os);

    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorInFunction
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }

    List<const List<Type>*> columns(valueSets.size());

    forAll(points, trackI)
    {
        // Collect sets into columns
        forAll(valueSets, i)
        {
            columns[i] = &valueSets[i][trackI];
        }

        this->writeTable(points[trackI], columns, os);
        os  << nl << nl;
    }
}


template<class Type>
void Foam::csvSetWriter<Type>::writeSeparator(Ostream& os) const
{
    os << token::COMMA;
}


namespace Foam
{
    // otherwise compiler complains about specialization
    template<>
    void csvSetWriter<scalar>::writeHeader
    (
        const coordSet& points,
        const wordList& valueSetNames,
        Ostream& os
    ) const
    {
        writeCoordHeader(points, os);

        forAll(valueSetNames, i)
        {
            if (i > 0)
            {
                writeSeparator(os);
            }
            os << valueSetNames[i];
        }

        os << nl;
    }
}


template<class Type>
void Foam::csvSetWriter<Type>::writeHeader
(
    const coordSet& points,
    const wordList& valueSetNames,
    Ostream& os
) const
{
    writeCoordHeader(points, os);

    forAll(valueSetNames, i)
    {
        for (label j=0; j<Type::nComponents; j++)
        {
            if (i>0 || j>0)
            {
                writeSeparator(os);
            }
            os << valueSetNames[i] << "_" << j;
        }
    }

    os << nl;
}


template<class Type>
void Foam::csvSetWriter<Type>::writeCoordHeader
(
    const coordSet& points,
    Ostream& os
) const
{
    const word axisName(points.axis());

    if (points.hasVectorAxis())
    {
        for
        (
            word::const_iterator iter = axisName.begin();
            iter != axisName.end();
            ++iter
        )
        {
            os << *iter;
            writeSeparator(os);
        }
    }
    else
    {
        os << axisName;
        writeSeparator(os);
    }
}


// ************************************************************************* //
