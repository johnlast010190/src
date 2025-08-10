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

#include "sampledSetWriters/gnuplot/gnuplotSetWriter.H"
#include "coordSet/coordSet.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::gnuplotSetWriter<Type>::gnuplotSetWriter()
:
    writer<Type>()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::gnuplotSetWriter<Type>::~gnuplotSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const char* Foam::gnuplotSetWriter<Type>::ext() const {
    return ".gplt";
}

template<class Type>
void Foam::gnuplotSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "set term postscript color" << nl
        << "set output \"" << points.name() << ".ps\"" << nl;

    // Set secondary Y axis if using two columns. Falls back to same
    // values if both on same scale. However, ignore if more columns.
    if (valueSetNames.size() == 2)
    {
        os  << "set ylabel \"" << valueSetNames[0] << "\"" << nl
            << "set y2label \"" << valueSetNames[1] << "\"" << nl
            << "set ytics nomirror" << nl << "set y2tics" << nl;
    }

    os  << "plot";

    forAll(valueSets, i)
    {
        if (i != 0)
        {
            os << ',';
        }

        os  << " \"-\" title \"" << valueSetNames[i] << "\" with lines";

        if (valueSetNames.size() == 2)
        {
            os  << " axes x1y" << (i+1) ;
        }
    }
    os  << nl;

    forAll(valueSets, i)
    {
        this->writeTable(points, *valueSets[i], os);
        os  << "e" << nl;
    }
}


template<class Type>
void Foam::gnuplotSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& trackPoints,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorInFunction
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }
    if (trackPoints.size() > 0)
    {
        os  << "set term postscript color" << nl
            << "set output \"" << trackPoints[0].name() << ".ps\"" << nl;

        forAll(trackPoints, trackI)
        {
            os  << "plot";

            forAll(valueSets, i)
            {
                if (i != 0)
                {
                    os << ',';
                }

                os  << " \"-\" title \"" << valueSetNames[i] << "\" with lines";
            }
            os << nl;

            forAll(valueSets, i)
            {
                this->writeTable(trackPoints[trackI], valueSets[i][trackI], os);
                os  << "e" << nl;
            }
        }
    }
}


// ************************************************************************* //
