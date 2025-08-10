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

\*---------------------------------------------------------------------------*/

#include "sampledSetWriters/xmgrace/xmgraceSetWriter.H"
#include "coordSet/coordSet.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/IOstreams/Fstreams/OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::xmgraceSetWriter<Type>::xmgraceSetWriter()
:
    writer<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::xmgraceSetWriter<Type>::~xmgraceSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const char* Foam::xmgraceSetWriter<Type>::ext() const {
    return ".agr";
}

template<class Type>
void Foam::xmgraceSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "@g0 on" << nl
        << "@with g0" << nl
        << "@    title \"" << points.name() << '"' << nl
        << "@    xaxis label " << '"' << points.axis() << '"' << nl;

    forAll(valueSets, i)
    {
        os  << "@    s" << i << " legend " << '"'
            << valueSetNames[i] << '"' << nl
            << "@target G0.S" << i << nl;

        this->writeTable(points, *valueSets[i], os);

        os  << '&' << nl;
    }
}


template<class Type>
void Foam::xmgraceSetWriter<Type>::write
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
        os  << "@g0 on" << nl
            << "@with g0" << nl
            << "@    title \"" << trackPoints[0].name() << '"' << nl
            << "@    xaxis label " << '"' << trackPoints[0].axis() << '"' << nl;

        // Data index.
        label sI = 0;

        forAll(trackPoints, trackI)
        {
            forAll(valueSets, i)
            {
                os  << "@    s" << sI << " legend " << '"'
                    << valueSetNames[i] << "_track" << i << '"' << nl
                    << "@target G0.S" << sI << nl;
                this->writeTable(trackPoints[trackI], valueSets[i][trackI], os);
                os  << '&' << nl;

                sI++;
            }
        }
    }
}

// ************************************************************************* //
