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

#include "interpolations/interpolationTable/tableReaders/openFoam/openFoamTableReader.H"
#include "global/fileOperations/fileOperation/fileOperation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::openFoamTableReader<Type>::openFoamTableReader(const dictionary& dict)
:
    tableReader<Type>(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::openFoamTableReader<Type>::~openFoamTableReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::openFoamTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, Type>>& data
)
{
    // Read data from file
    fileHandler().NewIFstream(fName)()() >> data;
}


template<class Type>
void Foam::openFoamTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& data
)
{
    // Read data from file
    fileHandler().NewIFstream(fName)()() >> data;
}


// ************************************************************************* //
