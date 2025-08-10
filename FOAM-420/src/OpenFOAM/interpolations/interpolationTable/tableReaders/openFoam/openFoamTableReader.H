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

Class
    Foam::openFoamTableReader

Description
    Reads an interpolation table from a file - OpenFOAM-format

SourceFiles
    openFoamTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef openFoamTableReader_H
#define openFoamTableReader_H

#include "interpolations/interpolationTable/tableReaders/tableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class openFoamTableReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class openFoamTableReader
:
    public tableReader<Type>
{

public:

    //- Runtime type information
    TypeName("openFoam");


    // Constructors

        //- Construct from dictionary
        openFoamTableReader(const dictionary& dict);

        //- Construct and return a copy
        virtual autoPtr<tableReader<Type>> clone() const
        {
            return autoPtr<tableReader<Type>>
            (
                new openFoamTableReader<Type>
                (
                    *this
                )
            );
        }


    //- Destructor
    virtual ~openFoamTableReader();


    // Member functions

        //- Read the table
        virtual void operator()(const fileName&, List<Tuple2<scalar, Type>> &);

        //- Read 2D table
        virtual void operator()
        (
            const fileName&,
            List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolations/interpolationTable/tableReaders/openFoam/openFoamTableReader.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
