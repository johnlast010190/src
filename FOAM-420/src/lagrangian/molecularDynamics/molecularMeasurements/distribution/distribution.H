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
    Foam::distribution

Description
    Accumulating histogram of values.  Specified bin resolution
    automatic generation of bins.

SourceFiles
    distributionI.H
    distribution.C

\*---------------------------------------------------------------------------*/

#ifndef distribution_H
#define distribution_H

#include "containers/HashTables/Map/Map.H"
#include "primitives/Pair/Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class distribution;

Ostream& operator<<(Ostream&, const distribution&);


/*---------------------------------------------------------------------------*\
                        Class distribution Declaration
\*---------------------------------------------------------------------------*/

class distribution
:
    public Map<label>
{
    // Private data

        scalar binWidth_;


public:

    //- Runtime type information

        TypeName("distribution");

    // Static functions

        //- Write to file

            static void write
            (
                const fileName& file,
                const List<Pair<scalar>>& pairs
            );


    // Constructors

        //- Construct null
        distribution();

        //- Construct from binWidth
        distribution(const scalar binWidth);

        //- Construct as copy
        distribution(const distribution&);


    //- Destructor
    virtual ~distribution();


    // Member Functions

        label totalEntries() const;

        scalar approxTotalEntries() const;

        scalar mean() const;

        scalar median();

        //- Add a value to the appropriate bin of the distribution.
        void add(const scalar valueToAdd);

        void add(const label valueToAdd);

        void insertMissingKeys();

        List<Pair<scalar>> normalised();

        List<Pair<scalar>> normalisedMinusMean();

        List<Pair<scalar>> normalisedShifted(scalar shiftValue);

        List<Pair<scalar>> raw();


        // Access

            inline scalar binWidth() const;


    // Member Operators

        void operator=(const distribution&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const distribution&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "distribution/distributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
