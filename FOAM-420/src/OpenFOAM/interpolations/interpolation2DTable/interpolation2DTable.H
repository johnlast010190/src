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
    (c) 2017-2023 Esi Ltd.

Class
    Foam::interpolation2DTable

Description
    2D table interpolation. The data must be in ascending order in both
    dimensions x and y.

SourceFiles
    interpolation2DTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolation2DTable_H
#define interpolation2DTable_H

#include "containers/Lists/List/List.H"
#include "primitives/Tuple2/Tuple2.H"
#include "interpolations/interpolationTable/tableReaders/tableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interpolation2DTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolation2DTable
:
    public List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          //!< Exit with a FatalError
            WARN,           //!< Issue warning and clamp value (default)
            CLAMP           //!< Clamp value to the start/end value
        };

        //- Convenience typedef
        typedef List<Tuple2<scalar, List<Tuple2<scalar, Type>>>> table;


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- File name
        fileName fileName_;

        //- The actual reader
        autoPtr<tableReader<Type>> reader_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        //- Return interpolated value in List
        Type interpolateValue
        (
            const List<Tuple2<scalar, Type>>& data,
            const scalar
        ) const;

        //- Return an X index from the matrix
        template<class BinaryOp>
        label Xi
        (
            const BinaryOp& bop,
            const scalar x,
            const bool reverse
        ) const;


public:

    // Constructors

        //- Construct null
        interpolation2DTable();

        //- Construct from components
        interpolation2DTable
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
            const boundsHandling bounds,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        interpolation2DTable(const fileName& fName);

        //- Construct by reading the fileName and boundsHandling from dictionary
        interpolation2DTable(const dictionary& dict);

        //- Construct copy
        interpolation2DTable(const interpolation2DTable& interpTable);


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        static boundsHandling wordToBoundsHandling(const word& bound);

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void checkOrder() const;

        //- Check if the values in the columns are matching
        void checkTableShapeAndValues() const;

        //- Write
        void write(Ostream& os) const;

        //- Scale table data by a scalar value
        inline void scale(const scalar& Scale);

        //- Interpolated daxpy operation
        //  table addition y = a*y + b*x where x is
        //  interpolated to y and added
        inline void idaxpy
        (
            const scalar& a,
            const scalar& b,
            const interpolation2DTable<Type>& x
        );

        //- Integrating over y keeping x constant
        inline Type integrateOverY
        (
            const scalar x,
            const scalar y0,
            const scalar y1
        ) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const List<Tuple2<scalar, Type>>& operator[](const label) const;

        //- Assignment operator
        //  (resizes this data table to size of input table)
        inline void operator=
        (
            const interpolation2DTable<Type>& interpolationTable
        );

        //- Return an interpolated value
        Type operator()(const scalar x, const scalar y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolations/interpolation2DTable/interpolation2DTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
