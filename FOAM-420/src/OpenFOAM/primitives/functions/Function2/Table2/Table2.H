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
    (c) 2022 Esi Ltd.

Class
    Foam::Function2s::Table

Description

    Calculating any tabulated properties:

    // Option 1
    {
        <func1Name>     table;
        file            "$FOAM_CASE/constant/mu.dat";
        outOfBounds     clamp;
    }

    // Option 2
    {
        <func1Name>      table
        (
            (
                1e5
                (
                    (300    A)
                    (450    B)
                    (600    C)
                )
            )
            (
                1.1e5
                (
                    (300    D)
                    (450    E)
                    (600    F)
                )
            )
            (
                2e5
                (
                    (300    G)
                    (450    H)
                    (600    I)
                )
            )
        );
        outOfBounds     clamp;
    }

    The tables are described in terms of
    the field1 and field2

    For example, the table:

                           field2
                       |  300   450   600
                  -----+-----------------
                  1e5  |   A     B     C
        field1  1.1e5  |   D     E     F
                  2e5  |   G     H     I


    Is specified by the following:

       (
            (
                1e5
                (
                    (300    A)
                    (450    B)
                    (600    C)
                )
            )
            (
                1.1e5
                (
                    (300    D)
                    (450    E)
                    (600    F)
                )
            )
            (
                2e5
                (
                    (300    G)
                    (450    H)
                    (600    I)
                )
            )
        )

    All field1 and field2 must be in ascending order
    in each list.

\*---------------------------------------------------------------------------*/

#ifndef Table2_H
#define Table2_H

#include "primitives/functions/Function2/Function2/Function2.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class Table Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Table
:
    public FieldFunction2<Type, Table<Type>>
{
    // Private member data

        //- Table values
        interpolation2DTable<Type> values_;

        //- Enumeration for handling out-of-bound values
        word bounds_;

        //- File name
        fileName fileName_;


    // Private member functions

        //- Explicit reader
        interpolation2DTable<Type> readValues
        (
            const word& name,
            const dictionary& dict
        );


public:

    //- Runtime type information
    TypeName("table");


    // Constructors

        //- Construct from name and dictionary
        Table(const word& name, const dictionary& dict);


    // Member Functions

        //- Return the non-uniform table of values
        const interpolation2DTable<Type>& values() const
        {
            return values_;
        }

        //- Evaluate the function and return the result
        virtual Type value(scalar x, scalar y) const;

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const
        {
            return tmp<Function2<Type>>(new Table(*this));
        }

        //- Write the function coefficients
        void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function2/Table2/Table2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
