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

Class
    Foam::BinSum

Description
    Sums into bins

SourceFiles
    BinSum.C

\*---------------------------------------------------------------------------*/

#ifndef BinSum_H
#define BinSum_H

#include "primitives/ops/ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                          Class BinSum Declaration
\*---------------------------------------------------------------------------*/

template
<
    class IndexType,
    class List,
    class CombineOp = plusEqOp<typename List::value_type>
>
class BinSum
:
    public List
{
    // Private data

        const IndexType min_;

        const IndexType max_;

        const IndexType delta_;


        //- Sum < lowest bin
        typename List::value_type lowSum_;

        //- Sum of >= highest bin
        typename List::value_type highSum_;

public:

    // Constructors

        //- Construct given min, max, delta
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta
        );

        //- Construct given min, max, delta and data
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta,
            const UList<IndexType>& indexVals,
            const List& vals,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );


        // Access

            //- Return the delta
            inline IndexType delta() const
            {
                return delta_;
            }

            //- Return the sum of all added elements < min
            inline const IndexType& lowSum() const
            {
                return lowSum_;
            }

            //- Return the sum of all added elements >= max
            inline const IndexType& highSum() const
            {
                return highSum_;
            }

            void add
            (
                const IndexType& indexVal,
                const typename List::const_reference val,
                const CombineOp& cop = plusEqOp<typename List::value_type>()
            );

            void add
            (
                const UList<IndexType>& indexVals,
                const List& vals,
                const CombineOp& cop = plusEqOp<typename List::value_type>()
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "containers/Lists/BinSum/BinSum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
