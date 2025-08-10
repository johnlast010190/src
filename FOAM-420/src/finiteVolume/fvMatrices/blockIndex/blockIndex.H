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
    (c) 2019 Esi Ltd.

Class
    Foam::blockIndex

Description
    Class used to find the indexes of the submatrices inside the blockMatrix

SourceFiles
    blockIndex.C

\*---------------------------------------------------------------------------*/

#ifndef blockIndex_H
#define blockIndex_H


#include "primitives/strings/lists/wordList.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blockIndex Declaration
\*---------------------------------------------------------------------------*/

class blockIndex
{
    // Private data

        //- Solving variable field names (U, p)
        wordList names_;

        //- Sizes of the fields
        labelList sizes_;

        //- Starting (matrix) index of each variable inside the matrix
        labelList startIndices_;

        //- Solving variable field names global style (U, U, U, p)
        wordList gNames_;

        // Return the index of a given field
        label index(const word&) const;

public:

    // Constructors

        //- Construct using names and sizes
        blockIndex(const wordList&, const labelList&);

        //- Copy constructor
        blockIndex(const blockIndex&);

        //- Null constructor
        blockIndex();

    // Member Functions

        // Check if the class holds valid information
        void check();

        // Return the index of a given field
        label findIndex(const word&) const;

        // Return the starting (matrix) index of a given field
        label findStartIndex(const word&) const;

        // Returns: true/false if change/not change rows
        // needed from the blockMatrixClass at the coupling sub-matrices
        bool increment(const label, const label) const;

        // provide global row index and return the name of the variable
        word variableName(const label) const;

    // Operators

        void operator=(const blockIndex&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
