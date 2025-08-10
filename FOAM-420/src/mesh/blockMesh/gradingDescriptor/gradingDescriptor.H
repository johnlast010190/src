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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::gradingDescriptor

Description
    Handles the specification for grading within a section of a block

    blockFraction: the fraction of the block the section occupies

    nDivFraction: the fraction of the divisions of the block allocated to
        the section

    expansionRatio: the expansions ratio for the grading with the section of
        block defined as the ratio of the size of the division at either and
        of the section.

SourceFiles
    gradingDescriptor.C

\*---------------------------------------------------------------------------*/

#ifndef gradingDescriptor_H
#define gradingDescriptor_H

#include "primitives/Scalar/scalar/scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class gradingDescriptor;
class gradingDescriptors;
Istream& operator>>(Istream&, gradingDescriptor&);
Ostream& operator<<(Ostream&, const gradingDescriptor&);

/*---------------------------------------------------------------------------*\
                     Class gradingDescriptor Declaration
\*---------------------------------------------------------------------------*/

class gradingDescriptor
{
    // Private data

        scalar blockFraction_;
        scalar nDivFraction_;
        scalar expansionRatio_;


public:

    friend class gradingDescriptors;


    // Constructors

        //- Default constructor
        gradingDescriptor();

        //- Construct from components
        gradingDescriptor
        (
            const scalar blockFraction,
            const scalar nDivFraction,
            const scalar expansionRatio
        );

        //- Construct from expansionRatio
        gradingDescriptor
        (
            const scalar expansionRatio
        );

        //- Construct from Istream
        gradingDescriptor(Istream&);


    //- Destructor
    ~gradingDescriptor();


    // Member Functions

        // Access

            scalar blockFraction() const
            {
                return blockFraction_;
            }

            scalar nDivFraction() const
            {
                return nDivFraction_;
            }

            scalar expansionRatio() const
            {
                return expansionRatio_;
            }


    // Member functions

        //- Return the inverse gradingDescriptor with 1/expansionRatio
        gradingDescriptor inv() const;


    // Member operators

        bool operator==(const gradingDescriptor&) const;
        bool operator!=(const gradingDescriptor&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, gradingDescriptor&);
        friend Ostream& operator<<(Ostream&, const gradingDescriptor&);

        friend Istream& operator>>(Istream&, gradingDescriptors&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
