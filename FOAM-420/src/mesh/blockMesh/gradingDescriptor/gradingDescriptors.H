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
    Foam::gradingDescriptors

Description
    List of gradingDescriptor for the sections of a block with additional IO
    functionality

SourceFiles
    gradingDescriptors.C

\*---------------------------------------------------------------------------*/

#ifndef gradingDescriptors_H
#define gradingDescriptors_H

#include "gradingDescriptor/gradingDescriptor.H"
#include "containers/Lists/List/List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class gradingDescriptors;
Istream& operator>>(Istream&, gradingDescriptors&);

/*---------------------------------------------------------------------------*\
                     Class gradingDescriptors Declaration
\*---------------------------------------------------------------------------*/

class gradingDescriptors
:
    public List<gradingDescriptor>
{

public:

    // Constructors

        //- Default constructor
        gradingDescriptors();

        //- Construct from a gradingDescriptor
        gradingDescriptors(const gradingDescriptor& gd);


    // Member functions

        //- Return the inverse gradingDescriptors with 1/expansionRatio
        gradingDescriptors inv() const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, gradingDescriptors&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
