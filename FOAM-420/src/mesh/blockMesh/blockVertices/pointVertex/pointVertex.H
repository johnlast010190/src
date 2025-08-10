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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::blockVertices::pointVertex

Description

SourceFiles
    pointVertex.C

\*---------------------------------------------------------------------------*/

#ifndef blockVertices_pointVertex_H
#define blockVertices_pointVertex_H

#include "blockVertices/blockVertex/blockVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockVertices
{

/*---------------------------------------------------------------------------*\
                           Class pointVertex Declaration
\*---------------------------------------------------------------------------*/

class pointVertex
:
    public blockVertex
{
protected:

    // Protected member data

        //- The vertex location
        point vertex_;


public:

    //- Runtime type information
    TypeName("point");


    // Constructors

        //- Construct from Istream setting pointsList
        pointVertex
        (
            const dictionary&,
            const label index,
            const searchableSurfaces& geometry,
            Istream&
        );


    //- Destructor
    virtual ~pointVertex()
    {}


    // Member Functions

        virtual operator point() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockVertices
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
