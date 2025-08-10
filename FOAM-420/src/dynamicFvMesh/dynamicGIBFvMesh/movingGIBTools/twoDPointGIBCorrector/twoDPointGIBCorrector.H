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
    (c) 2017 Esi Ltd.

Class
    Foam::twoDPointGIBCorrector

Description
    Derived class of twoDPointCorrector with additional functions for GIB

SourceFiles
    twoDPointGIBCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef twoDPointGIBCorrector_H
#define twoDPointGIBCorrector_H

#include "twoDPointCorrector/twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class twoDPointGIBCorrector Declaration
\*---------------------------------------------------------------------------*/

class twoDPointGIBCorrector
:
    public twoDPointCorrector
{
    // Private data

public:

    // Constructors

        //- Construct from components
        twoDPointGIBCorrector(const polyMesh& mesh);

    //- Destructor
    ~twoDPointGIBCorrector(){};

    // Member Functions

        //- Correct motion points based on original
        void correctPoints(pointField&, const pointField&) const;

        //- Correct points moving on planes
        void correctPlanes(pointField& p) const;

        //- Correct motion points
        void correctPointsPlus(pointField& p) const;

        //- Correct motion points
        void correctPointsPlus(pointField&, const pointField&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
