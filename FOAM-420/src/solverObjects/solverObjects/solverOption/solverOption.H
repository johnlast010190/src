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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::fv::solverOption

Description
    Pure virtual base classs for templated SolverOption class

\*---------------------------------------------------------------------------*/

#ifndef solverOption_H
#define solverOption_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "solverObject/solverObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class SolverOption Declaration
\*---------------------------------------------------------------------------*/

class solverOption
:
    public option
{

public:

    // Constructor

        //- Construct from components
        solverOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        )
        :
            option(name, modelType, dict, obr)
        {}

    // Destructor

        virtual ~solverOption()
        {}

    // Member Functions

        virtual class solverObject& solverObject() = 0;

        virtual void setSolveController(const foamSolve& solveController) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
