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
    (c) 2020-2023 Esi Ltd.

SourceFiles
    frameSwitches.C
\*---------------------------------------------------------------------------*/

#ifndef frameSwitches_H
#define frameSwitches_H

#include "primitives/bools/Switch/Switch.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class frameSwitches Declaration
\*---------------------------------------------------------------------------*/

class frameSwitches
{
    // Private Data Members

        //- Should frame use acceleration?
        const Switch acceleration_;

        //- Should frame use conservative formulation?
        const Switch conservative_;

        //- Should frame use ddtPhiCor?
        const Switch ddtPhiCorr_;


 public:


    // Constructors

        //- Construct from dictionary
        frameSwitches(const dictionary& dict);


    // Member Functions

        //- Should it use acceleration?
        inline Switch acceleration() const
        {
            return acceleration_;
        }

        //- Should it use conservative formulation?
        inline Switch conservative() const
        {
            return conservative_;
        }

        //- Should it use ddtPhiCor?
        inline Switch ddtPhiCorr() const
        {
            return ddtPhiCorr_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //