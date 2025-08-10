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
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::noPhaseChange

Description
    Dummy phase change model for 'none'

SourceFiles
    noPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef noPhaseChange_H
#define noPhaseChange_H

#include "submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class noPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class noPhaseChange
:
    public phaseChangeModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        noPhaseChange(const noPhaseChange&);

        //- Disallow default bitwise assignment
        void operator=(const noPhaseChange&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from surface film model
        noPhaseChange(surfaceFilmModel& film, const dictionary& dict);


    //- Destructor
    virtual ~noPhaseChange();


    // Member Functions

        // Evolution

            //- Return the model 'active' status - default active = true
            virtual bool active() const;

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
