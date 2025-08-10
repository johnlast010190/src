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
    (c) 2013-2017 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::solidification

Description
    Solidification phase change model where all film mass is converted when the
    local temperature > activation temperature.  The latent heat is
    assumed to be removed by heat-transfer to the wall.

SourceFiles
    solidification.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationPhaseChange_H
#define solidificationPhaseChange_H

#include "submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class solidification Declaration
\*---------------------------------------------------------------------------*/

class solidificationPhaseChange
:
    public phaseChangeModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        solidificationPhaseChange(const solidificationPhaseChange&);

        //- Disallow default bitwise assignment
        void operator=(const solidificationPhaseChange&);


protected:

    // Protected data

        //- Temperature at which solidification starts
        scalar T0_;

        //- Solidification limiter
        //  Maximum fraction of film which can solidify in a time-step
        scalar maxSolidificationFrac_;

        //- Solidification limiter
        //  Maximum rate at which the film can solidify
        dimensionedScalar maxSolidificationRate_;

        //- Accumulated solid mass [kg]
        volScalarField mass_;

        //- Accumulated solid thickness [m]
        volScalarField thickness_;


public:

    //- Runtime type information
    TypeName("solidification");


    // Constructors

        //- Construct from surface film model
        solidificationPhaseChange(surfaceFilmModel& film, const dictionary& dict);


    //- Destructor
    virtual ~solidificationPhaseChange();


    // Member Functions

        // Evolution

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
