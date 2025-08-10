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
    Foam::regionModels::surfaceFilmModels::thermocapillaryForce

Description
    Thermocapillary force

SourceFiles
    thermocapillaryForce.C

\*---------------------------------------------------------------------------*/

#ifndef thermocapillaryForce_H
#define thermocapillaryForce_H

#include "submodels/kinematic/force/force/force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                  Class thermocapillaryForce Declaration
\*---------------------------------------------------------------------------*/

class thermocapillaryForce
:
    public force
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermocapillaryForce(const thermocapillaryForce&);

        //- Disallow default bitwise assignment
        void operator=(const thermocapillaryForce&);


public:

    //- Runtime type information
    TypeName("thermocapillary");


    // Constructors

        //- Construct from surface film model
        thermocapillaryForce
        (
            surfaceFilmModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermocapillaryForce();


    // Member Functions

        // Evolution

            //- Correct
            virtual tmp<fvVectorMatrix> correct(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
