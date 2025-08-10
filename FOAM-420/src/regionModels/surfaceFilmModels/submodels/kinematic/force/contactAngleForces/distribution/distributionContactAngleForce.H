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
    (c) 2017 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::distributionContactAngleForce

Description
    PDF distribution based film contact angle force

See also
    Foam::regionModels::surfaceFilmModels::contactAngleForce
    Foam::distributionModel

SourceFiles
    distributionContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef distributionContactAngleForce_H
#define distributionContactAngleForce_H

#include "submodels/kinematic/force/contactAngleForces/contactAngleForce/contactAngleForce.H"
#include "distributionModel/distributionModel.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                Class distributionContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class distributionContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> distribution_;


    // Private member functions

        //- Disallow default bitwise copy construct
        distributionContactAngleForce(const distributionContactAngleForce&);

        //- Disallow default bitwise assignment
        void operator=(const distributionContactAngleForce&);


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("distributionContactAngle");


    // Constructors

        //- Construct from surface film model
        distributionContactAngleForce
        (
            surfaceFilmModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distributionContactAngleForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
