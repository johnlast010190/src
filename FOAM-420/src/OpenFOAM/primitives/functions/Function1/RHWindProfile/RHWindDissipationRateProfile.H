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
    (c) 2020 Esi Ltd.

Class
    Foam::Function1Types::RHWindDissipationRateProfile

Description
    RHWindDissipationRateProfile (RichardsHoxey) for ABL inlet specification

See also
    Foam::Function1

SourceFiles
    RHWindDissipationRateProfile.C

\*---------------------------------------------------------------------------*/

#ifndef RHWindDissipationRateProfile_H
#define RHWindDissipationRateProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/RHWindProfile/RHWindKineticEnergyProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                Class RHWindDissipationRateProfile Declaration
\*---------------------------------------------------------------------------*/

class RHWindDissipationRateProfile
:
    public RHWindKineticEnergyProfile
{
private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const RHWindDissipationRateProfile&);

public:

  // Runtime type information
    TypeName("RHWindDissipationRateProfile");

  // Constructors

        //- Construct from entry name and dictionary
        RHWindDissipationRateProfile
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        RHWindDissipationRateProfile
        (
            const RHWindDissipationRateProfile& wp
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return tmp<Function1<scalar>>
                (
                    new RHWindDissipationRateProfile(*this)
                );
        }


       //- Destructor
       virtual ~RHWindDissipationRateProfile();


    // Member Functions

        //- Return value for coordinate z
       virtual scalar value(const scalar z_) const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
