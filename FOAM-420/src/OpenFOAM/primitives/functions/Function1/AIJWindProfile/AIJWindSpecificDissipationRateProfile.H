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
    Foam::Function1Types::AIJWindSpecificDissipationRateProfile

Description
    AIJWindSpecificDissipationRateProfile for ABL inlet specification
    Usage:
    \verbatim
    <entryName> <AIJWindSpecificDissipationRateProfileFunction>;
    <entryName>Coeffs
    {
        Cmu 0.09;
        Uref 5;
        alpha 0.2;
        Href 10;
        zG 10;
    }
    \endverbatim

    Turbulent specific dissipation, \omega, is derived assuming local equilibrium
    of turbulent production and dissipation for the Power Law profile:

    \f[
    \omega = frac{epsilon}{k*Cmu_}
    \f]

    with Cmu = 0.09

See also
    Foam::Function1

SourceFiles
    AIJWindSpecificDissipationRateProfile.C

\*---------------------------------------------------------------------------*/

#ifndef AIJWindSpecificDissipationRateProfile_H
#define AIJWindSpecificDissipationRateProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/AIJWindProfile/AIJWindKineticEnergyProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
           Class AIJWindSpecificDissipationRateProfile Declaration
\*---------------------------------------------------------------------------*/

class AIJWindSpecificDissipationRateProfile
:
    public AIJWindKineticEnergyProfile
{
private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const AIJWindSpecificDissipationRateProfile&);

public:

  // Runtime type information
    TypeName("AIJWindSpecificDissipationRateProfile");

  // Constructors

        //- Construct from entry name and dictionary
        AIJWindSpecificDissipationRateProfile
        (
            const word& entryName,
            const dictionary& dict
        );


        //- Copy constructor
        AIJWindSpecificDissipationRateProfile
        (
            const AIJWindSpecificDissipationRateProfile& wp
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return tmp<Function1<scalar>>
                (
                    new AIJWindSpecificDissipationRateProfile(*this)
                );
        }


       //- Destructor
       virtual ~AIJWindSpecificDissipationRateProfile();


    // Member Functions

        //- Return value for coordinate z
       virtual scalar value(const scalar z) const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
