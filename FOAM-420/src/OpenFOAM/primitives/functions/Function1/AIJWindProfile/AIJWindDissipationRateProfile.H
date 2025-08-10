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
    Foam::Function1Types::AIJWindDissipationRateProfile

Description
    AIJWindDissipationRateProfile for ABL inlet specification
    Usage:
    \verbatim
    <entryName> <AIJWindDissipationRateProfileFunction>;
    <entryName>Coeffs
    {
        Cmu 0.09;
        Uref 5;
        alpha 0.2;
        Href 10;
        zG 10;
    }
    \endverbatim

    Turbulent dissipation, \epsilon, is derived assuming local equilibrium
    of turbulent production and dissipation for the Power Law profile:

    \f[
    \epsilon = Cmu^0.5*k*frac{U_s}{z_s}\alpha
    *\left(frac{z}{z_s}\right)^(\alpha-1)
    \f]

    with Cmu = 0.09

See also
    Foam::Function1

SourceFiles
    AIJWindDissipationRateProfile.C

\*---------------------------------------------------------------------------*/

#ifndef AIJWindDissipationRateProfile_H
#define AIJWindDissipationRateProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/AIJWindProfile/AIJWindKineticEnergyProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                Class AIJWindDissipationRateProfile Declaration
\*---------------------------------------------------------------------------*/

class AIJWindDissipationRateProfile
:
    public AIJWindKineticEnergyProfile
{
private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const AIJWindDissipationRateProfile&);

public:

  // Runtime type information
    TypeName("AIJWindDissipationRateProfile");

  // Constructors

        //- Construct from entry name and dictionary
        AIJWindDissipationRateProfile
        (
            const word& entryName,
            const dictionary& dict
        );


        //- Copy constructor
        AIJWindDissipationRateProfile
        (
            const AIJWindDissipationRateProfile& wp
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return tmp<Function1<scalar>>
                (
                    new AIJWindDissipationRateProfile(*this)
                );
        }


       //- Destructor
       virtual ~AIJWindDissipationRateProfile();


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
