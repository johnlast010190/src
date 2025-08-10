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
    Foam::Function1Types::RHWindProfile

Description
    RHWindProfile (RichardsHoxey) for ABL inlet specification
    Usage:
    \verbatim
        <entryName> <RHWindProfileFunction>;
    <entryName>Coeffs
    {
        Cmu 0.09;
        Kappa 0.41;
        Uref 5;
        Href 10;
        z0 0.2;
    }
    \endverbatim

See also
    Foam::Function1

SourceFiles
    RHWindProfile.C

\*---------------------------------------------------------------------------*/

#ifndef RHWindProfile_H
#define RHWindProfile_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class RHWindProfile Declaration
\*---------------------------------------------------------------------------*/

class RHWindProfile
:
    public Function1<scalar>
{
protected:
          //- Cmu constant
        scalar Cmu_;

        //- Kappa constant
        scalar Kappa_;

        //- Reference velocity
        scalar Uref_;

        //- Reference height
        scalar Href_;

        //- Surface roughness height
        scalar z0_;

    // Member Functions

        //- Return frictional velocity
        scalar Ufrictional() const;

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);

        //- Disallow default bitwise assignment
        void operator=(const RHWindProfile&);


public:

  // Runtime type information
    TypeName("RHWindProfile");

  // Constructors

        //- Construct from entry name and dictionary
        RHWindProfile
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        RHWindProfile
        (
            const RHWindProfile& wp
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const override
        {
            return tmp<Function1<scalar>>(new RHWindProfile(*this));
        }

  //- Destructor
        virtual ~RHWindProfile();


  // Member Functions

        //- Return value for coordinate z
        virtual scalar value(const scalar z_) const override;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
