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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::solidBodyMotionFunctions::SDA

Description
    Ship design analysis (SDA) 3DoF motion function.

    Comprising sinusoidal roll (rotation about x), heave (z-translation)
    and sway (y-translation) motions with changing amplitude and phase.

See also
    Foam::solidBodyMotionFunctions::tabulated6DoFMotion

SourceFiles
    SDA.C

\*---------------------------------------------------------------------------*/

#ifndef SDA_H
#define SDA_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class SDA Declaration
\*---------------------------------------------------------------------------*/

class SDA
:
    public solidBodyMotionFunction
{
    // Private data

        //- Center of gravity
        vector CofG_;

        //- Model scale ratio
        scalar lamda_;

        //- Max roll amplitude [rad]
        scalar rollAmax_;

        //- Min roll amplitude [rad]
        scalar rollAmin_;

        //- Heave amplitude [m]
        scalar heaveA_;

        //- Sway amplitude [m]
        scalar swayA_;

        //- Damping Coefficient [-]
        scalar Q_;

        //- Time Period for liquid [sec]
        scalar Tp_;

        //- Natural Period of Ship [sec]
        scalar Tpn_;

        //- Reference time step [sec]
        scalar dTi_;

        //- Incr. in Tp/unit 'dTi'[-]
        scalar dTp_;


    // Private Member Functions

        //- Disallow copy construct
        SDA(const SDA&);

        //- Disallow default bitwise assignment
        void operator=(const SDA&);


public:

    //- Runtime type information
    TypeName("SDA");


    // Constructors

        //- Construct from components
        SDA
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        SDA
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new SDA
                (
                    SBMFCoeffs_,
                    time_,
                    frameName_
                )
            );
        }


    //- Destructor
    virtual ~SDA();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
