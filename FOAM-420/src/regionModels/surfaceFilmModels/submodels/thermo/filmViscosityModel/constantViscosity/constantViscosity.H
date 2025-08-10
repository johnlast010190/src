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
    Foam::regionModels::surfaceFilmModels::constantViscosity

Description
    Constant viscosity model

SourceFiles
    constantViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef constantViscosity_H
#define constantViscosity_H

#include "submodels/thermo/filmViscosityModel/filmViscosityModel/filmViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class constantViscosity Declaration
\*---------------------------------------------------------------------------*/

class constantViscosity
:
    public filmViscosityModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        constantViscosity(const constantViscosity&);

        //- Disallow default bitwise assignment
        void operator=(const constantViscosity&);


protected:

    // Protected data

        //- Constant viscosity [Pa.s]
        dimensionedScalar mu0_;



public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from surface film model
        constantViscosity
        (
            surfaceFilmModel& film,
            const dictionary& dict,
            volScalarField& mu
        );


    //- Destructor
    virtual ~constantViscosity();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const volScalarField& p,
                const volScalarField& T
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
