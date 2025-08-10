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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::turbulentDispersionModels::Gosman

Description
    Turbulent dispersion model of Gosman et al.

    Reference:
    \verbatim
        "Multidimensional modeling of turbulent two-phase flows in stirred
         vessels"
        Gosman, A.D., Lekakou, C., Politis, S., Issa, R.I., and Looney, M.K.,
        AIChE Journal
        Volume 38, Issue 12, 1992, pp. 1946-1956
     \endverbatim

SourceFiles
    Gosman.C

\*---------------------------------------------------------------------------*/

#ifndef Gosman_H
#define Gosman_H

#include "disperseEulerian/forceModels/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace decoupledEulerian
{

class phase;
class dragModel;

/*---------------------------------------------------------------------------*\
                             Class Gosman Declaration
\*---------------------------------------------------------------------------*/

class Gosman
:
    public turbulentDispersionModel
{
    // Private data

        //- turbulent Schmidt number
        const dimensionedScalar Sct_;


public:

    //- Runtime type information
    TypeName("Gosman");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Gosman
        (
            const dictionary& dict,
            const phase& phase,
            const bool registerObject
        );


    //- Destructor
    virtual ~Gosman();


    // Member Functions

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decoupledEulerian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
