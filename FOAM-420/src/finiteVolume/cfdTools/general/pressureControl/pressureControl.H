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
    Foam::pressureControl

Description
    Provides controls for the pressure reference is closed-volume simulations
    and a general method for limiting the pressure during the startup of
    steady-state simulations.

SourceFiles
    pressureControlI.H
    pressureControl.C

\*---------------------------------------------------------------------------*/

#ifndef pressureControl_H
#define pressureControl_H

#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pressureControl Declaration
\*---------------------------------------------------------------------------*/

class pressureControl
{
    // Private data

        //- Optional cell in which the reference pressure is set
        label refCell_;

        //- Optional pressure reference level
        scalar refValue_;

        //- Pressure lower-limit
        dimensionedScalar pMax_;

        //- Pressure upper-limit
        dimensionedScalar pMin_;

        //- Pressure lower-limit
        bool limitMaxP_;

        //- Pressure upper-limit
        bool limitMinP_;


public:

    // Constructors

        //- Construct from the SIMPLE/PIMPLE sub-dictionary
        pressureControl
        (
            const volScalarField& p,
            const volScalarField& rho,
            const dictionary& dict,
            const scalar& abspOffset = 0.0,
            const bool pRefRequired = true
        );

        //- Construct from the SIMPLE/PIMPLE sub-dictionary
        //  Controlling pressure from other p fields
        pressureControl
        (
            const volScalarField& p,
            const volScalarField& pRef,
            const volScalarField& rho,
            const dictionary& dict,
            const scalar& abspOffset = 0.0,
            const bool pRefRequired = true
        );


    // Member Functions

        //- Return the cell in which the reference pressure is set
        inline label refCell() const;

        //- Return the pressure reference level
        inline scalar refValue() const;

        //- Limit the pressure if necessary and return true if so
        bool limit(volScalarField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/pressureControl/pressureControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
