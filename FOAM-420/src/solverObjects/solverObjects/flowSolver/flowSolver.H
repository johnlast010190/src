/*---------------------------------------------------------------------------* \
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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::fv::flowSolver

Description
    Helper class for flow solvers.

SourceFiles
    flowSolver.C

\*---------------------------------------------------------------------------*/
#ifndef flowSolver_H
#define flowSolver_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/include/fvCFD.H"
#include "fluidThermo/fluidThermo.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class flowSolver Declaration
\*---------------------------------------------------------------------------*/

class flowSolver
:
    public solverObject
{

    // Private member functions

        // Correct inactive GIB zones
        // Helper function to zero face field values inside the GIB zone
        // (for example, the face-based buoyancy term)
        void correctInactiveGIBZoneFaces(surfaceScalarField&) const;


protected:

    // Protected member data

        //- Interpolated density field
        autoPtr<surfaceScalarField> rhof_;

        //- Gravitational acceleration [m s^-2]
        autoPtr<uniformDimensionedVectorField> g_;

        //- Reference height [m]
        autoPtr<uniformDimensionedScalarField> hRef_;

        //- g*hRef
        autoPtr<volScalarField> gh_;

        //- g*hRef on the faces
        autoPtr<surfaceScalarField> ghf_;

    // Protected member functions

        //- Return the pressure
        virtual const volScalarField& p() const = 0;

        //- Return the thermo object
        virtual const fluidThermo& thermo() const = 0;

        //- Return the convective density (used for flux)
        virtual const volScalarField& rho() const = 0;

        //- Return the volumetric flux in transonic cases
        virtual const surfaceScalarField& phiv() const = 0;

        //- Interpolation for convective density
        virtual const tmp<surfaceScalarField> rhoInterpolation
        (
            const volScalarField& rho
        ) const;

        //- Interpolation for convected density
        virtual const tmp<surfaceScalarField> rhoInterpolation() const
        {
            return rhoInterpolation(rho());
        }

        //- Interpolation for convected density
        virtual void setOrComputeRhof();

        //- Return the density used for buoyancy (reference to rho() if
        //  distinct buoyancy not in use)
        virtual const volScalarField& buoyantRho() const = 0;

        //- Return solution control
        virtual const solutionControl& solnControl() const = 0;


public:

// Constructors

    // - Construct from components
    flowSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

    //- Destructor
    virtual ~flowSolver()
    {}


    // Public member functions

        //- Initialize gravity and hRef
        void initializeGravityHref();

        //- Calculate gh fields
        void calculateghFields(bool buoyant);

        //- Return face-based buoyancy source (rho*g term)
        virtual tmp<surfaceScalarField> faceBuoyancyForce
        (
            bool includeSnGradP=false
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
