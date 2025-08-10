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
    (c) 2015-2017 OpenFOAM Foundation

Class
    Foam::fv::patchRotatingVelocityForce

Description
    Calculates and applies the force necessary to maintain the specified mean
    velocity averaged over the specified patch.

    Note: Currently only handles kinematic pressure (incompressible solvers).

Usage
    Example usage:
    \verbatim
    selectionMode   all;                    // Apply force to all cells
    fields          (U);                    // Name of velocity field
    patch           inlet;                  // Name of the patch
    Ubar            (10.0 0 0);             // Desired mean velocity
    relaxation      0.2;                    // Optional relaxation factor
    \endverbatim

SourceFiles
    patchRotatingVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef patchRotatingVelocityForce_H
#define patchRotatingVelocityForce_H

#include "sources/derived/rotatingVelocityForce/rotatingVelocityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class patchRotatingVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class patchRotatingVelocityForce
:
    public rotatingVelocityForce
{
protected:

    // Protected data

        //- Patch name
        word patch_;

        //- Patch index
        label patchi_;


    // Protected Member Functions

        //- Calculate and return the magnitude of the mean velocity
        //  averaged over the specified patch
        virtual scalar omegaAve(const volVectorField& U) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchRotatingVelocityForce(const patchRotatingVelocityForce&);

        //- Disallow default bitwise assignment
        void operator=(const patchRotatingVelocityForce&);


public:

    //- Runtime type information
    TypeName("patchRotatingVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        patchRotatingVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
