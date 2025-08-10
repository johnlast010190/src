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
    (c) 2018 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::fv::rotatingVelocityForce

Group
    grpFvOptionsSources

Description
    Calculates and applies the force necessary to maintain the specified
    mean angular velocity.

Usage
    Example usage:
    \verbatim
    selectionMode   all;                    // Apply force to all cells
    fields          (U);                    // Name of velocity field
    origin          (0 0 0);                // Origin for rotation
    axis            (1 0 0);                // Axis of rotation
    omega           10;                     // Desired mean angular velocity
    relaxation      0.2;                    // Optional relaxation factor
    \endverbatim

SourceFiles
    rotatingVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingVelocityForce_H
#define rotatingVelocityForce_H

#include "memory/autoPtr/autoPtr.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "cellSetOption/cellSetOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class rotatingVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class rotatingVelocityForce
:
    public cellSetOption
{
protected:

    // Protected data

        //- Target average angular velocity in zone
        autoPtr<Function1<scalar>> omega_;

        //- rotating axis
        vector axis_;

         //- Rotional origin
        vector origin_;

        //- Angular force coefficient before correction
        scalar Comega0_;

        //- Change in angular force coefficient
        scalar dComega_;

        //- Relaxation factor
        scalar relaxation_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> rAPtr_;


    // Protected Member Functions

        //- Calculate and return the magnitude of the mean
        //  angular velocity averaged over the selected cellSet
        virtual scalar omegaAve(const volVectorField& U) const;

        //- Velocity field in forcing region based on target rotation
        virtual tmp<Foam::vectorField> targetVelocity() const;

        //- Write the source coefficient to file (for restarts etc)
        void writeProps(const scalar Crf) const;

        //- Correct driving force for a specified rotating velocity
        void update(fvMatrix<vector>& eqn);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        rotatingVelocityForce(const rotatingVelocityForce&);

        //- Disallow default bitwise assignment
        void operator=(const rotatingVelocityForce&);


public:

    //- Runtime type information
    TypeName("rotatingVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        rotatingVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        // Evaluate

            //- Correct the source coefficient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );


            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
