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
    Foam::fv::settlingVelocityTakacs

Description

    solverObject to create, store and compute a settling flux to be used
    in an additional passive scalar transport solverObject.

SourceFiles
    settlingVelocityTakacs.C

\*---------------------------------------------------------------------------*/
#ifndef settlingVelocityTakacs_H
#define settlingVelocityTakacs_H

#include "cellSetOption/cellSetOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class settlingVelocityTakacs Declaration
\*---------------------------------------------------------------------------*/

class settlingVelocityTakacs
:
    public cellSetOption
{
    // Private data

        //- Concentration field name
        word XName_;

        //- Flux field names
        word phiName_;
        word UName_;
        word phaseName_;

        //- Model parameter
        dimensionedVector V0_;
        dimensionedVector V0t_;
        dimensionedScalar rh_;
        dimensionedScalar rp_;
        dimensionedScalar rt_;
        dimensionedScalar Xns_;
        dimensionedScalar X34_;

        //- Scaling parameter to get X in [g/L]
        scalar fac_;
        scalar blend_;

        //- use constant drift y/n
        bool constant_;
        //- limit phasic drift
        bool limitBulk_;
        const scalar threshold_;
        const scalar limitHeight_;

        //- Flux field pointer
        autoPtr<surfaceScalarField> phiPtr_;
        autoPtr<volVectorField> UPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        settlingVelocityTakacs(const settlingVelocityTakacs&);

        //- Disallow default bitwise assignment
        void operator=(const settlingVelocityTakacs&);

        //- compute settling flux
        void computeFlux(const fvMesh& mesh);

        //- limit to selected cells
        void limitCellSet();


public:

    //- Runtime type information
    TypeName("settlingVelocityTakacs");


    // Constructors

        //- Construct from components
        settlingVelocityTakacs
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;

        // Evaluate

            //- Add settling advection term to scalar transport equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
