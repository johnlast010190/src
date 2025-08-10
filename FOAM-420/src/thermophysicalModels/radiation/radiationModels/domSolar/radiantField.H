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
    (c) Esi Ltd

Class
    Foam::radiantField

Description
    A storage class for basic radiant field.

    Example of the object specification:
    @verbatim
    radiantField
    {
        radiantInputType userDefined; //solarCalculator

        intensity constant 10; //radiation intensity for userDefined
        direction constant (0 -1 0); //radiation direction for userDefined

        //for solarCalculator input see
        // radiation/submodels/solarCalculator/solarCalculator.H

    }
    @endverbatim

SourceFiles
    radiantField.C

\*---------------------------------------------------------------------------*/

#ifndef radiantField_H
#define radiantField_H

#include "fields/volFields/volFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "submodels/solarCalculator/solarCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class radiantField Declaration
\*---------------------------------------------------------------------------*/

class radiantField
:
    public volScalarField
{

public:

        enum radiantInputType
        {
            ritUserDefined,
            ritSolarCalculator
        };

        static const NamedEnum<radiantInputType, 2> radiantInputTypeNames_;

private:
        //- direction/intensity input mode - valid entries are:
        //  userDefined, solarCalculater
        radiantInputType inputMode_;

        //For now just store both input options as autoPtrs
        //run time selection should be added if options expand

            //- radiation intensity
            autoPtr<Function1<scalar>> radiationFlux_;

            //- radiation direction vector
            autoPtr<Function1<vector>> radiationDirection_;

            //- solar calculater
            autoPtr<solarCalculator> solarCalc_;

        //- patch source coeffs
        const FieldField<Field, scalar>& patchSourceCoeffs_;

        //- current solar flux
        scalar currentFlux_;

        //- current solar direction
        vector currentDirection_;

        //- Update solar tracking index
        label updateTimeIndex_;

        //- initialisation switch
        bool updateAfterRead_;

        //- initialisation switch
        bool initializeField_;

        //- initial tolerance for explicit wave solve
        scalar initTol_;

        //- initialisation matrix solver switch
        label nSolverSolutions_;

        //- matrix solver switch
        bool solveField_;

        //- inlet outlet scaling
        scalar qScale_;


    // Private functions

        //- set fv dicts if matrix solvers to be used
        void setFvDicts();

        //- Provide boundary conditions for intensity field
        const dictionary radiantBoundaryConditions
        (
            const fvMesh& mesh,
            word
        ) const;

        //- create initial intensity field
        tmp<volScalarField> autoCreateRadiantField
        (
            const word& fieldName,
            const fvMesh& mesh
        );

        //- calculate scaling for incedent radiation to match incoming
        //- radiation transmission
        scalar inletOutletScaling();

        //- update radiant power flux and direction
        void updateRadiationProperties();

        //- add diffuse radiation over patches
        //  based on ASHRAE hanbook
        void updateSkyDiffusiveRadiation(volScalarField&);

        //- return vertical direction
        vector verticalDirection();

        //- initialise radiant field and input functions
        void initialise(const dictionary&);

        //- Disallow default bitwise copy construct
        radiantField(const radiantField&);

        //- Disallow default bitwise assignment
        void operator=(const radiantField&);


public:

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        radiantField
        (
            const word& fieldName,
            const fvMesh&,
            const dictionary&,
            const FieldField<Field, scalar>& transmissivity
        );


    // Destructor

        virtual ~radiantField();


    // Member Functions

        //- Update solver
        void calculate(volScalarField&);

        scalar getZenithAngle() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
