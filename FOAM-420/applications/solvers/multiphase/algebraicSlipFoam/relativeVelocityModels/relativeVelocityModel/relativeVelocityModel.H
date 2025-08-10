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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::relativeVelocityModel

Description

SourceFiles
    relativeVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVelocityModel_H
#define relativeVelocityModel_H

#include "cfdTools/general/include/fvCFD.H"
#include "db/dictionary/dictionary.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relativeVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class relativeVelocityModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        relativeVelocityModel(const relativeVelocityModel&);

        //- Disallow default bitwise assignment
        void operator=(const relativeVelocityModel&);

protected:

    // Protected data

        //- Mixture properties
        const immiscibleIncompressibleTwoPhaseMixture& mixture_;

        //- Name of the continuous phase
        const word continuousPhaseName_;

        //- Continuous phase fraction
        const volScalarField& alphac_;

        //- Dispersed phase fraction
        const volScalarField& alphad_;

        //- Continuous density
        const dimensionedScalar& rhoc_;

        //- Dispersed density
        const dimensionedScalar& rhod_;

        //- Relative phase velocity
        mutable volVectorField Urel_;


public:

    //- Runtime type information
    TypeName("relativeVelocityModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relativeVelocityModel,
        dictionary,
        (const dictionary& dict,
         const immiscibleIncompressibleTwoPhaseMixture& mixture),
        (dict, mixture)
    );


    // Constructors

        //- Construct from components
        relativeVelocityModel
        (
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        );


    // Selector
    static autoPtr<relativeVelocityModel> New
    (
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture
    );


    //- Destructor
    virtual ~relativeVelocityModel();


    // Member Functions

        //- Mixture properties
        const immiscibleIncompressibleTwoPhaseMixture& mixture() const
        {
            return mixture_;
        }

        //- Return the relative velocity (disperse phase - continuous phase)
        const volVectorField& Urel() const
        {
            return Urel_;
        }

        //- Return the stress tensor due to the phase transport
        tmp<volSymmTensorField> tauDm() const;

        //- Update the diffusion velocity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
