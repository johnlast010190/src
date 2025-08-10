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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::incompressibleThreePhaseMixture

Description

SourceFiles
    incompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleThreePhaseMixture_H
#define incompressibleThreePhaseMixture_H

#include "transportModel/transportModel.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "incompressible/physicalProperties/physicalProperties.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class incompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleThreePhaseMixture
:
    public IOdictionary,
    public transportModel
{
    // Private data

        word phase1Name_;
        word phase2Name_;
        word phase3Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;
        volScalarField alpha3_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField nu_;

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
        autoPtr<viscosityModel> nuModel3_;

        autoPtr<physicalProperties> physProp1_;
        autoPtr<physicalProperties> physProp2_;
        autoPtr<physicalProperties> physProp3_;

        //dimensionedScalar rho1_;
        //dimensionedScalar rho2_;
        //dimensionedScalar rho3_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    // Constructors

        //- Construct from components
        incompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~incompressibleThreePhaseMixture()
    {}


    // Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        const word phase3Name() const
        {
            return phase3Name_;
        }

        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        volScalarField& alpha1()
        {
            return alpha1_;
        }

        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        volScalarField& alpha2()
        {
            return alpha2_;
        }

        const volScalarField& alpha3() const
        {
            return alpha3_;
        }

        volScalarField& alpha3()
        {
            return alpha3_;
        }

        tmp<volScalarField> alpha1lim() const
        {
            return min(max(alpha1_, scalar(0)), scalar(1));
        }
        tmp<volScalarField> alpha2lim() const
        {
            return min(max(alpha2_, scalar(0)), scalar(1));
        }
        tmp<volScalarField> alpha3lim() const
        {
            return min(max(alpha3_, scalar(0)), scalar(1));
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return physProp1_->rho();
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return physProp2_->rho();
        };

        //- Return const-access to phase2 density
        const dimensionedScalar& rho3() const
        {
            return physProp3_->rho();
        };

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase3 viscosityModel
        const viscosityModel& nuModel3() const
        {
            return nuModel3_();
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- heat capacity
        virtual tmp<volScalarField> Cp() const;

        //- thermal conductivity
        virtual tmp<volScalarField> lambda() const;

        //- turbulent Prandtl number
        virtual tmp<volScalarField> Prt() const;

        //- laminar Prandtl number
        virtual tmp<volScalarField> Pr() const;

        //- thermal diffusivity
        virtual tmp<volScalarField> alphaLam() const;

        //- Correct the laminar viscosity
        void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
