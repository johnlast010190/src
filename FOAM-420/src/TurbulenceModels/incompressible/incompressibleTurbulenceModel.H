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
    (c) 2016 Esi Ltd.
    (c) 2013-2015 OpenFOAM Foundation

Class
    Foam::incompressibleTurbulenceModel

Description
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    incompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTurbulenceModel_H
#define incompressibleTurbulenceModel_H

#include "turbulenceModel.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                Class incompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTurbulenceModel
:
    public turbulenceModel
{

protected:

    // Protected data

        geometricOneField rho_;

        // Fields
        mutable volScalarField* alphatPtr_;

    // Protected member functions

        //- ***HGW Temporary function to be removed when the run-time selectable
        //  thermal transport layer is complete
        virtual void correctNut()
        {}


private:

    // Private Member Functions

        //- helper function for alphat construction
        tmp<volScalarField> alphatReadIfPresent
        (
            const word& fieldName,
            const fvMesh& mesh
        ) const;

        //- helper function for alphat creation
        tmp<volScalarField> alphatAutoCreate
        (
            const word& fieldName,
            const fvMesh& mesh
        ) const;

        //- Disallow default bitwise copy construct
        incompressibleTurbulenceModel(const incompressibleTurbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const incompressibleTurbulenceModel&);


public:

    //- Runtime type information
    TypeName("incompressibleTurbulenceModel");


    // Constructors

        //- Construct from components
        incompressibleTurbulenceModel
        (
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const word& propertiesName
        );


    //- Destructor
    virtual ~incompressibleTurbulenceModel();


    // Member Functions

        //- Return the density
        virtual tmp<volScalarField> rho() const = 0;

        //- Return the heat capacity
        virtual tmp<volScalarField> Cp() const = 0;

        //- Return the conductivity
        virtual tmp<volScalarField> lambda() const = 0;

        //- Return the turbulent Prandtl number
        virtual tmp<volScalarField> Prt() const = 0;

        //- Return the turbulent Prandtl number
        virtual tmp<volScalarField> Pr() const = 0;

        //- Return the laminar thermal diffusivity
        virtual tmp<volScalarField> alphaLam() const = 0;

        //- Return the turbulent thermal diffusivity
        virtual const volScalarField& alphat() const;

        //- Return the effective thermal diffusivity
        virtual const tmp<volScalarField> alphaEff() const;

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const;

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const;

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const;

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const;

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;

        virtual void correct();

        //- wall coefficient for block solvers
        virtual tmp<scalarField> uWallCoeffs
        (
            const label& patchI
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
