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
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::RASModels::EdwardsSA

Group
    grpRASTurbulence

Description
    Spalart-Allmaras one-eqn mixing-length model for incompressible and
    compressible external flows.

    Reference:
    \verbatim
        "Comparison of Eddy-Viscosity-Transport Turbulence Models for
        Three-Dimensional, Shock Separated Flowfields"
        J.R. Edwards,
        S. Chandra
        AAIA Journal, Vol 34, No. 4, 1996, pp. 756-763
    \endverbatim

SourceFiles
    EdwardsSA.C

\*---------------------------------------------------------------------------*/

#ifndef EdwardsSA_H
#define EdwardsSA_H

#include "RAS/RASModel/RASModel.H"
#include "RAS/SpalartAllmaras/SpalartAllmaras.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class EdwardsSA Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class EdwardsSA
:
    public SpalartAllmaras<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        EdwardsSA(const EdwardsSA&);
        EdwardsSA& operator=(const EdwardsSA&);


protected:

    // Protected data

        // Model coefficients


        // Fields



    // Protected Member Functions

        virtual tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        virtual tmp<volScalarField> fw(const volScalarField& Stilda) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("EdwardsSA");


    // Constructors

        //- Construct from components
        EdwardsSA
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~EdwardsSA()
    {}


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EdwardsSA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
