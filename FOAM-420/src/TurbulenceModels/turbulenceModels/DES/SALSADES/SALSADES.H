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
    (c) 2014 Esi Ltd.

Class
    Foam::LESModels::SALSADES

Group
    grpDESIurbulence

Description
    SALSADES for compressible flows

    Implemented according to
    \verbatim
        "A restatement of the Spalart-Allmaras Eddy-Viscosity Model in Strain-
        Adaptive Formulation"
        Rung, T., Bunge, U., Schatz, M., and Thiele, F.,
        AIAA Journal, Vol 41, No 7, 2003, pp. 1396-1399.
    \endverbatim

    Disables ashfordCorrection and sets sigmaMu to 1

SourceFiles
    SALSADES.C

\*---------------------------------------------------------------------------*/

#ifndef SALSADES_H
#define SALSADES_H

#include "DES/SpalartAllmarasDES/SpalartAllmarasDES.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       Class SALSADES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SALSADES
:
    public SpalartAllmarasDES<BasicTurbulenceModel>
{
    // Private Data

        //free stream stagnation density - to be added to coeff dict by user
        // if not specified density effect on "rd" will be ignored
        autoPtr<dimensionedScalar> rho0_;

    // Private Member Functions


        // Disallow default bitwise copy construct and assignment
        SALSADES(const SALSADES&);
        SALSADES& operator=(const SALSADES&);

protected:

    // Protected Member Functions

        virtual tmp<volScalarField> Cb1
        (
            const volScalarField& SStar
        ) const;

        virtual tmp<volScalarField> Omega(const volTensorField& gradU) const;

        virtual tmp<volScalarField> STilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volScalarField& SStar,
            const volScalarField& dTilda
        ) const;

        virtual tmp<volScalarField> r
        (
            const volScalarField& visc,
            const volScalarField& Stilda,
            const volScalarField& dTilda
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("SALSADES");


    // Constructors

        //- Constructor from components
        SALSADES
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
    virtual ~SALSADES()
    {}


    // Member Functions

        //- Read LESProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "SALSADES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
