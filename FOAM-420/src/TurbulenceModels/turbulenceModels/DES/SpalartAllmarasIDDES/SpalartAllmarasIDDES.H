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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Class
    Foam::LESModels::SpalartAllmarasIDDES

Group
    grpDESIurbulence

Description
    SpalartAllmaras IDDES turbulence model for incompressible and compressible
    flows

    Reference:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2008).
        A hybrid RANS-LES approach with delayed-DES and wall-modelled LES
        capabilities.
        International Journal of Heat and Fluid Flow, 29(6), 1638-1649.
    \endverbatim

SourceFiles
    SpalartAllmarasIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasIDDES_H
#define SpalartAllmarasIDDES_H

#include "DES/SpalartAllmarasDES/SpalartAllmarasDES.H"
#include "LES/LESdeltas/IDDESDelta/IDDESDelta.H"
#include "LES/LESdeltas/IDDESDelta2/IDDESDelta2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class SpalartAllmarasIDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasIDDES
:
    public SpalartAllmarasDES<BasicTurbulenceModel>
{
    // Private Member Functions

        //- Check that the supplied delta is an IDDESDelta2
        const IDDESDelta2& setDelta() const;

        tmp<volScalarField> alpha() const;
        tmp<volScalarField> ft(const volScalarField& magGradU) const;
        tmp<volScalarField> fl(const volScalarField& magGradU) const;

        tmp<volScalarField> rd
        (
            const volScalarField& nur,
            const volScalarField& magGradU
        ) const;

        //update alpha and derived variables
        void correctAlphas();

        //- Delay function
        tmp<volScalarField> fdt(const volScalarField& magGradU) const;

        // Disallow default bitwise copy construct and assignment
        SpalartAllmarasIDDES(const SpalartAllmarasIDDES&);
        void operator=(const SpalartAllmarasIDDES&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cdt1_;
            dimensionedScalar Cdt2_;
            dimensionedScalar Cl_;
            dimensionedScalar Ct_;

        // Fields

            const IDDESDelta2& IDDESDelta_;

            // re-use fields; alpha-derived (includes transcendentals)
            autoPtr<volScalarField> fB_;
            autoPtr<volScalarField> fe1_;


    // Protected Member Functions

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasIDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasIDDES
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
    virtual ~SpalartAllmarasIDDES()
    {}


    // Member Functions

        //- Correct nuTilda related properties and alphas (if mesh is moving)
        virtual void correct();

        //- Re-read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
