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
    (c) 2022 Esi Ltd.

Class
    Foam::porousPressureFvPatchScalarField

Description
    Porous pressure BC. Supported models:

    - porousMode = darcyForchheimer (DF)

    The porous BC introduces a pressure jump defined by:

        \f[
            \Delta p = -(d \mu U + 0.5 f \rho |U|^2 )L
        \f]

    where

    \vartable
        p      | pressure [Pa]
        \rho   | density [kg/m3]
        \mu    | laminar viscosity [Pa s]
        d      | Darcy coefficient
        f      | inertial coefficient
        L      | porous media length in the flow direction
    \endvartable

    - porousMode = alphaBeta (AB)

    The porous BC introduces a pressure jump defined by:

        \f[
            \Delta p = -rho( \alpha U + \beta |U|^2 )L
        \f]

    where

    \vartable
        p      | pressure [Pa]
        \rho   | density [kg/m3]
        \alpha  | Darcy coefficient
        \beta   | inertial coefficient
        L      | porous media length in the flow direction
    \endvartable


    - porousMode = powerLaw (PL)

    The porous BC introduces a pressure jump defined by:

        \f[
            \Delta p = -\rho C_0 |U|^{(C_1 - 1)} U L
        \f]

    where

    \vartable
        p      | pressure [Pa]
        \rho   | density [kg/m3]
        \C_0   | model linear coefficient
        \C_1   | model exponent coefficient
        L      | porous media length in the flow direction
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        porousMode   | porous model used       | no          | darcyForchheimer
        U            | velocity field name     | no          | U
        rho          | density field name      | no          | rho
        d            | Darcy coefficient in DF | yes         |
        f            | inertial coefficient in DF | yes         |
        alpha        | Darcy coefficient in AB  | yes         |
        beta         | inertial coefficient in AB| yes         |
        C0           | Darcy coefficient in PL  | yes         |
        C1           | inertial coefficient in PL | yes         |
        length       | porous media length in the flow direction | yes |
        temperatureDependence | Coeff correction for temperature | no | false
        T | name of T -- if temperatureDependence is true | no | T
        a1 | a1 coeff -- if temperatureDependence is true | no | 0.0
        a2 | a2 coeff -- if temperatureDependence is true | no | 0.0
        a3 | a3 coeff -- if temperatureDependence is true | no | 0.0
        b1 | b1 coeff -- if temperatureDependence is true | no | 0.0
        b2 | b2 coeff -- if temperatureDependence is true | no | 0.0
        b3 | b3 coeff -- if temperatureDependence is true | no | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            porousPressure;
        d               0.001;
        f               1000000;
        length          0.1;
        value           uniform 0;
    }
    \endverbatim


See also
    Foam::fixedValueFvPatchField
    Foam::coupledVelocityFvPatchVectorField
    Foam::basePressureFvPatchScalarField

SourceFiles
    porousPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef porousPressureFvPatchScalarField_H
#define porousPressureFvPatchScalarField_H

#include "derivedFvPatchFields/basePressure/basePressureFvPatchScalarField.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class porousPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class porousPressureFvPatchScalarField
:
    public basePressureFvPatchScalarField
{

public:

    // Public data

        //- porous jump mode

        enum porousMode
        {
            darcyForchheimer,
            alphaBeta,
            powerLaw
        };

        static const NamedEnum<porousMode, 3> porousModeNames_;

private:

    // Private data

        //- porous method
        porousMode porousMode_;

        //- Darcy pressure loss coefficient
        autoPtr<Function1<scalar>> d_;

        //- Inertia pressure lost coefficient
        autoPtr<Function1<scalar>> f_;

        //- Tangential component of inertia pressure lost coefficient
        autoPtr<Function1<scalar>> tangentf_;

        //- Alpha coefficient
        autoPtr<Function1<scalar>> alpha_;

        //- Beta coefficient
        autoPtr<Function1<scalar>> beta_;

        //- Tangential component of Beta coefficient
        autoPtr<Function1<scalar>> tangentBeta_;

        //- C0 coefficient
        autoPtr<Function1<scalar>> C0_;

        //- C1 coefficient
        autoPtr<Function1<scalar>> C1_;

        //- Porous media length
        const scalar length_;

        //- Relax
        const scalar relax_;


        // Variables related to the treatment of density on the boundary
        // In most cases density has to be relaxed

            //- Use average rho
            const Switch useAveragedRho_;

            //- Old field patch rho
            mutable scalarField oldRho_;

            //- Was old rho initialised?
            mutable bool isOldRhoInit_;

            //- Current relaxed field rho
            mutable scalarField rho_;

            //- Relaxation factor for rho
            const scalar rhoRelax_;


        //- name of the matrix diagonal inverse field
        word rAUName_;

        //- Temperature dependence activation
        Switch temperatureDependence_;

        //- Temperature dependent variables

            // Name of temperature field
            word TName_;

            //- Polynomial coefficient correction due to temperature
            // (aT1 + aT2*Tref + aT3*Tref^2)(Tref/T)
            // (bT1 + bT2*Tref + bT3*Tref^2)(Tref/T)^2
            scalar TRef_;
            scalar aT1_;
            scalar aT2_;
            scalar aT3_;
            scalar bT1_;
            scalar bT2_;
            scalar bT3_;

    // Private Functions


        // UpdateValue for segregated
        void updateSegBasedCoeffs();

        // calculate temperature dependent DF multipliers
        scalar calcLinearValue(const scalar&) const;
        scalar calcNonLinearValue(const scalar&) const;

        //- correct coefficients due to temperature dependency
        void correctTemperatureDepLinearTerm(scalarField&) const;
        void correctTemperatureDepNonLinearTerm(scalarField&) const;

        //- Multiply by rho in compressible formulation
        void multiplyByRho(scalarField& coeff) const;


public:

    //- Runtime type information
    TypeName("porousPressure");


    // Constructors

        //- Construct from patch and internal field
        porousPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porousPressureFvPatchScalarField
        //  onto a new patch
        porousPressureFvPatchScalarField
        (
            const porousPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        porousPressureFvPatchScalarField
        (
            const porousPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new porousPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        porousPressureFvPatchScalarField
        (
            const porousPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new porousPressureFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        virtual tmp<scalarField> C0Field() const;

        virtual tmp<scalarField> C1Field() const;

        virtual tmp<scalarField> C2Field() const;

        virtual tmp<scalarField> C3Field() const;

        virtual void updateCoeffs();

        //-
        virtual void boundaryRelaxMatrix(fvBlockMatrix<vector>&) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
