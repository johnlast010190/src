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
    (c) 2016 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::porousBafflePressureAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclicAMI
    condition as a base. There are three different porous jump models
    implemented.

    - porousMode = darcyForchheimer (DF)

    The porous baffle introduces a pressure jump defined by:

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

    The porous baffle introduces a pressure jump defined by:

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

    The porous baffle introduces a pressure jump defined by:

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

    For DF and AB models, the darcy and inertia coeffs can be corrected for
    temperature dependency using polynomial expression:

    D|alpha *= (aT1 + aT2*Tref + aT3*Tref^2)(Tref/T)
    F|beta  *= (bT1 + bT2*Tref + bT3*Tref^2)(Tref/T)^2

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclicAMI| yes |
        porousMode   | porous model used       | no          | darcyForchheimer
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        d            | Darcy coefficient in DF | yes         |
        f            | inertial coefficient in DF | yes         |
        alpha        | Darcy coefficient in AB  | yes         |
        beta         | inertial coefficient in AB| yes         |
        C0           | Darcy coefficient in PL  | yes         |
        C1           | inertial coefficient in PL | yes         |
        length       | porous media length in the flow direction | yes |
        uniformJump  | applies a uniform pressure drop on the patch based on the
        relax        | relaxation on pressure update | no    | 1.0
        minJump      | Clip on minimum pressure | no    | -GREAT
        maxJump      | Clip on maximum pressure | no    | GREAT
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
        type            porousBafflePressureAMI;
        patchType       cyclicAMI;
        jump            uniform 0;
        d               0.001;
        f               1000000;
        length          0.1;
        uniformJump     false;
        value           uniform 0;
    }
    \endverbatim

Note
     The underlying \c patchType should be set to \c cyclicAMI

SourceFiles
    porousBafflePressureAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBafflePressureAMIFvPatchField_H
#define porousBafflePressureAMIFvPatchField_H

#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchField.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class porousBafflePressureAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class porousBafflePressureAMIFvPatchField
:
    public fixedJumpAMIFvPatchField<scalar>,
    public referenceFrameFvPatch<vector>
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

        //- Name of flux field (default = phi)
        const word phiName_;

        //- Name of density field (default = rho)
        const word rhoName_;

        //- Darcy pressure loss coefficient
        autoPtr<Function1<scalar>> d_;

        //- Inertia pressure lost coefficient
        autoPtr<Function1<scalar>> f_;

        //- parameter 1 for alternative model
        autoPtr<Function1<scalar>> Ct1_;

        //- parameter 2 for alternative model
        autoPtr<Function1<scalar>> Ct2_;

        //- Alpha coefficient
        autoPtr<Function1<scalar>> alpha_;

        //- Beta coefficient
        autoPtr<Function1<scalar>> beta_;

        //- C0 coefficient
        autoPtr<Function1<scalar>> C0_;

        //- C1 coefficient
        autoPtr<Function1<scalar>> C1_;

        //- Porous media length
        scalar length_;

        //- Relaxation factor for updating the pressure
        scalar relax_;

        //- Min jump limiter
        scalar minJump_;

        //- Max jump limiter
        scalar maxJump_;

        //- Aplies uniform pressure drop
        Switch uniformJump_;

        //- use alternative formulation
        Switch alternativeModel_;

        //- use alternative formulation
        Switch fluxNormalVelocity_;

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

        // calculate temperature dependent DF multipliers
        scalar calcDcorrValue(const scalar&) const;
        scalar calcFcorrValue(const scalar&) const;

        //- limiters and relaxation of jump
        void stabiliseJump(const scalarField&);

        //- compute normal velocity
        tmp<scalarField> computeUn() const;

        //- make velocity uniform if uniformJump is true
        void computeUniformVelocity(scalarField&) const;

        //- correct coefficients due to temperature dependency
        void correctTemperatureDependency(scalarField&, scalarField&) const;

        //- DF jump
        void computeDFJump();

        //- AB jump
        void computeABJump();

        //- PL jump
        void computePLJump();

        //- print debugging information
        void printDebugInfo(const scalarField&) const;

public:

    //- Runtime type information
    TypeName("porousBafflePressureAMI");


    // Constructors

        //- Construct from patch and internal field
        porousBafflePressureAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousBafflePressureAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porousBafflePressureAMIFvPatchField
        // onto a new patch
        porousBafflePressureAMIFvPatchField
        (
            const porousBafflePressureAMIFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        porousBafflePressureAMIFvPatchField
        (
            const porousBafflePressureAMIFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new porousBafflePressureAMIFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        porousBafflePressureAMIFvPatchField
        (
            const porousBafflePressureAMIFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new porousBafflePressureAMIFvPatchField(*this, iF)
            );
        }


    // Member functions


        // Evaluation functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
