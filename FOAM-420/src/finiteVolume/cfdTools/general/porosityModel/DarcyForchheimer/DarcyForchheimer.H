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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

Class
    Foam::porosityModels::DarcyForchheimer

Description
    Darcy-Forchheimer law porosity model
    Two modes:
    -coeffsMode == standard

            \f[
                S = - (\mu d + \frac{\rho |U|}{2} f) U
            \f]

        where
        \vartable
            d        | Darcy coefficient [1/m2]
            f        | Forchheimer coefficient [1/m]
        \endvartable

    -coeffsMode == pure

            \f[
                S = - (alpha d + beta |U|) U
            \f]

    The coefficients can depend with temperature using the
    temperatureDependence switch.
    -standard mode:
        alpha(t) or d(t) = ddn(t)*(alpha or d)
        beta(t) or f(t) = ffn(t)* (beta or f)

        ddn and ffn are Function1 tables

    -pure mode:
        alpha(t) or d(t) = (a1 + a2*Tref + a3*Tref^2)(Tref / T) * (alpha or d)
        beta(t) or f(t) =  (b1 + b2*Tref + b3*Tref^2)(Tref / T)^2 * (beta or f)

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The orientation of the porous region is defined with the same notation as
    a co-ordinate system, but only a Cartesian co-ordinate system is valid.

SourceFiles
    DarcyForchheimer.C
    DarcyForchheimerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef DarcyForchheimer_H
#define DarcyForchheimer_H

#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "dimensionedTypes/dimensionedTensor/dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class DarcyForchheimer Declaration
\*---------------------------------------------------------------------------*/

class DarcyForchheimer
:
    public porosityModel
{

    // Private data

        //- Darcy coeffient XYZ components (user-supplied) [1/m2]
        dimensionedVector dXYZ_;

        //- Forchheimer coeffient XYZ components (user-supplied) [1/m]
        dimensionedVector fXYZ_;

        //- Darcy coefficient - converted from dXYZ [1/m2]
        mutable List<tensorField> D_;

        //- Forchheimer coefficient - converted from fXYZ [1/m]
        mutable List<tensorField> F_;

        //- Darcy coefficient correction D/Dn
        mutable List<scalarField> Dcorr_;

        //- Forchheimer coefficient correction F/Fn
        mutable List<scalarField> Fcorr_;

        //- Name of density field
        word rhoName_;

        //- Name of dynamic viscosity field
        word muName_;

        //- Name of kinematic viscosity field
        word nuName_;

        //- Modification mode for keeping positive tensor coefficients
        Switch defaultMode_;

        //- coefficient mode (standard OF or pure AU+BU^2)
        const word coeffsMode_;

        //- Method of computing the alpha and beta in pure
        //  userDefined/cylindrical
        const word pureCoeffsDef_;

        //- alpha and beta is now  (rho*alpha, rho*beta)
        const Switch pureCoeffsRhoMultiplier_;

        //- Temperature dependence activation
        const Switch temperatureDependence_;

        //- Temperature coeffs mode (table - polynomial)
        const word temperatureCoeffsMode_;

        //- Is it spatially dependent porosity
        const Switch spatialDependancy_;

        // Name of temperature field
        word TName_;

        //- Temperature dependent variables

            // coeff temperature correction d/dn=f(T)
            autoPtr<Function1<scalar>> ddn_;

            // coeff temperature correction f/fn=f(T)
            autoPtr<Function1<scalar>> ffn_;

            scalar TRef_;

            scalar aT1_;
            scalar aT2_;
            scalar aT3_;

            scalar bT1_;
            scalar bT2_;
            scalar bT3_;

        //- pureCoeffsDef_ == cylindrical -- curvefit coeffs
            List<scalar> dpList_;
            List<scalar> qList_;

        //- currectTimeIndex
            mutable label curTimeIndex_;

            // Spatial variation of d coeff
            autoPtr<Function1<vector>> dSpatialFun_;

            // Spatial variation of f coeff
            autoPtr<Function1<vector>> fSpatialFun_;

        //- Override fieldIndex to deal with spatial dependancy
        label fieldIndex(const label i) const;

        // calculate temperature dependent DF multipliers
        scalar calcDcorrValue(const scalar&) const;
        scalar calcFcorrValue(const scalar&) const;

        //- Writes coefficient fields
        void writeCoeffsFields() const;

        // initialise
        void initialise();

        // compute cylindrical coeffs
        void computeCylindricalCoeffs();

        // f = ax + bx^2 fitting
        void polynomialFit
        (
            scalar&,
            scalar&,
            const scalar&,
            const scalar&,
            const scalar&,
            const scalar&,
            const scalar&
        ) const;

        // calculate correction for f and d coeffs
        void updateCorrCoeffs() const;

        // set coeffs
        void setTemperatureDependenceCoeffs();

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            vectorField&,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Apply
        template<class RhoFieldType>
        void adjointApply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& Uprimal,
            const vectorField& U
        ) const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        DarcyForchheimer(const DarcyForchheimer&);

        //- Disallow default bitwise assignment
        void operator=(const DarcyForchheimer&);


public:

    //- Runtime type information
    TypeName("DarcyForchheimer");

    //- Constructor
    DarcyForchheimer
    (
        const word& name,
        const word& modelType,
        const objectRegistry& obr,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~DarcyForchheimer();


    // Member Functions

        //- Transform the model data wrt mesh changes
        void updateSpatialCoeffs() const;

        //- Transform on construction
        virtual void calcTransformModelData();

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const;

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;
        virtual void correct(fvBlockMatrix<vector>& UEqn) const;

        //- Add resistance
        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        virtual void correct
        (
            fvBlockMatrix<vector>& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;


        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volVectorField& Uprimal
        ) const;


        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        ) const;


        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volVectorField&
        ) const;


        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&,
            const volVectorField&
        ) const
        {
            NotImplemented;
        }

    // I-O

        //- Write
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/porosityModel/DarcyForchheimer/DarcyForchheimerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
