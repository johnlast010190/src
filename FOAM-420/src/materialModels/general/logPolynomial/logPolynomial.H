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
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd.

Class
    Foam::logPolynomial

Description
    Transport package using polynomial functions of \c ln(T):

        \f[
            ln(field)    = \sum_{i=1}^N \left( a[i] * ln(T)^{i-1} \right)
        \f]

Usage

    \table
        Property        | Description
        <fieldName>ModelCoeffs | Polynomial coefficients
    \endtable

    Example of the specification:
    \verbatim
    <fieldName>ModelCoeffs
    {
        <fieldName>Coeffs      8 ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expressions are evaluated as so (default field T):

        \f[
            \field    = 1000 - 0.05 ln(T) + 0.003 ln(T)^2
        \f]

SourceFiles
    logPolynomialTransportI.H
    logPolynomialTransport.C

\*---------------------------------------------------------------------------*/

#ifndef logPolynomial_H
#define logPolynomial_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/generalModel.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class logPolynomial Declaration
\*---------------------------------------------------------------------------*/

class logPolynomial
:
    public materialModel,
    public generalModel
{
    // Private data

        //- Dynamic field polynomial coefficients
        polynomialFunction coeffs_;

        //- Field name to be used
        word fieldName_;

        //- Pointer to the field offset by reference
        const refScalarField* field_;

        //- Is rho dependent on f(p)
        bool incompressible_;

        //- Is the equation of state is isochoric i.e. rho = const
        bool isochoric_;

        //- Is the model isotropic?
        bool isotropic_;


public:

    // Constructors

        //- Construct from dictionary
        logPolynomial
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<logPolynomial> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("logPolynomial");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return incompressible_;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            return isochoric_;
        }

        //- Is it isotropic model?
        virtual bool isotropic() const override
        {
            return isotropic_;
        }

        //- Whether the function is uniform
        virtual bool isConst() const override
        {
            return field_->isConst();
        }

        //- Definition of the function for single value
        virtual inline scalar general(const scalar val) const;

        //- Cell calculations
        virtual inline scalar generalCell(const label celli) const override;

        //- Internal field calculation
        virtual tmp<scalarField> generalInternal() const override;

        //- Calculate values for single patch
        virtual tmp<scalarField> generalPatch(const label patchi) const override;

        //- Calculation of whole volumetric field
        virtual tmp<volScalarField> generalGeometric() const override;

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override
        {}

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "logPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
