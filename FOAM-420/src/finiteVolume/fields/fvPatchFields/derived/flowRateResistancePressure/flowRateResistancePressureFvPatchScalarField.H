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
    (c) 2016 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::flowRateResistancePressureFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a flowRateResistancePressure condition for
    incompressible laminar flow. It is calculated as:

        \f[
            p = p0 + R*Q
        \f]

    where
    \vartable
        p   | patch pressue [m^2/s^2]
        p0  | far field pressure [m^2/s^2]
        R   | single scalar value [m^{-1}s^{-1}]
        Q   | flow rate across patch [m^3/s]
    \endtable

    \heading Patch usage

    \table
        Property   | Description                        | Req'd? | Default
        p0         | far field pressure                 | no     | 0.0
        R          | resistance                         | yes    |
        alpha      | dynamic pressure weight            | no     | 1.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type       flowRateResistancePressure;
        //p0        0.0;
        R          1e3;
        //alpha       1.0;
    }
    \endverbatim

Note
    - This BC uses a calculates the total pressure in the backflow, for
      generating a more realistic backflow profile. This can be turned off by
      using \f[ \alpha = 0.0 \f]. The speed by which the backflow profile
      converges can be increased by using \f[ \alpha > 1.0 \f]. Note that the
      shape of the velocity profile can also be changed by \f[ \alpha \f].
    - Backflow velocity profile is more realistic if the velocity BC is
      \verbatim pressureInletOutletVelocity \endverbatim.
    - Under-relaxation factor of the (P/S)IMPLE algorithm improves convergence
      of the soft coupling between Neumann velocity BC and Dirichlet pressure
      BC. Low relaxation factors prevent overshooting.
    - The resistance can be calculated from the Hagen-Poiseuilli equation as
      \f[
           R = \frac{128}{\pi} \frac{\mu \, L}{D^4}
      \f]
      , where \f[\mu\f] is the dynamic viscosity of the fluid, and \f[L\f] and
      \f[D\f] are the length and diameter of the pipe, respectively.

SourceFiles
    flowRateResistancePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateResistancePressureFvPatchScalarField_H
#define flowRateResistancePressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class flowRateResistancePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class flowRateResistancePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Far field pressure
        scalar p0_;

        //- Flow resistance
        scalar R_;

        //- Dynamic pressure weight
        scalar alpha_;

public:

    //- Runtime type information
    TypeName("flowRateResistancePressure");


    // Constructors

        //- Construct from patch and internal field
        flowRateResistancePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateResistancePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        flowRateResistancePressureFvPatchScalarField
        (
            const flowRateResistancePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateResistancePressureFvPatchScalarField
        (
            const flowRateResistancePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new flowRateResistancePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateResistancePressureFvPatchScalarField
        (
            const flowRateResistancePressureFvPatchScalarField&,
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
                new flowRateResistancePressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions


        // Access

            //- Return the far field pressure
            scalar p0() const
            {
                return p0_;
            }

            //- Return reference to the far field pressure to allow adjustment
            scalar& p0()
            {
                return p0_;
            }

            //- Return the downstream resistance constant
            scalar R() const
            {
                return R_;
            }

            //- Return reference to the downstream resistance constant
            scalar& R()
            {
                return R_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate the patch pressure
            virtual scalar pLumpedModel(scalar Q);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
