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
    (c) 2020 Rudolf Hellmuth - Vascular Flow Technologies
    (c) 2020 Esi Ltd.
    (c) 2016 OpenFOAM Fundation

Class
    Foam::windkesslPressureFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a windkesselPressure condition for
    incompressible laminar flow. It is calculated by the following ordinary
    differential equation (ODE):

        \f[
             \frac{dp}{dt} - \frac{p - p_0}{R_d\,C}
           =   R_p \frac{dQ}{dt}
             + \frac{1}{C} \left( 1 + \frac{R_p}{R_d} \right) Q
        \f]

    where
    \vartable
        p   | patch pressure [m^2 s^{-2}] or [kg m^{-1} s^{-2}]
        p_0 | far field pressure [m^2 s^{-2}] or [kg m^{-1} s^{-2}]
        Q   | volumetric flow rate across patch [m^3 s^{-1}]
        R_d | distal resistance [m^{-1} s^{-1}] or [kg m^{-4} s^{-1}]
        R_p | proximal resistance [m^{-1} s^{-1}] or [kg m^{-4} s^{-1}]
        C   | compliance [m s^2] or [kg^{-1} m^4 s^2]
    \endtable

    \heading Patch usage

    \table
        Property   | Description                        | Req'd? | Default
        p0         | far field pressure                 | no     | 0.0
        Rd         | distal resistance                  | yes    |
        Rp         | proximal resistance                | no     | 0.0
        C          | compliance                         | yes    |
        pOld       | pressure of previous time step     | no     | 0.0
        Qold       | flow rate of previous time step    | no     | 0.0
        alpha      | dynamic pressure weight            | no     | 1.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type       windkesselPressure;
        //p0        0.0;
        Rd         9.9e3;
        //Rp         0.1e3;
        //C           1.2e-4;
        //pOld        100.0;
        //Qold        0.01
        //alpha       5.0;
    }
    \endverbatim

Note
    - The characteristic time constant of the the ODE is given by
      \f[ \tau = R_d \, C \f].
    - The impedance is given by \f[ Z = R_p + \frac{R_d}{1 + j\,\omega \,\tau \f]
      , where \f[ j \f] is the imaginary number, and \f[ \omega \f] is the pulse
      frequency. \f[ \tau \f] is a cut-off number between high and low frequancy
      impedance values.
    - The ODE is implicitly integrated in time using backward differences. The
      resulting equation takes both flow rate and pressure of the previous time
      step -- \verbatim Qold \endverbatim and \verbatim pOld \endverbatim,
      respectively.
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
    - Backflow velocity profile is more realistic if the velocity BC is
      \verbatim pressureInletOutletVelocity \endverbatim.

SeeAlso
    Foam::flowResistancePressureFvPatchScalarField
    Foam::pressureInletOutletVelocity

SourceFiles
    windkesselPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef WindkesselPressureFvPatchScalarField_H
#define WindkesselPressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class WindkesselPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class WindkesselPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Far field (reference) pressure
        scalar p0_;

        //- Proximal flow resistance
        scalar Rp_;

        //- Distal flow resistance
        scalar Rd_;

        //- Compliance
        scalar C_;

        //- Flow rate of previous time step
        scalar Qold_;

        //- Pressure of previous time step
        scalar pOld_;

        //- A dummy variable that's updated at every run call. Qold_ takes the
        //  value stored in this variable at the first iteration of every time
        //  step.
        scalar Qdummy_;

        //- Current time index. This is used for checking when to update the
        //  variables Qold_ and pOld_ in an if-statement.
        label curTimeIndex_;

        //- Dynamic pressure weight
        scalar alpha_;

public:

    //- Runtime type information
    TypeName("windkesselPressure");


    // Constructors

        //- Construct from patch and internal field
        WindkesselPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        WindkesselPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        WindkesselPressureFvPatchScalarField
        (
            const WindkesselPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        WindkesselPressureFvPatchScalarField
        (
            const WindkesselPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new WindkesselPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        WindkesselPressureFvPatchScalarField
        (
            const WindkesselPressureFvPatchScalarField&,
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
                new WindkesselPressureFvPatchScalarField
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

            //- Return the proximal resistance
            scalar Rp() const
            {
                return Rp_;
            }

            //- Return reference to the proximal resistance
            scalar& Rp()
            {
                return Rp_;
            }

            //- Return the distal resistance
            scalar Rd() const
            {
                return Rd_;
            }

            //- Return reference to the distal resistance
            scalar& Rd()
            {
                return Rd_;
            }

            //- Return the compliance
            scalar C() const
            {
                return C_;
            }

            //- Return reference to the compliance
            scalar& C()
            {
                return C_;
            }

            //- Return the old flow rate
            scalar Qold() const
            {
                return Qold_;
            }

            //- Return reference to the compliance
            scalar& Qold()
            {
                return Qold_;
            }

            //- Return the patch pressure of the previous time step
            scalar pOld() const
            {
                return pOld_;
            }

            //- Return reference to the patch pressure of the previous time step
            scalar& pOld()
            {
                return pOld_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate the patch pressure
            virtual scalar pLumpedModel(scalar Q);

            //- Calculate the average of any variable over the patch
            scalar patchAvg(scalarField& var);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
