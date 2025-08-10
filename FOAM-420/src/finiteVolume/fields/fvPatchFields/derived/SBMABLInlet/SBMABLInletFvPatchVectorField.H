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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    SBMABLInletFvPatchVectorField


Description
    This boundary condition specifies a turbulent velocity inlet profile
    appropriate for LES/DES atmospheric boundary layers (ABL).  The mean
    velocity profile can be specified in 3 ways:

    1. Richards and Hoxey (1993)

        Mean velocity is derived from the friction velocity, flow direction and
        the direction of the parabolic co-ordinate \c z.

            \f[
                UMean = \frac{U^*}{K} ln\left(\frac{z + z_0}{z_0}\right)
            \f]

        where
        \vartable
            U^*     | frictional velocity
            K       | Karman's constant
            z       | distance from wall in vertical co-ordinate [m]
            z_0     | surface roughness length [m]
        \endvartable

        and:

            \f[
                U^* = K \frac{U_{ref}}{ln\left(\frac{Z_{ref} + z_0}{z_0}\right)}
            \f]

        where:
        \vartable
            U_{ref} | reference velocity at \f$Z_{ref}\f$ [m/s]
            Z_{ref} | reference height [m]
        \endvartable

        Assuming constant shear stress with height the turbulent kinetic energy
        can be found from:

            \f[
                k = \frac{(U^*)^2}{sqrt(Cmu)}
            \f]

        which is equivalent to a fixed intensity.
        The atmospheric BL dissipation profile, epsilon, is given by

            \f[
                \epsilon = \frac{(U^*)^3}{K(z + z_0)}
            \f]

        where
        \vartable
            U^*     | frictional velocity
            K       | Karman's constant
            z       | distance from ground in vertical co-ordinate [m]
            z_0     | surface roughness length [m]
        \endvartable


    2. AIJ (2004)

        Mean velocity is defined by a power law distribution

            \f[
                UMean = Us * \left(\frac{z}{zs}\right)^\alpha
            \f]


        where
        \vartable
            Us     | is the velocity at reference height
            z      | distance from wall in vertical co-ordinate [m]
            zs     | reference height [m]
            \alpha | power law exponent determined by terrain category
        \endvartable


        Turbulent intensity, k, is based on the estimation equation for
        vertical profile of intensity recommended by AIJ

            \f[
                i = \frac{\sigma(z)}{UMean} = 0.1*\frac{z}{z_G}^(-alpha-0.05)
            \f]
            \f[
                k = 3/2*sqr(i*|UMean|)
            \f]

        where z_G is the turbulent boundary layer height determined by terrain
        type

        Turbulent dissipation, \epsilon, is derived assuming local equilibrium
        of turbulent production and dissipation for the Power Law profile:

            \f[
                \epsilon = Cmu^0.5*k*frac{U_s}{z_s}\alpha
                *\left(frac{z}{z_s}\right)^(\alpha-1)
            \f]

        with Cmu = 0.09

    3. User defined tables for UMean and k or i(turbulent intensity) and
       epsilon and/or L(turbulent length scale). All table entries
       are a function of the local height above the ground level.

       If i is specified k is calculated from:

            \f[
                k = 3/2*sqr(i*|UMean|)
            \f]

        If epsilon and L is specified, k is calculated from:

            \f[
                k = (Cmu^-0.75*L * \epsilon)^(2/3)
            \f]

        If k and L is specified, epsilon is calculated from:

            \f[
                \epsilon =k*sqrt(k)/L *Cmu^0.75
            \f]


    Mean Reynolds stress for the turbulent perturbation prediction is derived
    from:

        \f[
            R = ((2.0/3.0)*I)*k - nut*twoSymm(fvc::grad(UMean))
        \f]


    nut the turbulent viscosity is derived from

        \f[
            nut = Cmu_*sqr(k)/epsilon
        \f]

    Unless explicitly defined, the turbulent length scale, L, is calculated
    from:

        \f[
            L =k*sqrt(k)/\epsilon *Cmu^0.75
        \f]



    Turbulent fluctuations are provided by the
    decayingTurbulenceInflowGenerator boundary. Currently only works for mean
    aligned with the X-axis.


    References:

        D.M. Hargreaves and N.G. Wright,  "On the use of the k-epsilon model
        in commercial CFD software to model the neutral atmospheric boundary
        layer", Journal of Wind Engineering and Industrial Aerodynamics
        95(2007), pp 355-369.

        Tominaga et.al., "AIJ guidelines for practical applications of CFD to
        pedestrian wind environment around buildings", J. of Wind Eng. and Ind.
        Aerodynamics, vol 96, pp 1749-1761, 2008

    \heading Patch usage

    \table
        Property     | Description             | Required   | Default value
        value        | Restart value       | yes               |

        fluctuations | Activate/deactivate turbulent content | no | true
        writeVortons | Output discrete vortex information, only needed on restart | no | false
        direction    | Direction of flow +1/-1 (x direction only | no | 1


        //////////////////////////////////////////////////////////////
        distanceType | Choose how to calculate z coordinate | no | wallDistance

        //distanceType == x|y|z
        origin | zero level for coordinate | no | 0.0

        //distanceType == vectorDistance
        axis | vertical direction | no | (0 0 1)
        origin | zero level for vector direction | no | (0 0 0)

        //distanceType == pointDistance
        origin | origin for distance calculation | no | (0 0 0)


        profileType | Input flow profile type id | yes

        PROFILE SPECIFIC ENTRIES
        profileType | AIJ
        Uref        | reference velocity at reference height | yes
        Href        | reference height | yes
        alpha       | power law velocity coefficient | yes
        zG          | boundary layer height | yes


        profileType | RichardsHoxey
        Uref        | reference velocity at reference height | yes
        Href        | reference height | yes
        z0          | roughness height | yes


        profileType | tabulated
        U | velocity profile as function of height | yes
        k | turbulence energy profile | no
        i | turbulence intensity profile | no
        epsilon | turbulence dissipation profile | no
        R | turbulent Reynolds stress profile | no
        L | turbulent length scale profile | no

        outOfBounds  | table bounds handling | no | clamp

    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch //mode AIJ
    {
        type            SBMABLInlet;
        value           uniform (11.3 0 0);
        delta           20;
        scaleVelocity   false; // if true, scales norm velocit
                               // to achieve constant flow rate

        distanceType    wallDistance;


        profileType     AIJ;
        AIJCoeffs
        {
            Uref            (11.3 0 0);
            Href            10.0;
            alpha           0.22;
            zG              20.0;
        }

    }
    myPatch //mode RichardsHoxey
    {
        type            SBMABLInlet;
        value           uniform (11.3 0 0);
        delta           20;

        distanceType    z;
        zCoeffs
        {
            origin (0 0 15);
        }

        profileType     RichardsHoxey;
        RichardsHoxeyCoeffs
        {
            Uref            (11.3 0 0);
            Href            10.0;
            z0              0.0;
        }
    }
    myPatch //userDefined
    {
        type            SBMABLInlet;
        value           uniform (11.3 0 0);
        delta           20;

        distanceType    vector;
        vectorCoeffs
        {
            axis       (0 0 1);
            origin     (0 0 0);
        }

        profileType     tabulated;
        tabulatedCoeffs
        {
            U table
            (
                (0 (0 0 0))
                (10 (11.3 0 0))
            );
            i table
            (
                (0 0.00)
                (1 0.1)
            );
            L table
            (
                (0 0.00)
                (200 100)
            );
            outOfBounds     clamp;
        }
    }
    \endverbatim

Note
    D.M. Hargreaves and N.G. Wright recommend Gamma epsilon in the
    k-epsilon model should be changed from 1.3 to 1.11 for consistency.
    The roughness height (Er) is given by Er = 20 z0 following the same
    reference.

SourceFiles
    SBMABLInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef SBMABLInletFvPatchVectorField_H
#define SBMABLInletFvPatchVectorField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/derived/inflowGenerator/decayingTurbulenceFvPatchVectorField.H"
#include "cfdTools/general/ABLProfile/ABLProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
       Class SBMABLInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class SBMABLInletFvPatchVectorField
:
    public decayingTurbulenceFvPatchVectorField
{

private:

    // Private data

        //- profile mode RichardsHoxey, AIJ, userDefined
        autoPtr<ABLProfile> profile_;

        //- turbulence structures on/off
        Switch fluctuations_;

        //- vorton writing
        //- writing vortons is necessary for smooth restarts
        Switch writeVortons_;

        // scale velocity to keep flow rate constant
        Switch scaleU_;

        bool updateTFR_;

        // flux name
        word phiName_;

        // density name
        word rhoName_;

        // target flow rate used for velocity scaling
        scalar targetFlowRate_;

    // Private member functions
private:

        //- scale provided flux to match flowRate
        void scaleFlux(const scalar, scalarField&) const;

        //- return target flow rate for inflow
        scalar calcFlowRate();

        //- scale normal velocity
        tmp<scalarField> scaleNormVelocity(const scalar) const;

public:

    //- Runtime type information
    TypeName("SBMABLInlet");


    // Constructors

        //- Construct from patch and internal field
        SBMABLInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SBMABLInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // SBMABLInletFvPatchVectorField onto a new patch
        SBMABLInletFvPatchVectorField
        (
            const SBMABLInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new SBMABLInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SBMABLInletFvPatchVectorField
        (
            const SBMABLInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new SBMABLInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
