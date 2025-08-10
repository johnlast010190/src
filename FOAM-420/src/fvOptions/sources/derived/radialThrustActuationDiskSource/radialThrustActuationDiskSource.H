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
    Foam::fv::radialThrustActuationDiskSource

Description
    Actuation disk source including radial thrust. It enables the definition of
    constant thrust or thrust based on open-water curves through the identification
    of the self-propulsion point. PID controller is included to adjust disk
    revolutions and match target force or speed.

Usage
    Example usage:
    \verbatim
		actuatorDisk
		{
		    type                    radialThrustActuationDiskSource;
		    active                  true;
		    selectionMode           cellZone;
		    cellZone                cellZoneName;
            diskType                owCurves; // fixed
            vesselSpeed             ( -2.196 0 0 );
            diskDir                 ( -1 0 0 );
            T                       59;
            SFC                     30.5;
            fieldNames              ( U );
            PoverD                  1;
            R                       0.127;
            rh                      0.024;
            rhoName                 rho;
            rhoInf                  998.63;
            propPosition            ( 3.48 0 -0.21181 );
            propOrientation         ( 1 ( 0 0 0 ) );
            coeffs                  ( 1 1 0.5 );
            J                       ( 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 );
            KT                      ( 0.4942 0.4763 0.4553 0.4321 0.4072 0.3814 0.3551 0.3286 0.3023 0.2763 0.2509 0.2259 0.2013 0.177 0.1526 0.128 0.1025 0.0757 0.047 );
            KQ                      ( 0.6986 0.672 0.6438 0.6142 0.5838 0.5527 0.5211 0.4894 0.4576 0.4258 0.3941 0.3625 0.3309 0.2992 0.2673 0.235 0.2021 0.1682 0.1331 );
            referenceFrame          shipFrame;
            propulsionTarget        force;  // speed
            wakeFraction            0.25;
            forces
            {
                type            forces;
                patches         ( hull hullTop );
                rho             rho;
                CofR            ( 0.10848 0 -0.11142 );
            }
            controller
            {
                coeffsType      calculated;
                tuningMode      nuutinen;
                startTime       14;
                estimatedThrust 59;
                estimatedRPM    9.5;
            }

            // If present, sample velocity field on the tri-surface
            sampleDisk
            {
                type            surfaces;
                libs            ( "libsampling.so" );
                writeControl    writeTime;
                fields          ( U );
                interpolationScheme cell;
                surfaceFormat   vtk;
                surfaces
                (
                    Disk
                    {
                        type            triSurfaceToMesh;
                        surface         "propellerDisk.stl";
                        source          cells;
                        referenceFrame  shipFrame;
                        moveSurface     true;
                        interpolate     false;
                    }
                );
            }


		    //shipResistanceThrust    true;      // Default is false. If true, forces dict needed and T is recalculated based on forces on patches
			//forces
			//{
			//	type                  forces;
			//    patches             ( patches );
			//    rho          	      rho;
			//    CofR                (5.93526 0 -0.17658);
		    //}
		    //thrustCorrection        10;		   // Default is 0.0 (for full-scale ships). If defined thrust = forceCalculated - thrustCorrection
		}
    \endverbatim

SourceFiles
    radialThrustActuationDiskSource.C
    radialThrustActuationDiskSourceTemplates.C

References
    \verbatim
        Hoekstra, M. (2006). A RANS-based analysis tool for ducted propeller
        systems in open water condition. International shipbuilding progress,
        53(3), 205-227.
    \endverbatim

    \verbatim
        Birk, L. (2019). Fundamentals of ship hydrodynamics: Fluid mechanics,
        ship resistance and propulsion. John Wiley & Sons.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef radialThrustActuationDiskSource_H
#define radialThrustActuationDiskSource_H

#include "sources/derived/constantThrustActuationDiskSource/constantThrustActuationDiskSource.H"
#include "sources/derived/radialThrustActuationDiskSource/controllers/controllers.H"
#include "sampledSurface/sampledSurfaces/sampledSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class radialThrustActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class radialThrustActuationDiskSource
:
    public constantThrustActuationDiskSource
{
public:

    enum diskModeType
    {
        mdFixed,
        mdOWCurves
    };

    static const NamedEnum<diskModeType, 2> diskModeTypeNames_;

    enum propulsionModeType
    {
        mdSpeed,
        mdForce
    };

    static const NamedEnum<propulsionModeType, 2> propulsionModeTypeNames_;

protected:

    // Protected data

        //- Mode for the actuator disk
        diskModeType diskType_;

        //- Mode for the propulsion if a controller is defined (target speed or
        //  force)
        propulsionModeType propulsionType_;

private:

    // Private data

        //- Sample surfaces to calculate wake fraction
        mutable sampledSurfaces sampledSurfaces_;

        //- Controllers
        mutable controllers controller_;

        //- Pitch to diameter ratio
        scalar PoverD_;

        //- Disk Radius
        scalar R_;

        //- Hub Radius
        scalar rh_;

        //- Propeller position
        vector propPosition_;

        //- Propeller orientation
        quaternion propOrientation_;

        //- Coeffcients for the radial distribution
        FixedList<scalar, 3> radialCoeffs_;

        //- Activate thrust based on forces calculation
        bool shipResistanceThrust_;

        //- Thrust correction applied to thrust-based forces calculation (if shipResistanceThrust_ true)
        scalar thrustCorrection_;

        //- Thrust coefficients in open water conditions
        List<scalar> thrustCoeffs_;

        //- Torque coefficients in open water conditions
        List<scalar> torqueCoeffs_;

        //- Advance coefficients in open water conditions
        List<scalar> advanceCoeffs_;

        //- Wake fraction
        mutable scalar w_;

        //- Skin friction correction
        mutable scalar SFC_;

        //- Ship speed
        mutable vector uShip_;

        //- Time to reach the target ship speed
        scalar tr_;

        //- Reference density
        scalar rho_;

        //- Thrust coefficient at self-propulsion point
        mutable scalar Kts_;

        //- Propeller revolutions self-propulsion point
        mutable scalar n_;

        //- Is wake fraction updated?
        mutable bool updateWake_;

        //- Is wake field to be sampled?
        mutable bool sampleWake_;

        //- Disk velocity, if sampleWake true
        mutable vector uDisk_;

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addRadialActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        radialThrustActuationDiskSource(const radialThrustActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const radialThrustActuationDiskSource&);


public:

    //- Runtime type information
    TypeName("radialThrustActuationDiskSource");


    // Constructors

        //- Construct from components
        radialThrustActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~radialThrustActuationDiskSource()
    {}


    // Member Functions

        //- Return PoverD
        scalar PoverD() const
        {
            return PoverD_;
        }

        //- Disk area
        scalar R() const
        {
            return R_;
        }

        //- Disk area
        scalar rh() const
        {
            return rh_;
        }

        //- Propeller position
        const vector& propPosition() const
        {
            return propPosition_;
        }

        //- Propeller orientation
        const quaternion& propOrientation() const
        {
            return propOrientation_;
        }

        //- Return Skin Friction Correction
        const scalar& SFC() const
        {
            return SFC_;
        }

        //- Return averaged velocity on the sampling disk
        const vector& uDisk() const
        {
            return uDisk_;
        }

        //- Update propeller position
        bool updatePropellerPosition(const vector& propPosition);

        //- Update propeller orientation
        bool updatePropellerOrientation(const quaternion& propOrientation);

        //- Calculate self-propulsion point for a given ship-speed based on
        //  thrust identity
        bool calcSelfPropulsionPoint() const;

        //vector

        //- Sample velocity field on a tri-surface
        void sampleDisk() const;

        //- Compute wake fraction
        void wakeFraction() const;

        //- Apply controller to match force or speed
        void controlDisk() const;

        //- Source term to fvMatrix<vector>
        virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);

        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Write data
        virtual void writeData(Ostream&) const;

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "radialThrustActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
