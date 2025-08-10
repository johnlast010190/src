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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::fv::panaRotorDiskSource

Description
    Cell based momemtum source

    Source approximates the mean effects of rotor forces on a cylindrical
    region within the domain

    \heading Source usage

    Example usage:
    \verbatim
    panaRotorDiskSourceCoeffs
    {
        fieldNames      (U);    // names of fields on which to apply source
        rhoName         rho;    // density field if compressible case
        nBlades         3;      // number of blades
        tipEffect       0.96;   // normalised radius above which lift = 0

        inletFlowType   local;  // inlet flow type specification

        geometryMode    auto;   // geometry specification

        refDirection    (-1 0 0); // reference direction
                                  // - used as reference for psi angle

        trimModel       fixed;  // fixed || targetForce

        flapCoeffs
        {
            beta0           0;  // coning angle [deg]
            beta1c          0;  // lateral flapping coeff (cos coeff)
            beta2s          0;  // longitudinal flapping coeff (sin coeff)
        }

        blade
        {
            ...
        }

        profiles
        {
            ...
        }
    }
    \endverbatim

    Where:
    Valid options for the \c geometryMode entry include:
    - auto          : determine rototor co-ord system from cells
    - specified     : specified co-ord system

    Valid options for the \c inletFlowType entry include:
    - fixed         : specified velocity
    - local         : use local flow conditions
    - surfaceNormal : specified normal velocity (positive towards rotor)


SourceFiles
    panaRotorDiskSource.C
    panaRotorDiskSourceTemplates.C

CAE/RGK 2017/05: for Pana, employ fan outlet flow rate Q
as opposed to normal velocity.

\*---------------------------------------------------------------------------*/

#ifndef panaRotorDiskSource_H
#define panaRotorDiskSource_H

//#include "cfdTools/general/fvOptions/fvOption.H"
#include "cellSetOption/cellSetOption.H"
#include "coordinate/systems/cylindricalCS.H"
#include "primitives/enums/NamedEnum.H"
#include "sources/derived/rotorDiskSource/bladeModel/bladeModel.H"
#include "sources/derived/rotorDiskSource/profileModel/profileModelList.H"
#include "volMesh/volMesh.H"
#include "dimensionSet/dimensionSet.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class trimModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class panaRotorDiskSource Declaration
\*---------------------------------------------------------------------------*/

class panaRotorDiskSource
:
    public cellSetOption
{
public:

    enum geometryModeType
    {
        gmAuto,
        gmSpecified
    };
    static const NamedEnum<geometryModeType, 2> geometryModeTypeNames_;

    enum inletFlowType
    {
        ifFixed,
        ifSurfaceNormal,
        ifLocal
    };
    static const NamedEnum<inletFlowType, 3> inletFlowTypeNames_;


protected:

    // Helper structures to encapsulate flap and trim data
    // Note: all input in degrees (converted to radians internally)

        struct flapData
        {
            scalar beta0;   // coning angle
            scalar beta1c;  // lateral flapping coeff (cos coeff)
            scalar beta2s;  // longitudinal flapping coeff (sin coeff)
        };


    // Protected data

        //- Name of density field
        word rhoName_;

        //- Reference density for rhoName = 'none'
        scalar rhoRef_;

        //- Rotational speed [rad/s]
        //  Positive anti-clockwise when looking along -ve lift direction
        scalar omega_;

        //- Number of blades
        label nBlades_;

        //- Inlet flow type
        inletFlowType inletFlow_;

        //- Flow Rate -  m^3/s,
        //flowRate0_ is the user-specified value for the fan
        //flowRate_ is the actual value computed at the fan
        scalar flowRate0_, flowRate_ ;

        // Fan Axis (same sense as flow rate)
        vector fanAxis_;

        //- Inlet velocity for specified iinflow
        vector inletVelocity_;

        // - When to apply feedback
        scalar  feedbackTime_;

        //- Tip effect [0-1]
        //  Ratio of blade radius beyond which lift=0
        scalar tipEffect_;

        //- Blade flap coefficients [rad/s]
        flapData flap_;

        //- Cell centre positions in local rotor frame
        //  (Cylindrical r, theta, z)
        List<point> x_;

        //- Rotation tensor for flap angle
        List<tensor> R_, Rcyl_;

        //- Inverse rotation tensor for flap angle
        List<tensor> invR_, invRcyl_;

        //- Area [m2] and volume of rotorDiskCell
        List<scalar> area_;

        // - Area Vectors
        List<vector>  areaSf_;

        // - rotorDisk Area
        scalar  rotorDiskA_;

        //- Rotor co-ordinate system (r, theta, z)
        cylindricalCS coordSys_;

        //- Max/min radii
        scalar rMax_;

        //- Volume of cells in rotor disk CAE/RGK 2017/09/12
        scalar  rotorDiskV_;

        //- Trim model
        autoPtr<trimModel> trim_;

        //- Blade data
        bladeModel blade_;

        //- Profile data
        profileModelList profiles_;


    // Protected Member Functions

        //- Check data
        void checkData();

        //- Set the face areas per cell, and optionally correct the rotor axis
        void setFaceArea(vector& axis, const bool correct);

        //- Create the co-ordinate system
        void createCoordinateSystem();

        //- Construct geometry
        void constructGeometry();

        //- Return the inlet flow field
        tmp<vectorField> inflowVelocity(const volVectorField& U) const;

        //- Helper function to write rotor values
        template<class Type>
        void writeField
        (
            const word& name,
            const List<Type>& values,
            const bool writeNow = false
        ) const;


public:

    //- Runtime type information
    TypeName("panaRotorDisk");


    // Constructors


        //- Construct from components
        panaRotorDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~panaRotorDiskSource();


    // Member Functions

        // Access

            //- Return the reference density for rhoName = 'none'
            inline scalar rhoRef() const;

            //- Return the rotational speed [rad/s]
            //  Positive anti-clockwise when looking along -ve lift direction
            inline scalar omega() const;

            //- Return the cell centre positions in local rotor frame
            // (Cylindrical r, theta, z)
            inline const List<point>& x() const;

            //- Return the rotor co-ordinate system (r, theta, z)
            inline const cylindricalCS& coordSys() const;

            //- Return true if solving a compressible case
            inline bool compressible() const;

            //- Return the density field [kg/m3]
            inline tmp<volScalarField> rho() const;

            //- Scalar CAE/RGK 2017/09/28 - this is called after the UEqn is solved, applies specific field values.
            virtual void correct(volVectorField& fld);


        // Evaluation

            //- Calculate forces
            void calculate
            (
                const vectorField& U,
                const scalarField& thetag,
                vectorField& force,
                const bool divideVolume = true,
                const bool output = true
            ) const;


        // Source term addition

            //- Source term to fvMatrix<vector>
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "panaRotorDiskSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "panaRotorDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
