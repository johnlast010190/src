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
    (c) 2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::functionObjects::propellerInfo

Group
    grpForcesFunctionObjects

Description
    Calculates propeller performance and wake field properties.

    Controlled by executeControl:
    - Propeller performance
      - Thrust coefficient, Kt
      - Torque coefficient, 10*Kq
      - Advance coefficient, J
      - Open water efficiency, etaO
      - Written to postProcessing/<name>/<time>/propellerPerformance.dat

    Controlled by writeControl:
    - Wake field text file
      - Wake: 1 - UzMean/URef
      - Velocity in cylindrical coordinates at xyz locations
      - Written to postProcessing/<name>/<time>/wake.dat
    - Axial wake field text file
      - 1 - Uz/URef at r/R and angle
      - Written to postProcessing/<name>/<time>/axialWake.dat
    - Velocity surface
      - Written to postProcessing/<name>/surfaces/time>/disk.<fileType>

Usage
    Example of function object specification:

    \verbatim
    propellerInfo1
    {
        type            propellerInfo;
        libs            ("libforces.so");
        writeControl    writeTime;

        patches         ("propeller.*");

        // Setting velocity reference back to wake average Uz velocity
        // computed by this function object
        URef            functionObjectValue;
        functionObjectResult UzMean;

        //// Function1 type = 'constant'
        //// URef            constant 5;

        rho             rhoInf;    // incompressible
        rhoInf          1.2;

        // Propeller data:

        // Optionally write propeller performance data
        writePropellerPerformance yes;

        radius          0.1;

        rotationMode    specified;    // specified | MRF

        // rotationMode = specified:
        CofR            (0 -0.1 0);
        n               25.15;
        axis            (0 1 0);

        // Optional reference direction for angle (alpha) = 0
        alphaAxis       (1 0 0);

        //// rotationMode = mrf:
        //// MRF             MRFNames;
        //// (CofR, n and axis retrieved from MRF model)

        // Optionally write wake text files
        // Note: controlled by writeControl
        writeWakeFields yes;

        // Sample plane (disk) properties
        // Note: controlled by writeControl
        sampleDisk
        {
            r1                   0.05;
            r2                   0.2;
            nTheta               36;
            nRadial              10;
            interpolationScheme  cellPoint;
            errorOnPointNotFound false;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property        | Description                   | Required | Default value
        type            | Type name: propellerInfo      | yes      |
        log             | Write to standard output      | no       | no
        patches         | Patches included in the forces calculation | yes |
        p               | Pressure field name           | no       | p
        U               | Velocity field name           | no       | U
        rho             | Density field name            | no       | rho
        URef            | Reference velocity            | yes      |
        rotationMode    | Rotation mode (see below)     | yes      |
        CofR            | Sample disk centre            | no*      |
        n               | Revolutions per second        | no*      |
        axis            | Propeller axis                | no*      |
        alphaAxis       | Axis that defines alpha=0 dir | no       |
        MRF             | Name of MRF zone              | no*      |
        originOffset    | Origin offset for MRF mode    | no       | (0 0 0)
        writePropellerPerformance| Write propeller performance text file | yes |
        writeWakeFields | Write wake field text files   | yes      |
        r1              | Sample disk inner radius      | no       | 0
        r2              | Sample disk outer radius      | no*      |
        nTheta          | Divisions in theta direction  | no*      |
        nRadial         | Divisions in radial direction | no*      |
        interpolationScheme | Sampling interpolation scheme | no* | cell
    \endtable

Notes:

- URef is either a scalar Function1 type or a scalar that is set as equal
  to some quantity computed by this function object
- rotationMode is used to set the center of rotation (CofR), axis
  and revolutions per second
  - if set to 'specified' all 3 entries are required
    - note: CofR is the sample disk origin
- if writePropellerPerformance is set to on|true:
  - propellerPerformance text file will be written
- if writeWakeFields is set to on|true:
  - wake and axialWake text files will be written

See also
    Foam::functionObject::forces
    Foam::Function1

SourceFiles
    propellerInfo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_propellerInfo_H
#define functionObjects_propellerInfo_H

#include "forces/forces.H"
#include "meshes/meshShapes/face/faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Function1;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class propellerInfo Declaration
\*---------------------------------------------------------------------------*/

class propellerInfo
:
    public forces
{

public:

    enum class rotationMode
    {
        SPECIFIED
    };

    static const Enum<rotationMode> rotationModeNames_;


protected:

    // Protected data

        //- Propeller radius
        scalar radius_;

        //- Pointer to reference velocity
        autoPtr<Function1<scalar>> URefPtr_;

        //- Reference velocity function type
        const word foType_;

        //- Function object result name
        word foResultName_;

        //- Default value when function object result is unavailable
        scalar defaultValue_;

        //- Flag to indicate that default value is available
        bool hasDefaultValue_;

        //- Rotation mode
        rotationMode rotationMode_;

        //- Name of MRF zone (if applicable)
        word MRFName_;

        //- Propeller speed (revolutions per second)
        scalar n_;

        //- Flag to write performance data
        bool writePropellerPerformance_;

        //- Propeller performance file
        autoPtr<OFstream> propellerPerformanceFilePtr_;

        //- Flag to write wake fields
        bool writeWakeFields_;


        // Wake field output

            //- Number of surface divisions in theta direction
            label nTheta_;

            //- Number of surface divisions in radial direction
            label nRadial_;

            //- Surface points
            pointField points_;

            //- Flag to raise an error if the sample point is not found
            //- in the mesh. Default = false to enable, e.g.
            //- reduced geometry/symmetric cases
            bool errorOnPointNotFound_;

            //- Surface faces
            faceList faces_;

            //- Surface point sample cell IDs
            labelList cellIds_;

            //- List of participating points (parallel reduced)
            boolList pointMask_;

            //- Interpolation scheme
            word interpolationScheme_;

            //- Wake field file
            autoPtr<OFstream> wakeFilePtr_;

            //- Axial wake field file
            autoPtr<OFstream> axialWakeFilePtr_;

            //- Default value when a sample point is not found;
            //- default = scalar::min
            scalar nanValue_;


    // Protected Member Functions

        //- Set the coordinate system
        void setCoordinateSystem(const dictionary& dict);

        //- Set the rotational speed
        void setRotationalSpeed();

        //- Create output files
        virtual void createFiles() override;

        //- Return the velocity field
        const volVectorField& U() const;

        //- Get the reference velocity value
        scalar getURefValue() const;


        // Sample surface functions

            //- Set the faces and points for the sample surface
            void setSampleDiskGeometry
            (
                const coordinateSystem& coordSys,
                const scalar r1,
                const scalar r2,
                const scalar nTheta,
                const label nRadius,
                faceList& faces,
                pointField& points
            ) const;

            //- Set the sample surface based on dictionary settings
            void setSampleDiskSurface(const dictionary& dict);

            //- Set the sample cells corresponding to the sample points
            void updateSampleDiskCells();

            //- Return the area average of a field
            scalar meanSampleDiskField(const scalarField& field) const;

            //- Interpolate from the mesh onto the sample surface
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>& psi,
                const Type& defaultValue
            ) const;


        // Propeller performance text file

            //- Write the wake fields
            void writePropellerPerformance();


        // Wake text files

            //- Write the wake text file
            void writeWake(const vectorField& U, const scalar URef);

            //- Write the axial wake text file
            void writeAxialWake(const vectorField& U, const scalar URef);

            //- Write the wake fields
            void writeWakeFields(const scalar URef);


        //- No copy construct
        propellerInfo(const propellerInfo&) = delete;

        //- No copy assignment
        void operator=(const propellerInfo&) = delete;


public:

    //- Runtime type information
    TypeName("propellerInfo");


    // Constructors

        //- Construct from Time and dictionary
        propellerInfo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        propellerInfo
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~propellerInfo() = default;


    // Member Functions

        //- Read the propellerInfo data
        virtual bool read(const dictionary&) override;

        //- Execute main functions
        virtual bool execute() override;

        //- Write the output
        virtual bool write() override;

        void UpdateMesh(const mapPolyMesh& mpm);

        virtual void movePoints(const polyMesh& mesh) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
