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
    (c) 2015-2016 OpenCFD Ltd.

Class
    Foam::functionObjects::forces

Group
    grpForcesFunctionObjects

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches, and the resistance
    from porous zones.

    Forces and moments are calculated, with optional co-ordinate system and
    writing of binned data, where force and moment contributions are collected
    into a user-defined number of bins that span the input geometries for a
    user-defined direction vector.

    Data is written into multiple files in the
    postProcessing/\<functionObjectName\> directory:
    - force.dat          : forces in global Cartesian co-ordinate system
    - moment.dat         : moments in global Cartesian co-ordinate system
    - forceBin.dat       : force bins in global Cartesian co-ordinate system
    - momentBin.dat      : moment bins in global Cartesian co-ordinate system
    - localForce.dat     : forces in local co-ordinate system
    - localMoment.dat    : moments in local co-ordinate system
    - localForceBin.dat  : force bins in local co-ordinate system
    - localMomentBin.dat : moment bins in local co-ordinate system

Usage
    Example of function object specification:
    \verbatim
    forces1
    {
        type        forces;
        libs        ("libforces.so");
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
        referenceFrame rotatingFrame;
        definedInFrame false;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: forces       | yes         |
        log          | Write force data to standard output | no | no
        writeFields  | Write the force and moment fields | no | no
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        referenceFrame | Used for origin (see below) | no   | no
        definedInFrame | Bin direction defined in frame? (see below) | no  | yes
        CofR         | Centre of rotation (see below) | no   |
        porosity     | flag to include porosity contributions | no | no
        directForceDensity | Force density supplied directly (see below)|no|no
        fD           | Name of force density field (see below) | no | fD
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with incresing distance | yes |
    \endtable

Note
  - For incompressible cases, set \c rho to \c rhoInf.  You will then be
    required to provide a \c rhoInf value corresponding to the free-stream
    constant density.
  - If the force density is supplied directly, set the \c directForceDensity
    flag to 'yes', and supply the force density field using the \c
    fDName entry
  - The centre of rotation should be specified using referenceFrame,
    because in new system it is updatable (advantageous for moving meshes)
    It is used for moment calculations can either be
    specified by as an origion of local coordinate system. When frame isn't
    defined either CofR or coordinateSystem must be specified.
  - Frame additionally allows to specify direction of binning local to frame
    (definedInFrame true;).

    Examples:
    \verbatim
        referenceFrame        rotatingFrame;
        definedInFrame        false;
    \endverbatim
    or
    \verbatim
        CofR        (0 0 0);
    \endverbatim
    or
    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            e3      (0 0 1);
            e1      (1 0 0);
        }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forces_H
#define functionObjects_forces_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "coordinate/systems/coordinateSystem.H"
#include "referenceFrames/coordinateFrame.H"
#include "fields/volFields/volFieldsFwd.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/Tuple2/Tuple2.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "primitives/bools/Switch/Switch.H"
#include "meshes/polyMesh/zones/cellZone/cellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class forces
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;

        // File streams

            //- Forces
            autoPtr<OFstream> forceFilePtr_;

            //- Moments
            autoPtr<OFstream> momentFilePtr_;

            //- Force bins
            autoPtr<OFstream> forceBinFilePtr_;

            //- Moment bins
            autoPtr<OFstream> momentBinFilePtr_;


        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- Is the force density being supplied directly?
            Switch directForceDensity_;

            //- The name of the force density (fD) field
            word fDName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Coordinate system used when evaluting forces/moments
            coordinateSystem coordSys_;

            //- Reference Frame Ptr
            coordinateFrame* coorFramePtr_;

            //- Are variables defined in the frame of reference?
            Switch definedInFrame_;

            //- Flag to indicate whether we are using a local co-ordinate sys
            bool localSystem_;

            //- Flag to include porosity effects
            bool porosity_;

            //- Names of zones that will be taken into account (if present)
            wordList porosityZones_;


            // Bin information

                //- Number of bins
                label nBin_;

                //- Direction used to determine bin orientation
                vector binDir_;

                //- Distance between bin divisions
                scalar binDx_;

                //- Minimum bin bounds
                scalar binMin_;

                //- Bin positions along binDir
                List<point> binPoints_;

                //- Should bin data be cumulative?
                bool binCumulative_;


            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;

            //- Write forces in old single file format
            bool combineFiles_;

    // Protected Member Functions

        //- Local coordinate system
        const coordinateSystem& csys() const;

        //- Frame bin direction
        virtual vector binDir() const;

        //- Create a field name
        word fieldName(const word& name) const;

        //- Create the output files
        virtual void createFiles();

        //- Write header for old style forces.dat
        void writeIntegratedHeader();

        //- Write header for integrated data
        virtual void writeIntegratedHeader
        (
            const word& header,
            Ostream& os
        ) const;

        //- Write header for binned data (old format)
        void writeBinHeader();

        //- Write header for binned data
        virtual void writeBinHeader(const word& header, Ostream& os) const;

        //- Initialise the fields
        void initialise();

        //- Initialise the collection bins
        void initialiseBins();

        //- Reset the fields prior to accumulation of force/moments
        void resetFields();

        //- Reset the fields prior to accumulation of force/moments
        bool checkZone(const cellZone&);

        //- Return the effective viscous stress (laminar + turbulent).
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- return p
        virtual const tmp<volScalarField> getp();

        //- Accumulate bin data
        void applyBins
        (
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP,
            const vectorField& d
        );

        //- Add patch contributions to force and moment fields
        void addToFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Add cell contributions to force and moment fields
        void addToFields
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Helper function to write integrated forces and moments (old format)
        void writeIntegratedForceMoment();

        //- Helper function to write integrated forces and moments
        void writeIntegratedForceMoment
        (
            const string& descriptor,
            const vectorField& fm0,
            const vectorField& fm1,
            const vectorField& fm2,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write force data
        void writeForces();

        //- Helper function to write binned forces and moments (old format)
        void writeBinnedForceMoment();

        //- Helper function to write binned forces and moments
        void writeBinnedForceMoment
        (
            const List<Field<vector>>& fm,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write binned data
        void writeBins();

        //- Disallow default bitwise copy construct
        forces(const forces&);

        //- Disallow default bitwise assignment
        void operator=(const forces&);


public:

    //- Runtime type information
    TypeName("forces");


    // Constructors

        //- Construct from Time and dictionary
        forces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool baseClass = false
        );

        //- Construct from objectRegistry and dictionary
        forces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool baseClass = false
        );


    //- Destructor
    virtual ~forces();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Calculate the forces and moments
        virtual void calcForcesMoment();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();


    // Access Functions

        //- Reference pressure
        const scalar& pRef()
        {
            return pRef_;
        }

        //- Reference density
        const scalar& rhoRef()
        {
            return rhoRef_;
        }

        //- Name of pressure field
        const word& pName()
        {
            return pName_;
        }

        //- Name of velocity field
        const word& UName()
        {
            return UName_;
        }

        //- Name of density field
        const word& rhoName()
        {
            return rhoName_;
        }

        //- Patch set
        labelHashSet& patchSet()
        {
            return patchSet_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
