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
    Foam::functionObjects::highSpeedFlowOutletMonitor

Group
    grpReportFunctionObjects

Description
    Prints the effective area and proccentage of the outlet with backflow.

Usage
    \table
        Property     | Description             | Required    | Default value
        inletName    | name of inlet patch     | yes         | none
        outletName   | name of outlet patch    | yes         | none
        R            | gas constant            | no          | RR/molWeight
        gamma        | Cp/Cv                   | no          | Cp/Cv
        rhoInlet     | averaged inlet density  | no          | calculated
        rhoOutlet    | averaged outlet density | no          | calculated
        T0           | upstream stagnation temperature | no  | calculated
        pRef         | reference pressure      | no          | looked-up
    \endtable

    Example of the boundary condition specification:
    \verbatim
    flowMonitor
    {
        type            highSpeedFlowOutletMonitor;
        libs            ("libreport.so");
        writeControl    writeTime;
        inletName       inlet;
        outletName      outlet;
        R               287;
        gamma           1.4;
        rhoInlet        1.2;
        rhoOutlet       $rhoInlet;
        T0              250;
        pRef            100000;
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    highSpeedFlowOutletMonitor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_highSpeedFlowOutletMonitor_H
#define functionObjects_highSpeedFlowOutletMonitor_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class highSpeedFlowOutletMonitor Declaration
\*---------------------------------------------------------------------------*/

class highSpeedFlowOutletMonitor
:
    public fvMeshFunctionObject,
    public writeFile
{
// Private member functions

    //- Return patch ids
    labelList patchIDs(const wordList& patchNames);

    //- Return zone ids
    labelList zoneIDs(const wordList& zoneNames);

    //- Patch area sum
    scalar sumPatchAreas(const labelList& patchIDs);

    //- Face zone area sum
    scalar sumZoneAreas(const labelList& zoneIDs);

    //- When zones are available interpolate fields on faces
    void faceInterpolate();

    //- Clear tmp fields
    void clearTmpFields();


protected:

   // Protected data

        //- Inlet patch names
        wordList inletPatches_;

        //- Inlet patch IDs
        labelList inPatchIDs_;

        //- Inlet face zone names
        wordList inletFaceZones_;

        //- Inlet face zone IDs
        labelList inFaceZoneIDs_;

        //- Outlet patch names
        wordList outletPatches_;

        //- Outlet patch IDs
        labelList outPatchIDs_;

        //- Outlet face zone names
        wordList outletFaceZones_;

        //- Outlet face zone IDs
        labelList outFaceZoneIDs_;

        //- Sum in areas
        scalar sumInAreas_;

        //- Sum in areas
        scalar sumOutAreas_;

        //- Pressure
        const volScalarField& p_;

        //- Container for face interpolated pressure
        tmp<surfaceScalarField> pFaces_;

        //- Temperature
        const volScalarField& T_;

        //- Container for face interpolated temperature
        tmp<surfaceScalarField> TFaces_;

        //- Velocity
        const volVectorField& U_;

        //- Container for face interpolated velocity
        tmp<surfaceVectorField> UFaces_;

        //- Mass/volumetric flow rate
        const surfaceScalarField& phi_;

        //- Density reference
        const volScalarField& rho_;

        //- Container for face interpolated density
        tmp<surfaceScalarField> rhoFaces_;

        //- Reference to magSf
        const surfaceScalarField& magSf_;

        //- Pointer to material tables
        const materialTables& mat_;


public:

    //- Runtime type information
    TypeName("highSpeedFlowOutletMonitor");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        highSpeedFlowOutletMonitor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~highSpeedFlowOutletMonitor(){};


    // Member Functions

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);

        //- Area weighted average
        scalar areaAverage
        (
            const volScalarField& field,
            const tmp<surfaceScalarField>& surfaceField,
            const bool isInlet = true
        );

        //- Area weighted average
        scalar areaAverage
        (
            const volVectorField& field,
            const tmp<surfaceVectorField>& surfaceField,
            const bool isInlet = true
        );

        //- Calc phi
        scalar calcPhi
        (
            const surfaceScalarField& field,
            const labelList& patchIds,
            const labelList& zoneIds
        );

        //- Calculate gamma surface average from both inlet and outlet
        scalar calcGamma();

        // Calc T0
        scalar calcT0();

        //- Calculates procentage of backflow
        // from input face flow rates
        scalar calcBackFlowProc();

        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct and assignment
        highSpeedFlowOutletMonitor(const highSpeedFlowOutletMonitor&);
        void operator=(const highSpeedFlowOutletMonitor&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
