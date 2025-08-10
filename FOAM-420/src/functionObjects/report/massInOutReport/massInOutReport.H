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
    (c) 2018 Esi Ltd.

Class
    Foam::functionObjects::massInOutReport

Group
    grpReportFunctionObjects

Description
    Prints the surface area and mass/volume flux of inlets and outlet faces of
    a fazeZone

Usage
    <FO>
    {
        type            massInOutReport;
        functionObjectLibs ( "libreport.so" );
        faceZoneName    porous;
        phiName         phi;
        logToFile       true;
    }

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    massInOutReport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massInOutReport_H
#define functionObjects_massInOutReport_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "convergenceTermination/convergenceTermination.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class massInOutReport Declaration
\*---------------------------------------------------------------------------*/

class massInOutReport
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        // Read from dictonary

            //- faceZone name
            word faceZoneName_;

            //- flux name
            word phiName_;

            //- compressible
            bool compressible_;

            //- MassFlow in
            scalar massIn_;
            //- MassFlow out
            scalar massOut_;

            //- SurfaceArea in
            scalar surfaceIn_;
            //- SurfaceArea out
            scalar surfaceOut_;


    // Private Member Functions
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct and assignment
        massInOutReport(const massInOutReport&);
        void operator=(const massInOutReport&);


public:

    //- Runtime type information
    TypeName("massInOutReport");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        massInOutReport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~massInOutReport(){};


    // Member Functions

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
