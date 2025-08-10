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
    Foam::functionObjects::surfaceAreaReport

Group
    grpReportFunctionObjects

Description
    Prints the surface area of a patch

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    surfaceAreaReport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceAreaReport_H
#define functionObjects_surfaceAreaReport_H

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
                           Class surfaceAreaReport Declaration
\*---------------------------------------------------------------------------*/

class surfaceAreaReport
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        // Read from dictonary

            word patchName_;

            label patchID_;

            scalar surfaceArea_;


    // Private Member Functions
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct and assignment
        surfaceAreaReport(const surfaceAreaReport&);
        void operator=(const surfaceAreaReport&);


public:

    //- Runtime type information
    TypeName("surfaceAreaReport");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        surfaceAreaReport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~surfaceAreaReport(){};


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
