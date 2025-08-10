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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.

Class
    Foam::functionObjects::ccnvReport

Group
    grpReportFunctionObjects

Description

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    ccnvReport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ccnvReport_H
#define functionObjects_ccnvReport_H

#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
//#include "db/IOstreams/Fstreams/OFstream.H"
#include <iostream>
#include <fstream>
#include "containers/Lists/PackedList/PackedList.H"
#include "fieldInstance/fieldInstance.H"
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class ccnvReport Declaration
\*---------------------------------------------------------------------------*/

class ccnvReport
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        // Read from dictonary
            //- CCNV precision (4 single and 8 double)
            long ccnv_pr_;

            //- velocity field name
            word UName_;

            //- phi field name
            word PhiName_;

            //- pressure field name
            word PName_;

            //- step number
            int stepNumber_;

            //- input
            Switch logToFile_;

            //- log file name
            fileName logFileName_;

            //- log file ptr
            std::ofstream logFilePtr_;

        //- Result values
            double maxDivergence_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        ccnvReport(const ccnvReport&);
        void operator=(const ccnvReport&);


public:

    //- Runtime type information
    TypeName("ccnvReport");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        ccnvReport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~ccnvReport();


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
