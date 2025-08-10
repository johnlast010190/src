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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::cpuInfo

Description
    General CPU characteristics.

    If the machine has multiple cpus/cores, only the characteristics
    of the first core are used.

Note
    Uses the information from /proc/cpuinfo

SourceFiles
    cpuInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cpuInfo_H
#define cpuInfo_H

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// forward declarations
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class cpuInfo Declaration
\*---------------------------------------------------------------------------*/

class cpuInfo
{
    // Private data

        // Various bits from /proc/cpuinfo

        std::string vendor_id;
        std::string model_name;
        int cpu_family;
        int model;
        float cpu_MHz;
        int siblings;
        int cpu_cores;

    // Private Member Functions

        //- Parse /proc/cpuinfo
        void parse();

        //- Disallow default bitwise assignment
        void operator=(const cpuInfo&) = delete;

        //- Disallow default copy constructor
        cpuInfo(const cpuInfo&) = delete;

public:

    // Constructors

        //- Construct null
        cpuInfo();


    //- Destructor
    ~cpuInfo();


    // Member Functions

        //- Write content as dictionary entries
        void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
