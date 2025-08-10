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
    (c) 2016-2018 OpenCFD Ltd.
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::TimePaths

Description
    A class for addressing time paths without using the Time class.

SourceFiles
    TimePaths.C

\*---------------------------------------------------------------------------*/

#ifndef TimePaths_H
#define TimePaths_H

#include "primitives/strings/fileName/fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TimePaths Declaration
\*---------------------------------------------------------------------------*/


class TimePaths
{
    // Private data

        bool processorCase_;
        const fileName rootPath_;
        bool distributed_;
        fileName globalCaseName_;
        fileName case_;
        const word system_;
        const word constant_;


    // Private Member functions

        //- Determine from case name whether it is a processor directory
        bool detectProcessorCase();


public:

    // Constructors

        //- Construct given database name, rootPath and casePath
        TimePaths
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName,
            const word& constantName
        );


        //- Construct given database name, rootPath and casePath
        TimePaths
        (
            const bool processorCase,
            const fileName& rootPath,
            const bool distributed,
            const fileName& globalCaseName,
            const fileName& caseName,
            const word& systemName,
            const word& constantName
        );


    // Member functions

            //- Return true if this is a processor case
            bool processorCase() const
            {
                return processorCase_;
            }

            //- Return root path
            const fileName& rootPath() const
            {
                return rootPath_;
            }

            //- Return global case name
            const fileName& globalCaseName() const
            {
                return globalCaseName_;
            }

            //- Return case name
            const fileName& caseName() const
            {
                return case_;
            }

            //- Return case name
            fileName& caseName()
            {
                return case_;
            }

            //- Return system name
            const word& system() const
            {
                return system_;
            }

            //- Return system name for the case
            //  which for parallel runs returns ../system()
            fileName caseSystem() const;

            //- Return constant name
            const word& constant() const
            {
                return constant_;
            }

            //- Return the input relative to the globalPath by stripping off
            //- a leading value of the globalPath
            //  This is required for compatability with ESI OpenFOAM v1906
            //  \param input the directory or filename to make case-relative
            //  \param caseTag replace globalPath with \<case\> for later
            //      use with expand(), or prefix \<case\> if the file name was
            //      not an absolute location
            inline fileName relativePath
            (
                const fileName& input,
                const bool caseTag = false
            ) const
            {
                return input.relative(rootPath()/globalCaseName(), caseTag);
            }

            //- Is case running with parallel distributed directories
            //  (i.e. not NFS mounted)
            bool distributed() const
            {
                return distributed_;
            }

            //- Return constant name for the case
            //  which for parallel runs returns ../constant()
            fileName caseConstant() const;

            //- Return path
            fileName path() const
            {
                return rootPath()/caseName();
            }

            //- Return system path
            fileName systemPath() const
            {
                return path()/system();
            }

            //- Return constant path
            fileName constantPath() const
            {
                return path()/constant();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
