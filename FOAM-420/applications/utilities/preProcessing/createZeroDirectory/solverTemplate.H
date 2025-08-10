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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::solverTemplate

Description
    Class to store solver template specifications

\*---------------------------------------------------------------------------*/

#ifndef solverTemplate_H
#define solverTemplate_H

#include "primitives/bools/lists/boolList.H"
#include "primitives/strings/lists/wordList.H"
#include "dimensionSet/dimensionSet.H"
#include "db/IOobject/IOobject.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                       Class solverTemplate Declaration
\*---------------------------------------------------------------------------*/

class solverTemplate
{
public:

    // Public enumerations

        //- Solver type
        enum solverType
        {
            stCompressible,
            stIncompressible,
            stBuoyant,
            stUnknown
        };

        //- Solver type names
        static const NamedEnum<solverType, 4> solverTypeNames_;


private:

    // Private data

        //- Solver type
        solverType solverType_;

        //- Multi-region flag
        bool multiRegion_;


        // Per-region info

            //- Region types
            wordList regionTypes_;

            //- Region names
            wordList regionNames_;

            //- Field names
            List<wordList> fieldNames_;

            //- Field types
            List<wordList> fieldTypes_;

            //- Field dimensions
            List<PtrList<dimensionSet>> fieldDimensions_;


    // Public member functions

        //- Read a word from a dictionary (offers some protection...)
        word readFromDict
        (
            IOobject& dictHeader,
            const word& entryName
        ) const;

        //- Read fluid region templates
        dictionary readFluidFieldTemplates
        (
            const word& regionName,
            const fileName& baseDir,
            const dictionary& solverDict,
            const Time& runTime
        ) const;

        //- Read solid region templates
        dictionary readSolidFieldTemplates
        (
            const word& regionName,
            const dictionary& solverDict
        ) const;

        //- Set the properties for region with index regionI
        void setRegionProperties
        (
            const dictionary& dict,
            const word& regionType,
            const word& regionName,
            const label regionI
        );


public:

    //- Constructor
    solverTemplate
    (
        const fileName& baseDir,
        const Time& runTime,
        const word& regionName
    );


    // Public member functions

        //- Solver type name
        word type() const;

        //- Return the multi-region flag
        bool multiRegion() const;

        //- Return the number of regions
        label nRegion() const;


        // Per-region info

            //- Return the region type
            const word& regionType(const label regionI) const;

            //- Return the region name
            const word& regionName(const label regionI) const;

            //- Return the field names
            const wordList& fieldNames(const label regionI) const;

            //- Return the field types
            const wordList& fieldTypes(const label regionI) const;

            //- Return the field dimensions
            const PtrList<dimensionSet>& fieldDimensions
            (
                const label regionI
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
