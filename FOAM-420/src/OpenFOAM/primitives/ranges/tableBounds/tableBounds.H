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
    (c) 2019 OpenCFD Ltd.

Namespace
    Foam::bounds

Description
    Namespace for bounding specifications.
    At the moment, mostly for tables.

SourceFiles
    tableBounds.C

\*---------------------------------------------------------------------------*/

#ifndef tableBounds_H
#define tableBounds_H

#include "primitives/enums/Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace bounds
{

  // Enumerations

    //- Enumeration for handling out-of-bound values
    enum class normalBounding
    {
        ERROR,  //!< Exit with a FatalError
        WARN,   //!< Issue warning and clamp value (this is a good default)
        CLAMP   //!< Clamp value to the start/end value
    };

    //- Enumeration for handling out-of-bound values that are repeatable
    enum class repeatableBounding
    {
        ERROR,  //!< Exit with a FatalError
        WARN,   //!< Issue warning and clamp value (this is a good default)
        CLAMP,  //!< Clamp value to the start/end value
        REPEAT  //!< Treat as a repeating list
    };


    //- Strings corresponding to the normalBounding
    extern const Foam::Enum<normalBounding> normalBoundingNames;

    //- Strings corresponding to the repeatableBounding
    extern const Foam::Enum<repeatableBounding> repeatableBoundingNames;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bounds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
