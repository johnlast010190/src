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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef makeCombustionTypes_H
#define makeCombustionTypes_H

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCombustionTypesThermo(CombModel, CombType, Thermo, Table)          \
                                                                               \
    typedef Foam::combustionModels::CombModel                                  \
        <Foam::combustionModels::CombType, Foam::Thermo>                       \
        CombModel##CombType##Thermo;                                           \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombModel##CombType##Thermo,                                           \
        #CombModel"<"#CombType","#Thermo">",                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace combustionModels                                             \
        {                                                                      \
            typedef CombModel<CombType, Thermo> CombModel##CombType##Thermo;   \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Table,                                                         \
                CombModel##CombType##Thermo,                                   \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeCombustionTypes(CombModel, CombType, Table)                        \
                                                                               \
    typedef Foam::combustionModels::CombModel                                  \
        <Foam::combustionModels::CombType>                                     \
        CombModel##CombType;                                                   \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CombModel##CombType,                                                   \
        #CombModel"<"#CombType">",                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace combustionModels                                             \
        {                                                                      \
            typedef CombModel<CombType> CombModel##CombType;                   \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Table,                                                         \
                CombModel##CombType,                                           \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
