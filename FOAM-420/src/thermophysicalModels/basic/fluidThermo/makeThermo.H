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
    (c) 2012-2016 OpenFOAM Foundation

InClass
    Foam::fluidThermo

Description
    Macros for creating basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "fluidThermo/fluidThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermoTypedefs(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)\
                                                                               \
typedef                                                                        \
    Transport                                                                  \
    <                                                                          \
        species::thermo                                                        \
        <                                                                      \
            Thermo                                                             \
            <                                                                  \
                EqnOfState                                                     \
                <                                                              \
                    Specie                                                     \
                >                                                              \
            >,                                                                 \
            Type                                                               \
        >                                                                      \
    > Transport##Type##Thermo##EqnOfState##Specie;                             \
                                                                               \
typedef                                                                        \
    Cthermo                                                                    \
    <                                                                          \
        BaseThermo,                                                            \
        Mixture<Transport##Type##Thermo##EqnOfState##Specie>                   \
    > Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie;           \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    (                                                                          \
        #Cthermo"<"#Mixture"<"                                                 \
      + Transport##Type##Thermo##EqnOfState##Specie::typeName()                \
      + ">>"                                                                   \
    ).c_str(),                                                                 \
    0                                                                          \
);


#define makeThermo(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)\
                                                                               \
makeThermoTypedefs                                                             \
(                                                                              \
    BaseThermo,                                                                \
    Cthermo,                                                                   \
    Mixture,                                                                   \
    Transport,                                                                 \
    Type,                                                                      \
    Thermo,                                                                    \
    EqnOfState,                                                                \
    Specie                                                                     \
)                                                                              \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    basicThermo,                                                               \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    objectRegistry                                                             \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    fluidThermo,                                                               \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    objectRegistry                                                             \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    objectRegistry                                                             \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
