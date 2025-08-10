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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef makeReactionThermo_H
#define makeReactionThermo_H

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "mixtures/SpecieMixture/SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactionMixtureThermo(BaseThermo,CThermo,MixtureThermo,Mixture,ThermoPhys)   \
                                                                               \
typedef MixtureThermo                                                          \
<                                                                              \
    CThermo,                                                                   \
    SpecieMixture                                                              \
    <                                                                          \
        Mixture                                                                \
        <                                                                      \
            ThermoPhys                                                         \
        >                                                                      \
    >                                                                          \
> MixtureThermo##Mixture##ThermoPhys;                                          \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    (#MixtureThermo"<"#Mixture"<" + ThermoPhys::typeName() + ">>").c_str(),    \
    0                                                                          \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    basicThermo,                                                               \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    objectRegistry                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    fluidThermo,                                                               \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    objectRegistry                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    objectRegistry                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    CThermo,                                                                   \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    objectRegistry                                                                     \
);


#define makeMatReactionThermo(BaseThermo,CThermo,MixtureThermo,Mixture,ThermoPhys,MatType)   \
                                                                               \
typedef MixtureThermo                                                          \
<                                                                              \
    CThermo,                                                                   \
    SpecieMixture                                                              \
    <                                                                          \
        Mixture                                                                \
        <                                                                      \
            ThermoPhys                                                         \
        >                                                                      \
    >                                                                          \
> MatType;                                                                     \
                                                                               \
defineTemplateTypeNameAndDebugWithName(MatType, #MatType, 0);                  \
addToRunTimeSelectionTable(basicThermo, MatType, objectRegistry);              \
addToRunTimeSelectionTable(fluidThermo, MatType, objectRegistry);              \
addToRunTimeSelectionTable(BaseThermo, MatType, objectRegistry);               \
addToRunTimeSelectionTable(CThermo, MatType, objectRegistry);


#define makeReactionThermo(BaseThermo,CThermo,MixtureThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
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
makeReactionMixtureThermo                                                      \
(                                                                              \
    BaseThermo,                                                                \
    CThermo,                                                                   \
    MixtureThermo,                                                             \
    Mixture,                                                                   \
    Transport##Type##Thermo##EqnOfState##Specie                                \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
