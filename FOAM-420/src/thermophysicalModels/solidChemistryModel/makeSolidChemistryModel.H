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

Description
    Macros for instantiating solid chemistry models

\*---------------------------------------------------------------------------*/

#ifndef makeSolidChemistryModel_H
#define makeSolidChemistryModel_H

#include "solidChemistryModel/solidChemistryModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidChemistryModel(sChemistry, SS, Comp, SThermo)                 \
                                                                               \
    typedef Foam::sChemistry<Foam::Comp, Foam::SThermo>                        \
        sChemistry##Comp##SThermo;                                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        sChemistry##Comp##SThermo,                                             \
        (Foam::word(sChemistry##Comp##SThermo::typeName_()) + "<"#Comp","      \
        + SThermo::typeName() + ">").c_str(),                                  \
        0                                                                      \
    );


#define makeSolidGasChemistryModel(sChemistry, SS, Comp, SThermo, GThermo)     \
                                                                               \
    typedef Foam::SS<Foam::Comp, Foam::SThermo, Foam::GThermo>                 \
        SS##Comp##SThermo##GThermo;                                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##SThermo##GThermo,                                            \
        (Foam::word(SS##Comp##SThermo##GThermo::typeName_()) + "<"#Comp","     \
        + SThermo::typeName() + "," + GThermo::typeName() + ">").c_str(),      \
        0                                                                      \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
