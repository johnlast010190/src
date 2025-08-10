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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef makeChemistryReductionMethods_H
#define makeChemistryReductionMethods_H

#include "chemistryModel/TDACChemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethod.H"

#include "chemistryModel/TDACChemistryModel/reduction/noChemistryReduction/noChemistryReduction.H"
#include "chemistryModel/TDACChemistryModel/reduction/DAC/DAC.H"
#include "chemistryModel/TDACChemistryModel/reduction/DRG/DRG.H"
#include "chemistryModel/TDACChemistryModel/reduction/DRGEP/DRGEP.H"
#include "chemistryModel/TDACChemistryModel/reduction/EFA/EFA.H"
#include "chemistryModel/TDACChemistryModel/reduction/PFA/PFA.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryReductionMethod(SS, Comp, Thermo)                         \
                                                                               \
    typedef chemistryReductionMethods::SS<Comp, Thermo> SS##Comp##Thermo;      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryReductionMethod<Comp, Thermo>::                                   \
        adddictionaryConstructorToTable<SS##Comp##Thermo>                      \
        add##chemistryReductionMethods##SS##Comp##Thermo##ConstructorToTable_;


#define makeChemistryReductionMethods(CompChemModel, Thermo)                    \
                                                                               \
    typedef chemistryReductionMethod<CompChemModel, Thermo>                    \
        chemistryReductionMethod##CompChemModel##Thermo;                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryReductionMethod##CompChemModel##Thermo,                       \
        "chemistryReductionMethod<"#CompChemModel","#Thermo">",                \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryReductionMethod##CompChemModel##Thermo,                       \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        none,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        DAC,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        DRG,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        DRGEP,                                                                 \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        EFA,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        PFA,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
