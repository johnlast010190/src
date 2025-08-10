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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef makeChemistrySolverTypes_H
#define makeChemistrySolverTypes_H

#include "chemistrySolver/chemistrySolver/chemistrySolver.H"

#include "chemistryModel/chemistryModel/chemistryModel.H"
#include "chemistryModel/TDACChemistryModel/TDACChemistryModel.H"

#include "chemistrySolver/noChemistrySolver/noChemistrySolver.H"
#include "chemistrySolver/EulerImplicit/EulerImplicit.H"
#include "chemistrySolver/ode/ode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistrySolverType(SS, Comp, Thermo)                              \
                                                                               \
    typedef SS<chemistryModel<Comp, Thermo>> SS##Comp##Thermo;                 \
    typedef SS<TDACChemistryModel<Comp, Thermo>> TDAC##SS##Comp##Thermo;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<chemistryModel<" + word(Comp::typeName_())                       \
      + "," + Thermo::typeName() + ">>").c_str(),                              \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Comp,                                                                  \
        SS##Comp##Thermo,                                                      \
        objectRegistry                                                                 \
    );                                                                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        TDAC##SS##Comp##Thermo,                                                \
        (#SS"<TDACChemistryModel<" + word(Comp::typeName_())                   \
      + "," + Thermo::typeName() + ">>").c_str(),                              \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Comp,                                                                  \
        TDAC##SS##Comp##Thermo,                                                \
        objectRegistry                                                                 \
    );


#define makeChemistrySolverTypes(CompChemModel,Thermo)                         \
                                                                               \
    makeChemistrySolverType                                                    \
    (                                                                          \
        noChemistrySolver,                                                     \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistrySolverType                                                    \
    (                                                                          \
        EulerImplicit,                                                         \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistrySolverType                                                    \
    (                                                                          \
        ode,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
