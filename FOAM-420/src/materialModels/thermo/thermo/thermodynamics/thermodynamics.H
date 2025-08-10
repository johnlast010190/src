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
    (c) 2021-2022 Esi Ltd.

Class
    Foam::thermodynamics

Description
    Class for selecting material models.

SourceFiles
    thermodynamics.C

\*---------------------------------------------------------------------------*/

#ifndef thermodynamics_H
#define thermodynamics_H

#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class thermodynamics Declaration
\*---------------------------------------------------------------------------*/

class thermodynamics
{

public:

    //- Runtime type information
    TypeName("thermodynamics");

    //- Lists models associated with enthalpy and entropy
    static const wordList listModels()
    {
        return
            wordList
            ({
                limitModel::typeName,
                CpModel::typeName,
                HaModel::typeName,
                HsModel::typeName,
                HcModel::typeName,
                SModel::typeName
            });
    }

    //- Destructor
    virtual ~thermodynamics();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
