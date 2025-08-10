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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "linearLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearLaw, 0);
    addToRunTimeSelectionTable(materialModel, linearLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearLaw::linearLaw
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    linearLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::linearLaw::incompressible() const
{
    if (p_->isConst())
    {
        return true;
    }
    return false;
}


bool Foam::linearLaw::isochoric() const
{
    if (p_->isConst())
    {
        return true;
    }
    return false;
}


Foam::baseModels<Foam::scalar>* Foam::linearLaw::castScalarModel
(
    const word& modelName
)
{
    if (modelName == rhoModel::typeName)
    {
        return dynamic_cast<rhoModel*>(this);
    }
    else if (modelName == SDepartureModel::typeName)
    {
        return dynamic_cast<SDepartureModel*>(this);
    }
    else if (modelName == psiModel::typeName)
    {
        return dynamic_cast<psiModel*>(this);
    }
    else if (modelName == ZModel::typeName)
    {
        return dynamic_cast<ZModel*>(this);
    }
    else if (modelName == CpMCvModel::typeName)
    {
        return dynamic_cast<CpMCvModel*>(this);
    }
    return nullptr;
}


bool Foam::linearLaw::read()
{
    p_ = constructOrReturnRefFieldPtr<scalar>("p");
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    psi_ = dict_->lookup<scalar>("psi");
    rho0_ = dict_->lookup<scalar>("rho0");

    return true;
}


// ************************************************************************* //
