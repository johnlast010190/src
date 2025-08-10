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

#include "icoPolynomialLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(icoPolynomialLaw, 0);
    addToRunTimeSelectionTable(materialModel, icoPolynomialLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::icoPolynomialLaw::icoPolynomialLaw
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
    icoPolynomialLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::icoPolynomialLaw::castScalarModel
(
    const word& modelName
)
{
    if (modelName == rhoModel::typeName)
    {
        return dynamic_cast<rhoModel*>(this);
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


bool Foam::icoPolynomialLaw::read()
{
    p_ = constructOrReturnRefFieldPtr<scalar>("p");
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    rhoCoeffs_.reset(dict_->lookup<scalarField>("rhoCoeffs"));
    rhoMin_ = dict_->lookupOrDefault<scalar>("rhoMin", 0.0);
    rhoMax_ = dict_->lookupOrDefault<scalar>("rhoMax", GREAT);

    return true;
}


// ************************************************************************* //
