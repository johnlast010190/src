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
    (c) 2010-2023 Esi Ltd.
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "explicitPorositySource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(explicitPorositySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        explicitPorositySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::explicitPorositySource::explicitPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    porosityPtr_(nullptr),
    adjointVelocityPrefix_("Ua_"),
    UName_("U")
{
    read(dict);

    if (selectionMode_ != smCellZone)
    {
        FatalErrorInFunction
            << "selection mode is " << selectionModeTypeNames_[selectionMode_]
            << exit(FatalError);
    }

    porosityPtr_.reset
    (
        porosityModel::New
        (
            name_,
            obr_,
            mesh_,
            coeffs_,
            cellSetName_
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::explicitPorositySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (fieldNames_[fieldi].find(adjointVelocityPrefix_, 0) == string::npos)
    {
        porosityPtr_->addResistance(porosityEqn);
    }
    else
    {
        const volVectorField& Uprimal =
            obr_.lookupObject<volVectorField>(UName_);

        porosityPtr_->addAdjointResistance(porosityEqn, Uprimal);
    }

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    //- Here Assuming the dimensions are the same inside the matrix which is
    //  true for the porosity classes. They are applied on the momentum block
    //  only. It is done to avoid coding more constructors in block
    //  Code can be easily extended if needed in the future
    fvBlockMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensionSets()[0]);

    porosityPtr_->addResistance(porosityEqn);

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (fieldNames_[fieldi].find(adjointVelocityPrefix_, 0) == string::npos)
    {
        porosityPtr_->addResistance(porosityEqn);
    }
    else
    {
        const volVectorField& Uprimal =
            obr_.lookupObject<volVectorField>(UName_);

        porosityPtr_->addAdjointResistance(porosityEqn, Uprimal);
    }

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    //- Here Assuming the dimensions are the same inside the matrix which is
    //  true for the porosity classes. They are applied on the momentum block
    //  only. It is done to avoid coding more constructors in block
    //  Code can be easily extended if needed in the future
    fvBlockMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensionSets()[0]);

    porosityPtr_->addResistance(porosityEqn);

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (fieldNames_[fieldi].find(adjointVelocityPrefix_, 0) == string::npos)
    {
        porosityPtr_->addResistance(porosityEqn);
    }
    else
    {
        const volVectorField& Uprimal =
            obr_.lookupObject<volVectorField>(UName_);

        porosityPtr_->addAdjointResistance(porosityEqn, Uprimal);
    }

    eqn -= alpha*porosityEqn;
}


Foam::label Foam::fv::explicitPorositySource::applyToField
(
    const word& fieldName, const word& regionName
) const
{
    //inject adjoint velocity fields into lists if they are called

    if
    (
        fieldName.find(adjointVelocityPrefix_, 0) != string::npos &&
        findIndex(fieldNames_, fieldName) == -1
    )
    {
        wordList& cfieldNames = const_cast<wordList&>(fieldNames_);
        wordList& cregionNames = const_cast<wordList&>(regionNames_);
        boolList& capplied = const_cast<boolList&>(applied_);
        label fsize = fieldNames_.size();

        cfieldNames.setSize(fsize+1);
        cfieldNames[fsize] = fieldName;

        if (cregionNames.size())
        {
            cregionNames.setSize(fsize+1);
            cregionNames[fsize] =
                (obr_.name() == word::null ? mesh_.name() : obr_.name());
        }

        capplied.setSize(fsize+1);
        capplied[fsize] = true;
    }

    return option::applyToField(fieldName, regionName);
}



bool Foam::fv::explicitPorositySource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (coeffs_.found("UNames"))
        {
            coeffs_.lookup("UNames") >> fieldNames_;
        }
        else if (coeffs_.found("U"))
        {
            UName_ = word(coeffs_.lookup("U"));
            fieldNames_ = wordList(1, UName_);
        }
        else
        {
            fieldNames_ = wordList(1, UName_);
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
