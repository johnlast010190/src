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
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fixedTemperatureConstraint.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(fixedTemperatureConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            fixedTemperatureConstraint,
            dictionary
        );
    }

    template<>
    const char* NamedEnum<fv::fixedTemperatureConstraint::temperatureMode, 3>::
    names[] =
    {
        "uniform",
        "lookup",
        "spatial"
    };
}

const Foam::NamedEnum<Foam::fv::fixedTemperatureConstraint::temperatureMode, 3>
    Foam::fv::fixedTemperatureConstraint::temperatureModeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperatureConstraint::fixedTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    coorFramePtr_(nullptr)
{
    fixedTemperatureConstraint::read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fixedTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const scalar t = mesh_.time().value();

    scalarField THEvalue(cells_.size(), 0.0);
    switch (mode_)
    {
        case tmUniform:
        {
            THEvalue = scalarField(cells_.size(), Tfun1_().value(t));
            break;
        }
        case tmLookup:
        {
            const volScalarField& T =
                mesh_.lookupObject<volScalarField>(TName_);
            THEvalue = scalarField(T, cells_);
            break;
        }
        case tmSpatial:
        {
            vectorField C(mesh_.C().internalField(), cells_);
            if (coorFramePtr_)
            {
                C = coorFramePtr_->coorSys().localPosition(C);
            }
            THEvalue = Tfun1_().value(C);
            break;
        }
    }

    if (thermoPtr_)
    {
        const scalarField pCells(thermoPtr_->p(), cells_);
        THEvalue = thermoPtr_->he(pCells, THEvalue, cells_);
    }

    // Optionally apply only a fraction of the constraint
    if (fraction_.valid())
    {
        eqn.setValues
        (
            cells_,
            THEvalue,
            scalarList(cells_.size(), fraction_->value(t))
        );
    }
    else
    {
        eqn.setValues(cells_, THEvalue);
    }
}


bool Foam::fv::fixedTemperatureConstraint::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        TName_ = coeffs_.lookupOrDefault<word>("T", "T");
        thermoPtr_ =
            obr_.lookupObjectRefPtr<basicThermo>(basicThermo::dictName);

        applied_.setSize(1, false);
        fieldNames_.setSize(1, thermoPtr_ ? thermoPtr_->he().name() : TName_);

        mode_ = temperatureModeNames_.read(coeffs_.lookup("mode"));

        if (mode_ != tmLookup)
        {
            Tfun1_ = Function1<scalar>::New("temperature", coeffs_);
        }

        if (coeffs_.found("referenceFrame"))
        {
            coorFramePtr_ = coordinateFrame::lookupNew(mesh_, coeffs_);
        }

        fraction_ =
            coeffs_.found("fraction")
          ? Function1<scalar>::New("fraction", coeffs_)
          : autoPtr<Function1<scalar>>();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
