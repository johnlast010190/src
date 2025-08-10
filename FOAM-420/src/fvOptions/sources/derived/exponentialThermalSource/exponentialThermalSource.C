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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "exponentialThermalSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "transportModel/transportModel.H"
#include "fluidThermo/fluidThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(exponentialThermalSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        exponentialThermalSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::exponentialThermalSource::exponentialThermalSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    eName_(coeffs_.lookupOrDefault<word>("eName", "h")),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    T0_(coeffs_.lookupOrDefault<scalar>("T0", 353.0 )),
    Cm_(coeffs_.lookupOrDefault<scalar>("Cm", 0.291e6 )),
    Ce_(coeffs_.lookupOrDefault<scalar>("Ce", 1.369 ))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::exponentialThermalSource::initialise()
{
    if (selectionMode_ != smCellZone)
    {
        FatalErrorInFunction
            << "The porosity region must be specified as a cellZone.  Current "
            << "selection mode is " << selectionModeTypeNames_[selectionMode_]
            << exit(FatalError);
    }

    fieldNames_.setSize(1, eName_);
    applied_.setSize(1, false);

    return true;
}


void Foam::fv::exponentialThermalSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    const volVectorField& U
    (
        obr_.lookupObject<volVectorField>(UName_)
    );

    scalarField& Esource = eqn.source();
    scalarField& Ediag = eqn.diag();
    const scalarField& V = mesh_.V();

    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
    // compressible
    Info<<"exponential thermal source - h eqn."<<endl;
    const fluidThermo& thermo
        (obr_.lookupObject<fluidThermo>(basicThermo::dictName));

    tmp<volScalarField> Cp(thermo.Cp());

    scalar Qtotal = 0;

    forAll(cells_, i)
    {
        label celli = cells_[i];

        scalar h = V[celli] * Cm_*(1.0-exp(-mag(U[celli])/Ce_));

        Qtotal += h*(T0_ - eqn.psi()[celli]/Cp()[celli]);

        Ediag[celli] -= h/Cp()[celli];
        Esource[celli] -= h*T0_;
    }

    reduce(Qtotal, sumOp<scalar>());

    Info<< "Heat transfer from " << name_ << ": "
         << Qtotal << " [W]" << endl;
    }
    else if (obr_.foundObject<transportModel>("transportProperties"))
    {
    // incompressible
    Info<<"exponential thermal source - T eqn."<<endl;
    const transportModel& transport =
            obr_.lookupObject<transportModel>("transportProperties");

    tmp<volScalarField> Cp(transport.Cp());
    tmp<volScalarField> rho(transport.rho());

    scalar Qtotal = 0;

    forAll(cells_, i)
    {
        label celli = cells_[i];

        scalar h = V[celli] * Cm_*(1.0-exp(-mag(U[celli])/Ce_));

        Qtotal += h*(T0_ - eqn.psi()[celli]);

        Ediag[celli] -= h/(Cp()[celli]*rho()[celli]);

        Esource[celli] -= h*T0_/(rho()[celli]*Cp()[celli]);
    }

    reduce(Qtotal, sumOp<scalar>());

    Info<< "Heat transfer from " << name_ << ": "
         << Qtotal << " [W]" << endl;

    }
    else
    {
        FatalErrorInFunction
            << "Could not access Cp field" << nl
            << "Supported objects to access Cp are: fluidThermo and transportProperties"
            << exit(FatalError);
    }
}

void Foam::fv::exponentialThermalSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn, fieldI);
}

void Foam::fv::exponentialThermalSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::exponentialThermalSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("UName", UName_);
        coeffs_.readIfPresent("eName", eName_);
        coeffs_.readIfPresent("T0", T0_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
