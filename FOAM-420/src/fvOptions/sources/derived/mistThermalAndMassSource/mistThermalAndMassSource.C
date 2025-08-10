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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2021 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "mistThermalAndMassSource.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(mistThermalAndMassSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        mistThermalAndMassSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::mistThermalAndMassSource::mistThermalAndMassSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr),
    mistName_(coeffs_.lookupOrDefault<word>("mistName", "mist")),
    humidityName_(coeffs_.lookupOrDefault<word>("humidityName", "w")),
    TName_(coeffs_.lookupOrDefault<word>("TName", "T")),
    eName_(coeffs_.lookupOrDefault<word>("eName", "h")),
    pName_(coeffs_.lookupOrDefault<word>("pName", "p")),
    phiName_(coeffs_.lookupOrDefault<word>("phiName", "phi")),
    Mvap_(18.02),
    Mair_(28.96),
    Dc_
    (
        dimensionedScalar
        (
            "Dc",
            dimArea/dimTime,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>("Dc", 0.0)
        )
    ),
    rhod_
    (
        dimensionedScalar
        (
            "rhod",
            dimDensity,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>("liqDensity", 1000)
        )
    ),
    diam_
    (
        dimensionedScalar
        (
            "diam",
            dimLength,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>("diam", 0.0005)
        )
    ),
    dHevap_
    (
        dimensionedScalar
        (
            "dHevap",
            dimLength*dimLength/dimTime/dimTime,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>("dHevap", 2.257e+6)
        )
    ),
    Urel_
    (
        dimensionedScalar
        (
            "terminalVelocity",
            dimLength/dimTime,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>
            (
                "terminalVelocity",
                0.0
            )
        )
    ),
    humSat_
    (
        dimensionedScalar
        (
            "humSat",
            dimless,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>("humSat", 0.0)
        )
    ),
    Sc_
    (
        dimensionedScalar
        (
            "Schmidt",
            dimless,
            coeffs_.subDict("mistProperties").lookupOrDefault<scalar>("Schmidt", 0.7)
        )
    )
{
    coeffs_.lookup("fields") >> fieldNames_;

    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::mistThermalAndMassSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == TName_)
    {
        if (debug)
        {
            Info<<"Add thermal source terms into transport Equation for "
                << eqn.psi().name() << endl;
        }
        eqn -= mDot() * dHevap_/getCp();
    }
    else if (eqn.psi().name() == mistName_)
    {
        if (debug)
        {
            Info<<"Add mass source term into transport Equation for "
                << eqn.psi().name() << endl;
            Info<<"min mass source: " << min(mDot()) << "\n"
                <<"max mass source: " << max(mDot()) << endl;
        }
        eqn -= mDot();
    }
    else if (eqn.psi().name() == humidityName_)
    {
        if (debug)
        {
            Info<<"Inserting mass source term into transport Equation for "
                << eqn.psi().name() << endl;
            Info<<"min mass source: " << min(mDot()) << "\n"
                <<"max mass source: " << max(mDot()) << endl;
        }
        eqn += mDot();
    }
}

void Foam::fv::mistThermalAndMassSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == eName_)
    {
        if (debug)
        {
            Info<<"Add thermal source terms into "<<
                "transport Equation for field " << eqn.psi().name() << endl;
        }
        eqn -= mDot() * dHevap_*getRho();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::iad() const
{
    const volScalarField& mist =
        obr_.lookupObject<volScalarField>(mistName_);

    return 6 * max(mist, scalar(0.0)) / diam_ * getRho()/rhod_;
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::Sh() const
{
    return (scalar(2.0) + 0.6*sqrt(Re())*cbrt(Sc_));
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::deltaC() const
{
    const volScalarField& humidity =
        obr_.lookupObject<volScalarField>(humidityName_);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "deltaC",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            max(humSat() - humidity, scalar(0.0))
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::mDot() const
{
    return Sh()*massDCoeff()/diam_*deltaC()*iad();
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::getRho() const
{
    if (obr_.foundObject<volScalarField>("rho"))
    {
        tmp<volScalarField> rho(obr_.lookupObject<volScalarField>("rho"));

        return rho;
    }
    else if (obr_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            obr_.lookupObject<transportModel>("transportProperties");

       tmp<volScalarField> rho(transport.rho());

       return rho;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get rho()" << nl
            << "Supported objects are: basicThermo and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::getCp() const
{
    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo
            (obr_.lookupObject<fluidThermo>(basicThermo::dictName));

        tmp<volScalarField> Cp(thermo.Cp());

        return Cp;
    }
    else if (obr_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            obr_.lookupObject<transportModel>("transportProperties");

        tmp<volScalarField> Cp(transport.Cp());

        return Cp;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get Cp()" << nl
            << "Supported objects are: basicThermo and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::getPabs() const
{
    tmp<volScalarField> pField(obr_.lookupObject<volScalarField>(pName_));

    if (obr_.foundObject<transportModel>("transportProperties"))
    {
        const IOdictionary& transportProperties
            = obr_.lookupObject<IOdictionary>("transportProperties");

        dimensionedScalar pRef("pRef",transportProperties.lookup("pRef"));

        return pField()*getRho() + pRef;
    }
    else if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const dimensionedScalar pRef
        (
            obr_.lookupObject<fluidThermo>(basicThermo::dictName).pRef()
        );
        return pField() + pRef;
    }
    else
    {
        FatalErrorInFunction
            << "Cannot calculate the absolute pressure"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::humSat() const
{
    if (humSat_.value()!=0.0)
    {
        return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "humSat",
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                     ),
                    mesh_,
                    humSat_
                 )
             );
    }
    else
    {

        tmp<volScalarField> T(obr_.lookupObject<volScalarField>(TName_));

        // Tetens'equation for saturation pressure in Pa
        // requires temperature in degree C

        dimensionedScalar oneK("oneK", dimTemperature, 1.0);
        dimensionedScalar cPa("cPa", dimMass/(dimLength*dimTime*dimTime), 610.78);
        tmp<volScalarField> pSat
        (
            cPa*Foam::exp((17.27*(T()/oneK-273.15))/((T/oneK-273.15)+237.3))
        );

        return (Mvap_/Mair_)*pSat/getPabs();
    }

}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::Re() const
{
    const surfaceScalarField& phi =
        obr_.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        const incompressible::turbulenceModel& turbModel =
            obr_.lookupObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        return Urel_*diam_/turbModel.nu();
    }
    else
    {
        const compressible::turbulenceModel& turbModel =
            obr_.lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
             );

        return Urel_*diam_/turbModel.nu();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::fv::mistThermalAndMassSource::massDCoeff() const
{
    if (Dc_.value()!=0.0)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Dc",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                Dc_
            )
        );
    }
    else
    {
        if (obr_.foundObject<transportModel>("transportProperties"))
        {
            const IOdictionary& transportProperties
                = obr_.lookupObject<IOdictionary>("transportProperties");

            dimensionedScalar Dlw("Dw",transportProperties.lookup("Dw"));

            // tmp<volScalarField> Dtw(obr_.lookupObject<volScalarField>("Dtw"));
            // return Dlw + Dtw;

            // return laminar diffusion in concentration boundary layer
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "Dlw",
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    Dlw
                )
            );
        }
        else if (obr_.foundObject<IOdictionary>(basicThermo::dictName))
        {
            const IOdictionary& thermophysicalProperties =
                obr_.lookupObject<IOdictionary>(basicThermo::dictName);

            // return laminar diffusion in concentration boundary layer
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "Dlw",
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    (basicThermo::dictName == basicThermo::matDictName)
                  ? dimensionedScalar
                    (
                        "D",
                        thermophysicalProperties.optionalSubDict
                        (
                            humidityName_ + "concentrationTransport"
                        ).lookup("D")
                    )
                  : dimensionedScalar
                    (
                        "Dw",thermophysicalProperties.lookup("Dw")
                    )
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "No access to Dw from mistThermalAndMassSource"
                << exit(FatalError);

            return tmp<volScalarField>(nullptr);
        }
    }
}

// ************************************************************************* //
