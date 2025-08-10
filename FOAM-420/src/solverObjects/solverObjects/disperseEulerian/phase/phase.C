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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "phase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{
    defineTypeNameAndDebug(phase, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decoupledEulerian::phase::phase
(
    const word& name,
    const dictionary& dict,
    const objectRegistry& obr,
    const fvMesh& mesh,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            "phase" + Foam::word(dict.name()),
            mesh.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    name_(name),
    dict_(dict),
    obr_(obr),
    mesh_(mesh),
    rhoc_(),
    muc_(),
    rhocVal_(dimensionedScalar("rhoc", dimDensity, -1)),
    mucVal_(dimensionedScalar("muc", dimensionSet(1,-1,-1,0,0), -1)),
    kappac_(dimensionedScalar("kappac", dimArea/dimTime, 0)),
    Pr_(0),
    Tsat_(dimensionedScalar("Tsat", dimTemperature, 0)),
    L_(dimensionedScalar("L", dimVelocity*dimVelocity, 0)),     // given in kJ/kg = 1000 m^2/s^2, needs to be converted to m^2/s^2
    Dc_(dimensionedScalar("Dc", dimArea/dimTime, 0)),
    alphad_(),
    Ud_(),
    rhod_(dict.lookup("rhod")),
    diam_(dict.lookup("diam")),
    Td_(dimensionedScalar("Td", dimTemperature, 0)),
    Md_(dimensionedScalar("Md", dimensionSet(1,0,0,0,-1), 0)),  // in kg/kmol //Note that R is given in J/kmol/K = kg m^2/(s^2 kmol K)
    dHevap_(dimensionedScalar("dHevap", dimTemperature, 0))     // dH/rho/cp
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decoupledEulerian::phase::~phase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::decoupledEulerian::phase::read(const dictionary& dict)
{
    if (dict.found("rhoc"))
    {
        rhocVal_ = dict.lookup("rhoc");
    }
    if (dict.found("muc"))
    {
        mucVal_ = dict.lookup("muc");
    }

    kappac_ = dimensionedScalar("kappac", dimArea/dimTime, dict.lookupOrDefault<scalar>("kappac", 0));
    Pr_ = dict.lookupOrDefault<scalar>("Pr", 0);
    Tsat_ = dimensionedScalar("Tsat", dimTemperature, dict.lookupOrDefault<scalar>("Tsat", 0));
    L_ = dimensionedScalar("L", dimVelocity*dimVelocity, dict.lookupOrDefault<scalar>("L", 0)*1000);
    Dc_ = dimensionedScalar("Dc", dimArea/dimTime, dict.lookupOrDefault<scalar>("Dc", 0));
    rhod_ = dict.lookup("rhod");
    diam_ = dict.lookup("diam");
    Td_ = dimensionedScalar("Td", dimTemperature, dict.lookupOrDefault<scalar>("Td", 0));
    Md_ = dimensionedScalar("Md", dimensionSet(1,0,0,0,-1), dict.lookupOrDefault<scalar>("Md", 0));
    dHevap_ = dimensionedScalar("dHevap", dimTemperature, dict.lookupOrDefault<scalar>("dHevap", 0));

    rhoc_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rhoc",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            rhocVal_
        )
    );
    muc_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "muc",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            mucVal_
        )
    );

    alphad_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "alpha" + dict.name(),
                mesh_.time().timeName(),
                obr_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Ud_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U" + dict.name(),
                mesh_.time().timeName(),
                obr_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );
}

void Foam::decoupledEulerian::phase::updateRhoc() const
{
    if (rhocVal_.value() < 0)
    {
        typedef compressible::turbulenceModel cmpTurbModel;
        typedef incompressible::turbulenceModel icoTurbModel;

        if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
        {
            const cmpTurbModel& turb =
                obr_.lookupObject<cmpTurbModel>
                (
                    cmpTurbModel::propertiesName
                );
            rhoc_() = turb.rho();
        }
        else if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
        {
            const icoTurbModel& turb =
                obr_.lookupObject<icoTurbModel>
                (
                    icoTurbModel::propertiesName
                );
            rhoc_() = turb.rho();
        }
    }
}

void Foam::decoupledEulerian::phase::updateMuc() const
{
    if (mucVal_.value() < 0)
    {
        const turbulenceModel& turb =
            obr_.lookupObject<turbulenceModel>(turbulenceModel::propertiesName);

        muc_() = turb.nu()*rhoc();
    }
}

void Foam::decoupledEulerian::phase::getRhod()
{
    rhod_ = dict_.lookup("rhod");
}

void Foam::decoupledEulerian::phase::getDiam()
{
    diam_ = dict_.lookup("diam");
}

bool Foam::decoupledEulerian::phase::writeData(Ostream& os) const
{
    return os.good();
}

Foam::tmp<Foam::volScalarField> Foam::decoupledEulerian::phase::Re() const
{
    return rhoc()*magUr()*diam_/muc();
}

Foam::tmp<Foam::volScalarField> Foam::decoupledEulerian::phase::magUr() const
{
    return mag(Uc() - Ud());
}

Foam::tmp<Foam::volVectorField> Foam::decoupledEulerian::phase::Ur() const
{
    return Uc() - Ud();
}

// ************************************************************************* //
