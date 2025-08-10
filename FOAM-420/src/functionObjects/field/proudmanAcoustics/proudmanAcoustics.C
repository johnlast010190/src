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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "proudmanAcoustics/proudmanAcoustics.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(proudmanAcoustics, 0);
    addToRunTimeSelectionTable(functionObject, proudmanAcoustics, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::proudmanAcoustics::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Proudman's formula");
    writeCommented(os, "Time");
    writeDelimited(os, "min");
    writeDelimited(os, "max");
    os << endl;
}


void Foam::functionObjects::proudmanAcoustics::calcProudmansFormula
(
    volScalarField& acPower,
    const volScalarField& kTurb,
    const volScalarField& eTurb
)
{
    tmp<volScalarField> Mt (sqrt(2*kTurb)
            / dimensionedScalar( "vel", dimLength/dimTime, alpha0_)
                            );

    acPower = alphae_  * eTurb * pow(Mt,5)
                * dimensionedScalar( "rho", dimMass/sqr(dimLength)/dimLength, rhoInf_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::proudmanAcoustics::proudmanAcoustics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    rhoInf_(dict.lookupOrDefault<scalar>("rhoInf",1.0)),
    alpha0_(dict.lookupOrDefault<scalar>("alpha0",343)),
    alphae_(dict.lookupOrDefault<scalar>("alphae",0.1)),
    writedB_(dict.lookupOrDefault<Switch>("writedB",true)),
    refAcPower_(dict.lookupOrDefault<scalar>("refAcPower",1e-12))
{
    read(dict);

    writeFileHeader(file());

    volScalarField* proudmanAcousticsPtr
    (
        new volScalarField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "0",
                dimMass/dimLength/sqr(dimTime)/dimTime,
                0.0
            )
        )
    );

    mesh_.objectRegistry::store(proudmanAcousticsPtr);

    if (writedB_)
    {
        volScalarField* acPowerdBPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "acPowerdB",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "0",
                    dimless,
                    0.0
                )
            )
        );

        mesh_.objectRegistry::store(acPowerdBPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::proudmanAcoustics::~proudmanAcoustics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::proudmanAcoustics::read(const dictionary& dict)
{
    Log << type() << " " << name() <<  " read:" << nl;

    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::proudmanAcoustics::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    volScalarField& proudmanAcoustics =
        const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>(type())
        );

    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    if (mesh_.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model =
            mesh_.lookupObject<cmpModel>(turbulenceModel::propertiesName);

        const tmp<volScalarField> kTurb = model.k();
        const tmp<volScalarField> eTurb = model.epsilon();

        calcProudmansFormula(proudmanAcoustics, kTurb(), eTurb());
    }
    else if (mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model =
            mesh_.lookupObject<icoModel>(turbulenceModel::propertiesName);

        const tmp<volScalarField> kTurb = model.k();
        const tmp<volScalarField> eTurb = model.epsilon();

        calcProudmansFormula(proudmanAcoustics, kTurb(), eTurb());

    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    scalar minPa = gMin(proudmanAcoustics.primitiveField());
    scalar maxPa = gMax(proudmanAcoustics.primitiveField());

    if (Pstream::master())
    {
        writeTime(file());

        file()
            << token::TAB << minPa
            << token::TAB << maxPa
            << endl;
    }

    Log << "    min/max (W/m3) = "
            << minPa << ", " << maxPa << endl;


    if (writedB_)
    {
        volScalarField& acPowerdB =
            const_cast<volScalarField&>
            (
                mesh_.lookupObject<volScalarField>("acPowerdB")
            );

        acPowerdB = 10 * log10(proudmanAcoustics/
                      dimensionedScalar("pv",dimMass/dimLength/sqr(dimTime)/dimTime,refAcPower_)
                            );

        scalar minPadB = gMin(acPowerdB.primitiveField());
        scalar maxPadB = gMax(acPowerdB.primitiveField());

        if (Pstream::master())
        {
            writeTime(file());

            file()
                << token::TAB << minPadB
                << token::TAB << maxPadB
                << endl;
        }

        Log << "    min/max (dB) = "
                << minPadB << ", " << maxPadB << endl;

    }

    return true;
}


bool Foam::functionObjects::proudmanAcoustics::write()
{
    const volScalarField& proudmanAcoustics =
        obr_.lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << proudmanAcoustics.name() << endl;

    proudmanAcoustics.write();

    if (writedB_)
    {
        const volScalarField& acPowerdB =
            obr_.lookupObject<volScalarField>("acPowerdB");

        Log << type() << " " << name() << " write:" << nl
            << "    writing field acPowerdB" << endl;

        acPowerdB.write();

    }


    return true;
}


// ************************************************************************* //
