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
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "massInOutReport/massInOutReport.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massInOutReport, 0);
    addToRunTimeSelectionTable(functionObject, massInOutReport, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::massInOutReport::writeFileHeader(Ostream& os) const
{
    Log << "    Logging surface statistics to file: " << fileName_
        << endl;

    writeCommented(os, "Time");

    writeDelimited(os,"surfaceAreaIn");

    writeDelimited(os,"surfaceAreaOut");

    if (compressible_)
    {
        writeDelimited(os,"massFluxIn");

        writeDelimited(os,"massFluxOut");
    }
    else
    {
        writeDelimited(os,"volumeFluxIn");

        writeDelimited(os,"volumeFluxOut");
    }

    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionObjects::massInOutReport::massInOutReport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    faceZoneName_("none"),
    phiName_("none"),
    compressible_(false),
    massIn_(0),
    massOut_(0),
    surfaceIn_(0),
    surfaceOut_(0)
{
    read(dict);
    writeFileHeader(file());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool functionObjects::massInOutReport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    //determine compressible/incompressible
    {
        if (dict.found("faceZoneName"))
        {
            dict.lookup("faceZoneName") >> faceZoneName_;
            label faceZoneID = mesh_.faceZones().findZoneID(faceZoneName_);
            if (faceZoneID != -1)
            {
               // ideally check if faceZone is oriented and give warning
            }
            else
            {
                WarningInFunction
                    << faceZoneName_ << " faceZone does not exist"
                    << endl;
            }
        }
        else
        {
            WarningInFunction
                << faceZoneName_ << " could not be found in the database, "
                << "deactivating."
                << endl;
        }

        if (dict.found("phi"))
        {
            dict.lookup("phi") >> phiName_;
        }
        else
        {
            phiName_="phi";
        }
        if
        (
            !foundObject<surfaceScalarField>(phiName_)
        )
        {
            WarningInFunction
                << phiName_ << " could not be found in the database, "
                << "deactivating."
                << endl;
        }
        else if
        (
            lookupObject<surfaceScalarField>(phiName_).dimensions()
            == dimDensity*dimVelocity*dimArea
        )
        {
            compressible_ = true;
        }

    }

    return true;
}


bool Foam::functionObjects::massInOutReport::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    //calculate
    calculate();

    //Write to screen
    Log <<     "SurfaceArea [m^2]     In | Out = ";
    Log << surfaceIn_ << " | " << surfaceOut_ << endl;
    if (compressible_)
    {
        Log << "Massflux [kg/s]       In | Out = ";
    }
    else
    {
        Log << "VolumeFlux [m^3/s]    In | Out = ";
    }
    Log << massIn_ << " | " << massOut_ << endl;

    file() << obr_.time().timeName() << tab
           << surfaceIn_ << tab << surfaceOut_ << tab
           << massIn_ << tab << massOut_
           << endl;

    Log << endl;

    return true;
}


void functionObjects::massInOutReport::calculate()
{
    label faceZoneID = mesh_.faceZones().findZoneID(faceZoneName_);

    const surfaceScalarField& phi =
        obr_.lookupObject<surfaceScalarField>(phiName_);

    surfaceIn_ = 0;
    surfaceOut_ = 0;
    massIn_ = 0;
    massOut_ = 0;

    if (faceZoneID != -1)
    {
        const faceZone& fz = mesh_.faceZones()[faceZoneID];

        forAll(fz, zfI)
        {

            label gfI =  fz[zfI];

            scalar phiTmp = 0;
            scalar surfaceTmp = 0;
            if (gfI < mesh_.nInternalFaces())
            {
                phiTmp = phi[gfI];
                surfaceTmp = mesh_.magSf()[gfI];
            }
            else
            {
                label patchIndex = mesh_.boundaryMesh().whichPatch(gfI);
                label prfI = gfI-mesh_.boundaryMesh()[patchIndex].start();
                phiTmp = phi.boundaryField()[patchIndex][prfI];
                surfaceTmp = mesh_.magSf().boundaryField()[patchIndex][prfI];
            }
            if (fz.flipMap()[zfI])
            {
                if (phiTmp<0)
                {
                    massIn_ -= phiTmp;
                    surfaceIn_ += surfaceTmp;
                }
                else
                {
                    massOut_ += phiTmp;
                    surfaceOut_ += surfaceTmp;
                }
            }
            else
            {
                if (phiTmp<0)
                {
                    massOut_ -= phiTmp;
                    surfaceOut_ += surfaceTmp;
                }
                else
                {
                    massIn_ += phiTmp;
                    surfaceIn_ += surfaceTmp;
                }
            }
        }
    }
    else
    {
        WarningInFunction
            << faceZoneName_ << " faceZone does not exist"
            << endl;
    }
    reduce(surfaceIn_, sumOp<scalar>());
    reduce(surfaceOut_, sumOp<scalar>());
    reduce(massIn_, sumOp<scalar>());
    reduce(massOut_, sumOp<scalar>());
}

bool functionObjects::massInOutReport::write()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
