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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "flowCharacteristics/flowCharacteristics.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(flowCharacteristics, 0);
    addToRunTimeSelectionTable(functionObject, flowCharacteristics, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::flowCharacteristics::writeFileHeader(Ostream& os) const
{
    Info<< "    Logging flow charasteristics to file: " << fileName_
        << endl;

    writeCommented(os, "Time");
    writeDelimited(os,"Lvol");
    writeDelimited(os,"Lext");
    writeDelimited(os,"Uvol");
    writeDelimited(os,"UBC");
    writeDelimited(os,"UDP");
    writeDelimited(os,"dtc");
    writeDelimited(os,"dta");

    os << endl;
}


void Foam::functionObjects::flowCharacteristics::calculate()
{
    // Characteristic lenghts
    const scalarField& V = mesh_.V();
    totalVol_ = gSum(V);
    label nDirs = mesh_.nSolutionD();
    const Vector<label>& vectorSol = mesh_.solutionD();

    boundBox bb(mesh_.points(), true);
    point pMin(bb.min());
    point pMax(bb.max());
    point dp = pMax - pMin;
    lExt_ = cmptMax(dp);


    if (nDirs==2)
    {
        scalar depth = 0.0;
        if (vectorSol[0]==-1) depth = pMax[0]-pMin[0];
        if (vectorSol[1]==-1) depth = pMax[1]-pMin[1];
        if (vectorSol[2]==-1) depth = pMax[2]-pMin[2];

        lVol_ = sqrt(totalVol_/depth);
    }
    else if (nDirs==3)
    {
        lVol_ = pow(totalVol_, 1.0/3.0);
    }


    // Velocity scale
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    UVol_ = gMax(mag(U.internalField())());
    UBC_ = gMax(mag(U.boundaryField())());

    const volScalarField& p = mesh_.lookupObject<volScalarField>(pName_);

    const volScalarField::Boundary& pB = p.boundaryField();
    scalar pBCMax = -GREAT;
    scalar pBCMin = GREAT;

    forAll(pB, pI)
    {
        if
        (
            mesh_.boundary()[pI].patch().type() == "patch"
          ||mesh_.boundary()[pI].patch().physicalType() == "patch"
        )
        {
            pBCMax = max(pBCMax, gMax(p.boundaryField()[pI]));
            pBCMin = min(pBCMin, gMin(p.boundaryField()[pI]));
        }
    }

    if (p.dimensions() == dimDensity)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName_);
        scalar rhoAv = gAverage(rho.internalField());
        UDP_ = sqrt((pBCMax - pBCMin)/rhoAv);
    }
    else
    {
        UDP_ = sqrt(pBCMax - pBCMin);
    }


    dtc_ = f_*min(lVol_, lExt_)/max(max(UVol_,UBC_), UDP_);
    dta_ = f_*max(lVol_, lExt_)/max(max(UVol_,UBC_), UDP_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::flowCharacteristics::flowCharacteristics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    rhoName_(dict.lookupOrDefault<word>("rhoName", "rho")),
    f_(dict.lookupOrDefault<scalar>("frequency", 0.3)),
    totalVol_(0.0),
    lVol_(0.0),
    lExt_(0.0),
    UVol_(0.0),
    UBC_(0.0),
    UDP_(0.0),
    dtc_(0.0),
    dta_(0.0)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::flowCharacteristics::~flowCharacteristics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::flowCharacteristics::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    calculate();

    Ostream& rpInfo(Info);

    rpInfo.precision(6);


    Info<< "L_{vol, ext}: " << lVol_ << ", " << lExt_ << " [m]" << endl;

    Info<< "U_{vol, bc, dp}: " << UVol_ << ", "
         << UBC_ << ", " << UDP_ << " [m/s]" << endl;
    Info<< "Dt_{cons, aggr}: " << dtc_ << ", " << dta_ << " [s]" << endl;

    rpInfo.precision(IOstream::defaultPrecision());

    file() << obr_.time().timeName() << tab
           << lVol_ << tab << lExt_ << tab
           << UVol_ << tab << UBC_ << tab << UDP_ << tab
           << dtc_ << tab << dta_
           << endl;

    return true;
}


bool Foam::functionObjects::flowCharacteristics::write()
{
    return true;
}


bool Foam::functionObjects::flowCharacteristics::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
