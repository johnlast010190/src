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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2013-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "DESModelRegions/DESModelRegions.H"
#include "fields/volFields/volFields.H"
#include "DES/DESModel/DESModelBase.H"
#include "turbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(DESModelRegions, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        DESModelRegions,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::DESModelRegions::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "DES model region coverage (% volume)");

    writeCommented(os, "Time");
    writeDelimited(os, "LES");
    writeDelimited(os, "RAS");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::DESModelRegions::DESModelRegions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    resultName_(name)
{
    read(dict);

    tmp<volScalarField> tDESModelRegions
    (
        (
            new volScalarField
            (
                IOobject
                (
                    resultName_,
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0", dimless, 0.0)
            )
        )
    );

    store(resultName_, tDESModelRegions);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::DESModelRegions::~DESModelRegions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::DESModelRegions::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    dict.readIfPresent("result", resultName_);

    return true;
}


bool Foam::functionObjects::DESModelRegions::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    volScalarField& DESModelRegions =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(resultName_)
        );


    if (foundObject<DESModelBase>(turbulenceModel::propertiesName))
    {
        const DESModelBase& model =
            lookupObject<DESModelBase>
            (
                turbulenceModel::propertiesName
            );

        DESModelRegions.forceAssign(model.LESRegion());

        scalar prc =
            gSum(DESModelRegions.primitiveField()*mesh_.V())
           /gSum(mesh_.V())*100.0;

        file() << time_.timeName()
            << token::TAB << prc
            << token::TAB << 100.0 - prc
            << endl;

        Log << "    LES = " << prc << " % (volume)" << nl
            << "    RAS = " << 100.0 - prc << " % (volume)" << nl
            << endl;
    }
    else
    {
        Log << "    No DES turbulence model found in database" << nl
            << endl;
    }

    return true;
}


bool Foam::functionObjects::DESModelRegions::write()
{
    const volScalarField& DESModelRegions =
        lookupObject<volScalarField>(resultName_);

    Log << type() << " " << name() <<  " output:" << nl
        << "    writing field " << DESModelRegions.name() << nl
        << endl;

    DESModelRegions.write();

    return true;
}


// ************************************************************************* //
