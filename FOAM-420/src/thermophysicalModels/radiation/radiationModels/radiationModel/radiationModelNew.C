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

#include "radiationModels/radiationModel/radiationModel.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiation::radiationModel>
Foam::radiation::radiationModel::New
(
    const volScalarField& T
)
{
    // check if radiation model already exists
    if (T.mesh().foundObject<IOdictionary>(radiationModel::typeName))
    {
        FatalErrorInFunction
            << "Radiation model already exists in registry, please check your setup!"
            << nl << "You have probably added a radiation model via fvOption"
            << " to a solver that already uses a radiation model."
            << exit(FatalError);
    }

    IOobject radIO
    (
        dictName,
        T.db().time().constant(),
        T.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    word modelType("none");
    if (radIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary(radIO).lookup("radiationModel") >> modelType;
    }
    else
    {
        Info<< "Radiation model not active: radiationProperties not found"
            << endl;
    }

    Info<< "Selecting radiationModel " << modelType << endl;

    const auto ctor = ctorTableLookup("radiationModel type", TConstructorTable_(), modelType);
    return autoPtr<radiationModel>(ctor(T));
}


Foam::autoPtr<Foam::radiation::radiationModel>
Foam::radiation::radiationModel::New
(
    const dictionary& dict,
    const volScalarField& T
)
{
    const word modelType(dict.lookup("radiationModel"));

    Info<< "Selecting radiationModel " << modelType << endl;

    const auto ctor = ctorTableLookup("radiationModel type", dictionaryConstructorTable_(), modelType);
    return autoPtr<radiationModel>(ctor(dict, T));
}


Foam::radiation::radiationModel& Foam::radiation::radiationModel::lookupOrCreate
(
    const volScalarField& T
)
{
    if (!T.db().foundObject<radiationModel>(dictName))
    {
        if (radiationModel::debug)
        {
            InfoInFunction
                << "constructing radiation model " << dictName
                << " for region " << T.db().name() << endl;
        }
        regIOobject::store(New(T).ptr());
    }

    return T.db().lookupObjectRef<radiationModel>(dictName);
}


// ************************************************************************* //
