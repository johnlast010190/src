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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "objectiveFunctionObject/objectiveFunctionObject.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/fvPatchFields/derived/resistivePressure/resistivePressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(objectiveFunctionObject, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::objectiveFunctionObject::setInletPatches(const dictionary& dict)
{
    if (dict.found("inletPatches"))
    {
        inletPatchPtr_.reset(new boolList(mesh_.boundary().size(), false));
        boolList& inletPatch = inletPatchPtr_();

        const dictionary& inletPatches =
            dict.subDict("inletPatches");
        if (inletPatches.found("exactNamed"))
        {
            wordList namedPatches(inletPatches.lookup("exactNamed"));
            forAll(namedPatches, nI)
            {
                label patchI =
                    mesh_.boundaryMesh().findPatchID(namedPatches[nI]);

                if
                (
                    (patchI != -1)
                )
                {
                    inletPatch[patchI] = true;
                }
            }
        }
        if (inletPatches.found("partialNamed"))
        {
            wordList partialNamePatches(inletPatches.lookup("partialNamed"));
            const fvPatchList& patches = mesh_.boundary();

            forAll(partialNamePatches, nI)
            {
                word substring = partialNamePatches[nI];

                forAll(patches, pI)
                {
                    word name = patches[pI].name();

                    if
                    (
                        name.find(substring, 0) != string::npos
                    )
                    {
                        inletPatch[pI] = true;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::objectiveFunctionObject::P() const
{
    if (compressible())
    {
        return {p()};
    }
    else
    {
        return tmp<volScalarField>(rho()*p());
    }
}


Foam::tmp<Foam::volScalarField>
Foam::objectiveFunctionObject::rho() const
{
    if (foundObject<volScalarField>(rhoName_))
    {
        return (lookupObject<volScalarField>(rhoName_));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return laminarT.rho();
    }
    else
    {
        FatalErrorInFunction
            << "Couldn't find a valid rho field" << exit(FatalError);

        return
            tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rho",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::objectiveFunctionObject::mu() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.mu();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return turb.mu()*rho();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu", transportProperties.lookup("nu"));

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);
        ::abort();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::objectiveFunctionObject::muEff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.muEff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return turb.muEff()*rho();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu", transportProperties.lookup("nu"));

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);
        ::abort();
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::objectiveFunctionObject::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U())));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U())));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu", transportProperties.lookup("nu"));

        return -rho()*nu*dev(twoSymm(fvc::grad(U())));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);
        ::abort();
    }
}


Foam::boolList
Foam::objectiveFunctionObject::getObjectivePatches
(
    const dictionary& objectiveDict
) const
{
    boolList objectivePatch(mesh_.boundary().size(), false);
    label nObjectivePatches = 0;

    if (objectiveDict.found("objectivePatches"))
    {
        const dictionary& patchDict = objectiveDict.subDict("objectivePatches");
        objectivePatch = false;

        if (patchDict.found("exactNamed"))
        {
            wordList namedPatches(patchDict.lookup("exactNamed"));

            forAll(namedPatches, nI)
            {
                label patchI =
                    mesh_.boundaryMesh().findPatchID(namedPatches[nI]);

                if (patchI != -1)
                {
                    objectivePatch[patchI] = true;
                    nObjectivePatches++;
                }
            }
        }

        if (patchDict.found("partialNamed"))
        {
            wordList partialNamePatches(patchDict.lookup("partialNamed"));
            const fvPatchList& patches = mesh_.boundary();

            forAll(partialNamePatches, nI)
            {
                word substring = partialNamePatches[nI];

                forAll(patches, pI)
                {
                    word name = patches[pI].name();

                    if (name.find(substring, 0) != string::npos)
                    {
                        objectivePatch[pI] = true;
                        nObjectivePatches++;
                    }
                }
            }
        }
    }
    else
    {
        objectivePatch = true;
        nObjectivePatches = objectivePatch.size();
    }

    if (nObjectivePatches == 0)
    {
        FatalErrorInFunction
            <<"No objective patches set for objective " << name()
            << ". Check `objectivePatches` syntax." << exit(FatalError);
    }

    return objectivePatch;
}


Foam::scalar
Foam::objectiveFunctionObject::objectivePatchInletArea() const
{
    scalar area = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (p.type() == "inlet" || p.patch().physicalType() == "inlet")
         && objectivePatch_[patchI]
        )
        {
            area += sum(p.magSf());
        }
    }

    reduce(area, sumOp<scalar>());

    return area;
}


Foam::scalar
Foam::objectiveFunctionObject::objectivePatchOutletArea() const
{
    scalar area = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (p.type() == "outlet" || p.patch().physicalType() == "outlet")
         && objectivePatch_[patchI]
        )
        {
            area += sum(p.magSf());
        }
    }

    reduce(area, sumOp<scalar>());

    return area;
}


Foam::scalar
Foam::objectiveFunctionObject::inletFlowRate() const
{
    scalar ifr = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if (!inletPatchPtr_.valid())
        {
            if (p.type() == "inlet" || p.patch().physicalType() == "inlet")
            {
                ifr += sum(phi().boundaryField()[patchI]);
            }
        }
        else
        {
            const boolList& inletPatch = inletPatchPtr_();
            if (inletPatch[patchI])
            {
                ifr += sum(phi().boundaryField()[patchI]);
            }
        }
    }

    reduce(ifr, sumOp<scalar>());

    // to shore things up if the initial inlet velocity is negative or zero
    ifr = max(mag(ifr), SMALL);

    return ifr;
}


Foam::scalar
Foam::objectiveFunctionObject::calculatePowerLoss() const
{
    scalar powerLoss = 0;

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb =
        tdevRhoReff().boundaryField();

    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    forAll(mesh_.boundary(), patchI)
    {
        if (objectivePatch_[patchI])
        {
            const fvPatch& patch = mesh_.boundary()[patchI];

            const scalarField& rhop = rho()().boundaryField()[patchI];
            const vectorField& Vp = U().boundaryField()[patchI];

            scalarField sw( (Sfb[patchI] & devRhoReffb[patchI]) & Vp );

            scalarField pp( p().boundaryField()[patchI] );
            scalarField phip( phi().boundaryField()[patchI] );

            //grab outside pressure from resistive outlet
            if
            (
                isA<resistivePressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
            )
            {
                pp =
                    dynamic_cast<const resistivePressureFvPatchScalarField&>
                    (
                        p().boundaryField()[patchI]
                    ).p0();
            }

            //different formulation based on flow solver
            if (compressible())
            {
                if
                (
                    patch.type() == "inlet"
                 || patch.patch().physicalType() == "inlet"
                 || patch.type() == "outlet"
                 || patch.patch().physicalType() == "outlet"
                )
                {
                    scalar tmpPower = 0;

                    tmpPower = sum(phip*(pp/rhop + 0.5*magSqr(Vp)));

                    powerLoss -= tmpPower;
                }
            }
            else
            {
                if
                (
                    patch.type() == "inlet"
                 || patch.patch().physicalType() == "inlet"
                 || patch.type() == "outlet"
                 || patch.patch().physicalType() == "outlet"
                )
                {
                    scalar tmpPower = sum
                    (
                        rhop*phip*(pp + 0.5*magSqr(Vp))
                    );

                    powerLoss -= tmpPower;
                }
            }
        }
    }

    reduce(powerLoss, sumOp<scalar>());

    return powerLoss;
}


void Foam::objectiveFunctionObject::writeFileHeader
(
    Ostream& os
)
{
    writeCommented(os, "Time");
    writeDelimited(os, "objectiveValue");
    writeDelimited(os, "normalizedWeighting");
    os << endl;
}


void Foam::objectiveFunctionObject::logToFile()
{
    if (logToFile_ && (Pstream::master() || !Pstream::parRun()))
    {
        fileName logFile;

        if (Pstream::parRun())
        {
            logFile =
                mesh_.time().path()/".."/"postProcessing"/"adjoint"
               /mesh_.time().timeName();
        }
        else
        {
            logFile =
                mesh_.time().path()/"postProcessing"/"adjoint"
               /mesh_.time().timeName();
        }

        mkDir(logFile);

        fileName objFileName
        (
            word(name()) + (".dat")
        );

        objFilePtr_.reset
        (
            new OFstream(logFile/objFileName)
        );
    }
}


void Foam::objectiveFunctionObject::createFiles
(
    bool useAdjointFileFormat
)
{
    if (useAdjointFileFormat)
    {
        logToFile();
    }
    else
    {
        objFilePtr_ = createFile(this->type());
    }

    if (objFilePtr_.valid())
    {
        writeFileHeader(objFilePtr_());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::objectiveFunctionObject::objectiveFunctionObject
(
    const Time& time,
    dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    fvMeshFunctionObject
    (
        objectiveDict.lookup<word>("objectiveName"),
        time,
        objectiveDict
    ),
    writeFile(mesh_, objectiveDict.lookup<word>("objectiveName")),
    UName_
    (
        objectiveDict.lookupOrDefault<word>("UName", "U")
    ),
    pName_
    (
        objectiveDict.lookupOrDefault<word>("pName", "p")
    ),
    phiName_
    (
        objectiveDict.lookupOrDefault<word>("phiName", "phi")
    ),
    rhoName_
    (
        objectiveDict.lookupOrDefault<word>("rhoName", "rho")
    ),
    psiName_
    (
        objectiveDict.lookupOrDefault<word>("psiName", "psi")
    ),
    TName_
    (
        objectiveDict.lookupOrDefault<word>("TName", "T")
    ),
    objectiveValue_(0.0),
    objectivePatch_
    (
        getObjectivePatches(objectiveDict)
    ),
    baseGeometryName_("geometry"),
    zoneNames_
    (
        objectiveDict.lookupOrDefault<wordList>("volumeZonesToOptimize",wordList(0))
    ),
    inletPatchPtr_(),
    logToFile_
    (
        objectiveDict.lookupOrDefault<bool>("logToFile", true)
    ),
    objFilePtr_()
{
    setInletPatches(objectiveDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectiveFunctionObject::~objectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::objectiveFunctionObject::compressible() const
{
    if (phi().dimensions() == dimVelocity*dimArea)
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool
Foam::objectiveFunctionObject::read
(
    const dictionary& objectiveDict
)
{
    fvMeshFunctionObject::read(objectiveDict);

    objectivePatch_ = getObjectivePatches(objectiveDict);

    zoneNames_ = objectiveDict.lookupOrDefault<wordList>
    (
        "volumeZonesToOptimize",
        wordList(0)
    );
    return true;
}

void
Foam::objectiveFunctionObject::writeOutputToFile()
{
    if (writeToFile() && Pstream::master())
    {
        writeTime(objFilePtr_());
        objFilePtr_()
        << tab << objectiveValue_
        << endl;
    }
}

bool Foam::objectiveFunctionObject::writeOutputToAdjointFile() const
{
    return true;
}

void
Foam::objectiveFunctionObject::setBaseGeometryName(const word& geometryName)
{
    baseGeometryName_ = geometryName;
}


// ************************************************************************* //
