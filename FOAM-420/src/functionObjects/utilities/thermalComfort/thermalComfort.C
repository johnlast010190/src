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
    (c) 2013 OpenFOAM Foundation
    (c) 2017, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "thermalComfort/thermalComfort.H"
#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "global/constants/constants.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "radiationModels/domSolar/domSolar.H"
#include "radiationModels/domSolar/radiantField.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "thermalComfort/comfortFunctions.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(thermalComfort, 0);
    addToRunTimeSelectionTable(functionObject, thermalComfort, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::thermalComfort::thermalComfort
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    calculateUTCI_(dict.lookupOrDefault<Switch>("utci", true)),
    correctUTCI_(dict.lookupOrDefault<Switch>("correctUTCI", false)),
    calculateISO7730_(dict.lookupOrDefault<Switch>("iso7730", true)),
    calculateASHRAE55_(dict.lookupOrDefault<Switch>("ashrae55", false)),
    calculateASHRAE55PMV_(dict.lookupOrDefault<Switch>("ashrae55PMV", true)),
    calculateWetBulbTemp_
    (
        dict.lookupOrDefault<Switch>("wetBulbTemperature", false)
    ),
    calculateISO7933_(dict.lookupOrDefault<Switch>("iso7933", false)),
    calculateISO7243_(dict.lookupOrDefault<Switch>("iso7243", false))
{
    thermalComfort::read(dict);
    thermalComfort::createFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::thermalComfort::~thermalComfort()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::functionObjects::thermalComfort::pRef()
{
    if (!pRef_.valid())
    {
        const surfaceScalarField& phi =
            lookupObject<surfaceScalarField>("phi");
        if (phi.dimensions() == dimVelocity*dimArea)
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
             );
            pRef_.reset
            (
                new dimensionedScalar
                (
                    "pRef",
                    transportProperties.lookup("pRef")
                )
            );
        }
        else
        {
            const fluidThermo& thermo =
                mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
            pRef_.reset
            (
                new dimensionedScalar("pRef", dimPressure, thermo.pRefValue())
            );
        }
    }
    return pRef_();
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::thermalComfort::getCp() const
{
    if (mesh_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo
            (mesh_.lookupObject<fluidThermo>(basicThermo::dictName));

        tmp<volScalarField> Cp(thermo.Cp());

        return Cp;
    }
    else if (mesh_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            mesh_.lookupObject<transportModel>("transportProperties");

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
Foam::functionObjects::thermalComfort::getNu() const
{
    if (mesh_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo
            (mesh_.lookupObject<fluidThermo>(basicThermo::dictName));

        tmp<volScalarField> nu(thermo.nu());

        return nu;
    }
    else if (mesh_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            mesh_.lookupObject<transportModel>("transportProperties");

        tmp<volScalarField> nu(transport.nu());

        return nu;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get nu()" << nl
            << "Supported objects are: basicThermo and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::thermalComfort::getKappa() const
{
    if (mesh_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo
            (mesh_.lookupObject<fluidThermo>(basicThermo::dictName));

        tmp<volScalarField> kappa(thermo.kappa());

        return kappa;
    }
    else if (mesh_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            mesh_.lookupObject<transportModel>("transportProperties");

        tmp<volScalarField> kappa(transport.lambda());

        return kappa;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get kappa()" << nl
            << "Supported objects are: basicThermo and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::thermalComfort::getRho() const
{
    if (mesh_.foundObject<volScalarField>("rho"))
    {
        tmp<volScalarField> rho(mesh_.lookupObject<volScalarField>("rho"));

        return rho;
    }
    else if (mesh_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            mesh_.lookupObject<transportModel>("transportProperties");

       tmp<volScalarField> rho(transport.rho());

       return rho;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get rho()" << nl
            << "Supported objects are: basicThermo  and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


void Foam::functionObjects::thermalComfort::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    const scalar initValue
)
{
    if (!foundObject<volScalarField>(fieldName))
    {
        mesh_.objectRegistry::store
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(fieldName, dims, initValue),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}


void Foam::functionObjects::thermalComfort::createFields()
{
    // Create relative humidity field
    storeField("relativeHumidity", dimless, RH_);

    // Create radiant Temperature-Field
    storeField("Tmrt", dimTemperature);

    if (calculateISO7730_)
    {
        storeField("DR", dimensionSet(0, 1, -1, 0, 0, 0, 0));
        storeField("TU");
        storeField("PMV");
        storeField("PPD");
        storeField("TOp", dimTemperature);
    }

    if (calculateASHRAE55_)
    {
        // Create field for ASHRAE-55  equivalent temperature
        storeField("SET", dimTemperature);

        if (calculateASHRAE55PMV_)
        {
            // Create field for ASHRAE-55 PMV
            storeField("PMVASHRAE");

            // Create field for ASHRAE-55 PPD
            storeField("PPDASHRAE");
        }
    }

    if (calculateWetBulbTemp_ || calculateISO7243_)
    {
        storeField("wetBulbTemperature", dimTemperature);
    }

    if (calculateUTCI_)
    {
        storeField("UTCI", dimTemperature);

        if (!foundObject<volScalarField>("wallDistCorrect"))
        {
            // wall distance field used for correcting local velocity
            // to height of 10m
            storeField("wallDistCorrect", dimless, 1.0);
            volScalarField& wDistCorrect =
                mesh_.lookupObjectRef<volScalarField>("wallDistCorrect");

            if (correctUTCI_)
            {
                wallDist wDist(mesh_, true);
                const volScalarField& wdist = wDist.y();
                const scalar a1 = log10(10/0.01);
                scalarField corrWallDict(wdist.primitiveField());
                forAll(corrWallDict, celli)
                {
                    if (wdist[celli] > 0.1)
                    {
                        corrWallDict[celli] = (a1/log10(wdist[celli]/0.01));
                    }
                    else
                    {
                        corrWallDict[celli] = 3.0;
                    }
                }
                wDistCorrect.primitiveFieldRef() = corrWallDict;
            }
        }
    }

    if (calculateISO7933_)
    {
        // Create final rectal temperature field
        storeField("Tre", dimTemperature);

        // Create total water loss field
        storeField("SWtotg", dimMass);

        // Create 'time when limit for rectal temperature is reached' field
        storeField("DlimTre", dimTime);

        // Create 'time when limit for water loss Dmax50 is reached' field
        // Dmax50: 7.5 percent of body mass of an average person
        storeField("DlimLoss50", dimTime);

        // Create 'time when limit for water loss Dmax95 is reached' field
        // Dmax95: 5 percent of body mass of 95 percent of the working people
        storeField("DlimLoss95", dimTime);
    }

    if (calculateISO7243_)
    {
        // Create wet bulb global temperature field
        storeField("WBGT", dimTemperature);
    }
}


bool Foam::functionObjects::thermalComfort::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << thermalComfort::type() << " "
        << thermalComfort::name() <<  " read:" << nl;

    clo_ = dict.lookupOrDefault<scalar>("clo", 1.0);

    if (clo_ > 1.5)
    {
        WarningInFunction
            << "Clothing insulation : " << clo_ << " exceeds maximum of 1.5 "
            << " recommended by ASHRAE-55-2017"
            <<endl;
    }

    if ((clo_ < 0.1) || (clo_ > 1.0))
    {
        WarningInFunction
            << "Clothing insulation : " << clo_ << " exceeds limits "
            << " of 0.1 (min) or 1.0 (max) recommended by ISO7933-2004"
            <<endl;
    }

    //convert clo to m2K/W
    clo_ *= 0.155;

    met_ = dict.lookupOrDefault<scalar>("met", 1.0);

    if (met_ > 2)
    {
        WarningInFunction
            << "Metabablic rate : " << met_ << " exceeds maximum of 2 "
            << "recommended by ASHRAE-55-2017"
            <<endl;
    }

    if ((met_ < 1.7) || (met_ > 7.7))
    {
        WarningInFunction
            << "Metabolic rate : " << met_ << " exceeds limits "
            << "of 1.7 (min) or 7.7 (max) recommended by ISO7933-2004"
            <<endl;
    }

    work_ = dict.lookupOrDefault<scalar>("work", 0.0);
    //convert mets to w/m2
    met_ *= 58.15;
    work_ *= 58.15;

    // PHS Model parameters
    weight_ = dict.lookupOrDefault<scalar>("weight", 75.0);
    height_ = dict.lookupOrDefault<scalar>("height", 1.8);
    drink_ =  dict.lookupOrDefault<Switch>("drink", true);
    duration_ = dict.lookupOrDefault<scalar>("duration", 480.0);
    posture_ = dict.lookupOrDefault<scalar>("posture", 2);
    accl_ =  dict.lookupOrDefault<Switch>("accl", true);
    thetaW_ = dict.lookupOrDefault<scalar>("thetaW", 0.0);
    walkSpeed_ = dict.lookupOrDefault<scalar>("walkSpeed", 0.0);

    // misc defaults
    occupantControl_ =
        dict.lookupOrDefault<Switch>("occupantControl", true);

    RH_ = dict.lookupOrDefault<scalar>("RH", 100.0);

    Tu_ = dict.lookupOrDefault<scalar>("Tu", 40.0);

    calculateTu_ = dict.lookupOrDefault<Switch>("calculateTu", false);

    upperLimitTu_ = dict.lookupOrDefault<scalar>("upperLimitTu", 100.0);
    lowerLimitTu_ = dict.lookupOrDefault<scalar>("lowerLimitTu", 0.0001);

    UName_ = dict.lookupOrDefault<word>("UName", "U");
    pName_ = dict.lookupOrDefault<word>("pName", "p");

    if (dict.found("solarSources"))
    {
        const wordList sources = wordList(dict.lookup("solarSources"));

        forAll(sources, sI)
        {
            solarFields_.append("I" + sources[sI]);
        }
    }
    else
    {
        solarFields_ = wordList(1, "Isolar");
    }

    return true;
}


bool Foam::functionObjects::thermalComfort::execute()
{
    // Do nothing - only valid on write
    return true;
}


bool Foam::functionObjects::thermalComfort::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    const volScalarField& T = lookupObject<volScalarField>("T");
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volScalarField& k = lookupObject<volScalarField>("k");
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>("phi");

    comfortFunctions comfort(clo_, met_, work_, Tu_);

    volScalarField& relativeHumidity =
        lookupObjectRef<volScalarField>("relativeHumidity");

    if (foundObject<volScalarField>("w"))
    {
        const volScalarField& w = lookupObject<volScalarField>("w");
        dimensionedScalar Mvap("Mwv", dimensionSet(1, 0, 0, 0, 1), 18.02e-3);
        // Check "Mwv" - shouldn't it be Mair?
        dimensionedScalar Mair("Mwv", dimensionSet(1, 0, 0, 0, 1), 28.96e-3);

        volScalarField Psat
        (
            IOobject
            (
                "Psat",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            ),
            mesh_,
            dimensionedScalar("Psat", dimPressure, 0.0)
        );

        Psat.primitiveFieldRef() =
            2337*Foam::exp(6879*(1/293.15 - 1./T.primitiveField())
          - 5.031*Foam::log(max(T.primitiveField(), SMALL)/293.15));

        forAll(Psat.boundaryField(), patchi)
        {
            Psat.boundaryFieldRef()[patchi] =
                2337*Foam::exp
                (
                    6879*(1/293.15 - 1/T.boundaryField()[patchi])
                  - 5.031*Foam::log(T.boundaryField()[patchi]/293.15)
                );
        }

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
             );
            dimensionedScalar rho("rho", transportProperties.lookup("rho"));

            relativeHumidity =
                100*w
               *(
                    max
                    (
                        pRef() + p*rho,
                        dimensionedScalar("pRef", pRef().dimensions(), 1e4)
                    ))*Mair/(Psat*(Mvap + w*(Mair - Mvap))
                );
        }
        else
        {
            relativeHumidity =
                100*w
               *max
                (
                    pRef() + p,
                    dimensionedScalar("p", p.dimensions(), 1e4)
                )*Mair/(Psat*(Mvap + w*(Mair - Mvap)));
        }

        relativeHumidity.boundaryFieldRef().evaluate();
        relativeHumidity.write();
    }

    volScalarField Qrs
    (
        IOobject
        (
            "Qrs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qrs", dimMass/pow3(dimTime), 0.0)
    );

    const radiation::radiationModel& radiation =
        mesh_.lookupObject<radiation::radiationModel>("radiationProperties");

    volScalarField& radiantTemperature =
        lookupObjectRef<volScalarField>("Tmrt");

    if (isA<radiation::fvDOM>(radiation))
    {
        const radiation::fvDOM& dom =
            const_cast<radiation::fvDOM&>
            (
                refCast<const radiation::fvDOM>(radiation)
            );

        for (int i = 0; i < dom.nRay(); i++)
        {
            Qrs += dom.IRay(i).cellQr()();
        }

        //Setup for solar radiation
        forAll(solarFields_, sFI)
        {
            word solarField = solarFields_[sFI];

            if (foundObject<volScalarField>(solarField))
            {
                const volScalarField& Isolar =
                    lookupObject<volScalarField>(solarField);

                forAll(Isolar.primitiveField(), celli)
                {
                    Qrs[celli] += 0.25*Isolar[celli];
                }
            }
        }
        radiantTemperature =
            Foam::pow((mag(Qrs)/constant::physicoChemical::sigma), 0.25);
        radiantTemperature.boundaryFieldRef().evaluate();
    }
    else
    {
        const fvPatchList& patches = mesh_.boundary();
        DynamicList<boundBox> patchBB(patches.size());
        DynamicList<scalar> aveWallPatchTemperature(patches.size());
        DynamicList<scalar> wallPatchAreas(patches.size());

        forAll(patches, patchi)
        {
            const fvPatch& fvp = patches[patchi];
            if (isType<wallFvPatch>(fvp))
            {
                const polyPatch& pp = fvp.patch();
                scalar Apatch  = gSum(fvp.magSf());
                if (Apatch > SMALL)
                {
                    scalar Pavg =
                        gSum(T.boundaryField()[patchi]*fvp.magSf())/Apatch;

                    const pointField& lpoints = pp.localPoints();
                    boundBox bb(lpoints, true);

                    patchBB.append(bb);
                    aveWallPatchTemperature.append(Pavg);
                    wallPatchAreas.append(Apatch);
                }
            }
        }

        forAll(radiantTemperature, celli)
        {
            if (aveWallPatchTemperature.size() > 0)
            {
                scalar sumTemp(0);
                scalar sumWeights(0);
                point cc = mesh_.cellCentres()[celli];

                forAll(aveWallPatchTemperature, i)
                {
                    point nearest = patchBB[i].nearest(cc);
                    scalar pDist = mag(cc-nearest);
                    scalar weight = wallPatchAreas[i]/(pDist + SMALL);
                    sumTemp += (weight * aveWallPatchTemperature[i]);
                    sumWeights += weight;
                }
                radiantTemperature[celli] = sumTemp/sumWeights;
            }
            else
            {
                radiantTemperature[celli] = T[celli];
            }
        }
    }

    radiantTemperature.write();

    if (calculateUTCI_)
    {
        volScalarField& UTCI = lookupObjectRef<volScalarField>("UTCI");
        volScalarField& wallDistCorrect =
            lookupObjectRef<volScalarField>("wallDistCorrect");

        forAll(UTCI, celli)
        {
            scalar magCellU = mag(U[celli]);

            if (correctUTCI_)
            {
                magCellU *= wallDistCorrect[celli];
            }

            UTCI[celli] =
                comfort.UTCI
                (
                    magCellU,
                    T[celli],
                    radiantTemperature[celli],
                    relativeHumidity[celli]
                );
        }
        UTCI.boundaryFieldRef().evaluate();
        UTCI.write();
    }


    if (calculateISO7730_)
    {
        volScalarField& PMV = lookupObjectRef<volScalarField>("PMV");
        volScalarField& PPD = lookupObjectRef<volScalarField>("PPD");
        volScalarField& TOp = lookupObjectRef<volScalarField>("TOp");
        volScalarField& DR = lookupObjectRef<volScalarField>("DR");
        volScalarField& TU = lookupObjectRef<volScalarField>("TU");

        forAll(PMV, celli)
        {
            scalar cellT = T[celli];
            scalar cellRT = radiantTemperature[celli];
            scalar cellRH = relativeHumidity[celli];
            vector cellU = U[celli];
            scalar cellK = k[celli];

            scalar va = mag(cellU);

            if (calculateTu_)
            {
                TU[celli] =
                    comfort.TU(va, cellK, upperLimitTu_, lowerLimitTu_);
                Tu_ = TU[celli];
            }
            else
            {
                TU[celli] = Tu_;
            }

            DR[celli] = comfort.DR(va, cellT, Tu_);
            PMV[celli] = comfort.PMV(va, cellT, cellRT, cellRH);
            PPD[celli] = comfort.PPD(PMV[celli]);
            TOp[celli] = comfort.TOp(va, cellT, cellRT);
        }

        if (calculateTu_)
        {
            TU.boundaryFieldRef().evaluate();
            TU.write();
        }

        DR.boundaryFieldRef().evaluate();
        PMV.boundaryFieldRef().evaluate();
        PPD.boundaryFieldRef().evaluate();
        TOp.boundaryFieldRef().evaluate();

        DR.write();
        PMV.write();
        PPD.write();
        TOp.write();
    }

    if (calculateASHRAE55_)
    {
        volScalarField& SET = lookupObjectRef<volScalarField>("SET");

        scalar rho = 1.0;

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
             );
            dimensionedScalar transRho(transportProperties.lookup("rho"));

            rho = transRho.value();
        }

        if (calculateASHRAE55PMV_)
        {
            volScalarField& PMV = lookupObjectRef<volScalarField>("PMVASHRAE");
            volScalarField& PPD = lookupObjectRef<volScalarField>("PPDASHRAE");
            const scalar pMin = 1e4;
            forAll(PMV, celli)
            {
                //Ensure pressures are correctly scaled
                scalar cellP = pMin;
                if (p.dimensions() == dimPressure)
                {
                    cellP = max(p[celli] + pRef().value(), pMin);
                }
                else
                {
                    cellP = max(p[celli]*rho + pRef().value(), pMin);
                }
                scalar cellT = T[celli];
                scalar cellRT = radiantTemperature[celli];
                const scalar cellRH = relativeHumidity[celli];
                scalar cellMagU = mag(U[celli]);

                SET[celli] =
                    comfort.correctSET
                    (
                        cellP,
                        cellRH,
                        cellMagU,
                        cellT,
                        cellRT,
                        occupantControl_
                    );

                PMV[celli] = comfort.PMV(cellMagU, cellT, cellRT, cellRH);
                PPD[celli] = comfort.PPD(PMV[celli]);
            }

            SET.boundaryFieldRef().evaluate();
            SET.write();

            PMV.boundaryFieldRef().evaluate();
            PMV.write();

            PPD.boundaryFieldRef().evaluate();
            PPD.write();
        }
        // if PMV and PPD (ASHRAE) not necessary, only calculate SET
        else
        {
            const scalar pMin = 1e4;
            forAll(SET, celli)
            {
                //Ensure pressures are correctly scaled
                scalar pAbsolute = pMin;
                if (p.dimensions() == dimPressure)
                {
                    pAbsolute = max(p[celli] + pRef().value(), pMin);
                }
                else
                {
                    pAbsolute = max(p[celli]*rho + pRef().value(), pMin);
                }

                // Using uncorrected SET function
                SET[celli] =
                    comfort.SET
                    (
                        mag(U[celli]),
                        T[celli],
                        pAbsolute,
                        radiantTemperature[celli],
                        relativeHumidity[celli]
                    );
            }

            SET.boundaryFieldRef().evaluate();
            SET.write();
        }
    }

    if (calculateWetBulbTemp_ && !calculateISO7243_)
    {
        volScalarField& wetBulb =
            lookupObjectRef<volScalarField>("wetBulbTemperature");

        forAll(wetBulb, celli)
        {
            wetBulb[celli] =
                comfort.wetBulb(T[celli], relativeHumidity[celli]);
        }
        wetBulb.boundaryFieldRef().evaluate();
        wetBulb.write();
    }

    if (calculateISO7933_)
    {
        volScalarField& Tre = lookupObjectRef<volScalarField>("Tre");
        volScalarField& SWtotg = lookupObjectRef<volScalarField>("SWtotg");
        volScalarField& DlimTre = lookupObjectRef<volScalarField>("DlimTre");
        volScalarField& DlimLoss50 =
            lookupObjectRef<volScalarField>("DlimLoss50");
        volScalarField& DlimLoss95 =
            lookupObjectRef<volScalarField>("DlimLoss95");

        forAll(Tre, celli)
        {
            scalar cellT = T[celli];
            scalar cellRT = radiantTemperature[celli];
            scalar cellRH = relativeHumidity[celli];
            scalar va = mag(U[celli]);

            scalarField iso7933Metrics =
                comfort.PHSiso7933
                (
                    va,
                    cellT,
                    cellRT,
                    cellRH,
                    weight_,
                    height_,
                    drink_,
                    duration_,
                    posture_,
                    accl_,
                    thetaW_,
                    walkSpeed_
                );

            Tre[celli] = iso7933Metrics[0];
            SWtotg[celli] = iso7933Metrics[1];
            DlimTre[celli] = iso7933Metrics[2];
            DlimLoss50[celli] = iso7933Metrics[3];
            DlimLoss95[celli] = iso7933Metrics[4];
        }

        Tre.boundaryFieldRef().evaluate();
        Tre.write();

        SWtotg.boundaryFieldRef().evaluate();
        SWtotg.write();

        DlimTre.boundaryFieldRef().evaluate();
        DlimTre.write();

        DlimLoss50.boundaryFieldRef().evaluate();
        DlimLoss50.write();

        DlimLoss95.boundaryFieldRef().evaluate();
        DlimLoss95.write();
    }

    if (calculateISO7243_)
    {
        volScalarField& WBGT = lookupObjectRef<volScalarField>("WBGT");
        volScalarField& wetBulb =
            lookupObjectRef<volScalarField>("wetBulbTemperature");

        scalar rhoRef = 1.0;
        if (phi.dimensions() == dimVelocity*dimArea)
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            rhoRef =
                dimensionedScalar(transportProperties.lookup("rho")).value();
        }

        // Access to properties values
        const volScalarField kappa(getKappa());
        const volScalarField rho(getRho());
        const volScalarField nu(getNu());
        const volScalarField Cp(getCp());

        // Radiation Parameters

        // Flag to defined indoor/outdoor condition
        Switch flagSolarRad = false;

        // Solar Heat Flux field
        volScalarField solarFlux
        (
            IOobject
            (
                "solarFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("solarFlux", dimMass/pow3(dimTime), 0.0)
        );

        // Need to check if more than one solar source is defined
        forAll(solarFields_, sFI)
        {
            word solarField = solarFields_[sFI];

            if (foundObject<volScalarField>(solarField))
            {
                flagSolarRad = true;

                const volScalarField& Isolar =
                    lookupObject<volScalarField>(solarField);

                forAll(solarFlux, celli)
                {
                    solarFlux[celli] += Isolar[celli];
                }
            }
        }

        scalar zenithAngle = 0.0;

        if (flagSolarRad)
        {
            const word modelType
            (
                IOdictionary
                (
                    IOobject
                    (
                        "radiationProperties",
                        this->mesh_.time().constant(),
                        this->mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                ).lookup("radiationModel")
            );

            if (modelType != "fvDOM")
            {
                FatalErrorInFunction
                    << "fvDOM radiation model is the only one supported "
                    << " with ISO7243 thermal comfort standard "
                    << exit(FatalError);
            }

            // Access to solarCalculator/userDefined zenith angle
            const radiation::radiationModel& domRadObj =
                lookupObject<radiation::radiationModel>("radiationProperties");

            if (isA<radiation::fvDOM>(domRadObj))
            {
                const radiation::fvDOM& fvDOMobj =
                    refCast<const radiation::fvDOM>(domRadObj);

                zenithAngle =
                    fvDOMobj.getDomSolarObj()[0].getZenithAngle();
            }
        }

        const scalar pMin = 1e4;
        forAll(WBGT, celli)
        {
            scalar pAbsolute = pMin;
            if (p.dimensions() == dimPressure)
            {
                pAbsolute = max(p[celli] + pRef().value(), pMin);
            }
            else
            {
                pAbsolute = max(p[celli]*rhoRef + pRef().value(), pMin);
            }

            wetBulb[celli] =
                comfort.wetBulb(T[celli], relativeHumidity[celli]);

            WBGT[celli] =
                comfort.WBGTiso7243
                (
                    pAbsolute,
                    mag(U[celli]),
                    T[celli],
                    relativeHumidity[celli],
                    wetBulb[celli],
                    rho[celli],
                    nu[celli],
                    Cp[celli],
                    kappa[celli],
                    solarFlux[celli],
                    flagSolarRad,
                    zenithAngle
                );
        }

        WBGT.boundaryFieldRef().evaluate();
        WBGT.write();

        wetBulb.boundaryFieldRef().evaluate();
        wetBulb.write();
    }

    return true;
}


// ************************************************************************* //
