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
    (c) 2010-2019, Esi Ltd

\*---------------------------------------------------------------------------*/

#include "zoneForces/zoneForces.H"
#include "cfdTools/general/porosityModel/porosityModel/IOporosityModelList.H"

#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(zoneForces, 0);

    addToRunTimeSelectionTable(functionObject, zoneForces, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::zoneForces::zoneForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    nAveragingSteps_(1),
    verbose_(false),
    logToFile_(false),
    porousZoneNames_(wordList(0)),
    averagingIndex_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::zoneForces::~zoneForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::functionObjects::zoneForces::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
             lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if
    (
        foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        return getRho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu", transportProperties.lookup("nu"));

        return getRho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);
        ::abort();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::zoneForces::getRho() const
{
    if (foundObject<volScalarField>("rho"))
    {
        return(lookupObject<volScalarField>("rho"));
    }
    else
    {
        dimensionedScalar rho
        (
            "rho",
            lookupObject<dictionary>("transportProperties")
                .lookup("rho")
        );

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                rho
            )
        );
    }
}

bool Foam::functionObjects::zoneForces::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    //calculate
    averagingIndex_++;
    calculate();

    const HashTable<const porosityModel*> models =
        obr_.lookupClass<porosityModel>();

    if (averagingIndex_ == nAveragingSteps_)
    {
        averagingIndex_ = 0;
    }

    //write
    if (logToFile_ && (Pstream::master() || !Pstream::parRun()))
    {
        massFluxFilePtr_()
            << mesh_.time().timeName();
        pressureDropFilePtr_()
            << mesh_.time().timeName();

        forAll(porousZoneNames_, zoneI)
        {
            massFluxFilePtr_()
                << tab << massFlux_[zoneI];
            pressureDropFilePtr_()
                << tab << pressureDrop_[zoneI];
        }
        pressureDropFilePtr_() << endl;
        massFluxFilePtr_() << endl;
    }

    return true;
}

Foam::scalar Foam::functionObjects::zoneForces::zoneMassFlux
(
    const labelList& faces
)
{
    scalar fsum = 0;
    const surfaceScalarField& phi
         = mesh_.lookupObject<surfaceScalarField>("phi");

    surfaceScalarField rhof(fvc::interpolate(getRho()));

    bool incomp = false;
    if (phi.dimensions() == dimVelocity * dimArea)
    {
        incomp = true;
    }

    forAll(faces, i)
    {
       label meshFaceI = faces[i];
       if (meshFaceI < mesh_.nInternalFaces())
       {
           if (incomp)
           {
               fsum += rhof[meshFaceI]*mag(phi[meshFaceI]);
           }
           else
           {
               fsum += mag(phi[meshFaceI]);
           }
       }
       else
       {
          label patchI = mesh_.boundaryMesh().whichPatch(meshFaceI);
          const polyPatch& patch = mesh_.boundaryMesh()[patchI];
          if (!isA<emptyPolyPatch>(patch))
          {
              label bfI = meshFaceI - patch.start();

              if (incomp)
              {
                  fsum += rhof.boundaryField()[patchI][bfI]
                      *mag(phi.boundaryField()[patchI][bfI]);
              }
              else
              {
                   fsum += mag(phi.boundaryField()[patchI][bfI]);
              }
          }
       }
    }
    reduce(fsum, sumOp<scalar>());
    return fsum;
}


Foam::scalar Foam::functionObjects::zoneForces::zoneNormalForce
(
    const labelList& faces,
    const vector& force
)
{
    scalar sumForce = 0;

    const surfaceVectorField& Sf = mesh_.Sf();

    forAll(faces, i)
    {
        label meshFaceI = faces[i];

        if (meshFaceI < mesh_.nInternalFaces())
        {
            sumForce += mag(force & Sf[meshFaceI]);
        }
        else
        {
            label patchI = mesh_.boundaryMesh().whichPatch(meshFaceI);
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
           if (!isA<emptyPolyPatch>(patch))
           {
               sumForce += mag
               (
                   force & Sf.boundaryField()[patchI][meshFaceI - patch.start()]
               );
           }
        }
    }
    reduce(sumForce, sumOp<scalar>());
    if (sumForce > SMALL)
    {
        sumForce = 2.0*magSqr(force)/sumForce;
    }
    else
    {
        sumForce = 0;
    }

    return sumForce;
}

void Foam::functionObjects::zoneForces::calculate()
{
    const HashTable<const porosityModel*> models =
        obr_.lookupClass<porosityModel>();

    const volVectorField& U = lookupObject<volVectorField>("U");
    const volScalarField rho(this->getRho());
    const volScalarField mu(this->mu());

    //averaging coefficients
    scalar alpha = scalar(averagingIndex_ - 1.0)/scalar(averagingIndex_);
    scalar beta = 1.0/scalar(averagingIndex_);

    List<scalar> massFlux(massFlux_.size(), 0.);
    List<vector> force(force_.size(), vector::zero);
    List<scalar> pressureDrop(pressureDrop_.size(), 0.);

    vectorField porousForceField(mesh_.nCells(), vector::zero);
    vectorField porousMomentField(mesh_.nCells(), vector::zero);

    forAll(porousZoneNames_, zoneI)
    {
        porosityModel& pm =
            const_cast<porosityModel&>(*models[porousZoneNames_[zoneI]]);
        vectorField fPTot(pm.force(U, rho, mu));

        // calculate average force direction
        vector aveForce = gSum(fPTot);
        if (mag(aveForce) > SMALL)
        {
           aveForce /= mag(aveForce);
        }
        forAll(fPTot, cellI)
        {
           if ((fPTot[cellI] & aveForce) < 0.)
           {
              fPTot[cellI] = vector::zero;
           }
        }

        force[zoneI] = gSum(fPTot);
        force_[zoneI] = alpha*force_[zoneI] + beta*force[zoneI];

        pressureDrop[zoneI] = zoneNormalForce
            (boundaryFaces_[zoneI], force[zoneI]);
        pressureDrop_[zoneI] = alpha*pressureDrop_[zoneI]
            + beta*pressureDrop[zoneI];

        massFlux[zoneI] = 0.5 * zoneMassFlux(boundaryFaces_[zoneI]);
        massFlux_[zoneI] = alpha*massFlux_[zoneI] + beta*massFlux[zoneI];

        if (verbose_)
        {
            Info<< pm.name()
                 << " pressure drop [Pa]  = " << pressureDrop_[zoneI]
                 << " mass flux [kg/s] = " << massFlux_[zoneI]
                 << " force [N] = " << force_[zoneI]
                 << endl;
        }
    }
}


bool Foam::functionObjects::zoneForces::write()
{
    return true;
}


bool Foam::functionObjects::zoneForces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    const HashTable<const porosityModel*> models =
        obr_.lookupClass<porosityModel>();
    porousZoneNames_ = models.sortedToc();

    verbose_ = false;
    if (dict.found("verbose"))
    {
        verbose_= dict.lookup("verbose");
    }

    logToFile_= dict.lookupOrDefault<bool>("logToFile", true);

    averagingIndex_ = 0;
    nAveragingSteps_ = 1;
    if (dict.found("nAveragingSteps"))
    {
        nAveragingSteps_= readLabel(dict.lookup("nAveragingSteps"));
    }

    //file update
    if (logToFile_ && (Pstream::master() || !Pstream::parRun()))
    {
        massFluxFilePtr_ = createFile("porous-massFlux");
        pressureDropFilePtr_ = createFile("porous-pressureDrop");

        //add headers
        writeCommented(massFluxFilePtr_(), "Time");
        writeCommented(pressureDropFilePtr_(), "Time");

        forAll(porousZoneNames_, zoneI)
        {
            const porosityModel& pm = *models[porousZoneNames_[zoneI]];
            writeDelimited(massFluxFilePtr_(),pm.name());
            writeDelimited(pressureDropFilePtr_(),pm.name());
        }
        massFluxFilePtr_() <<endl;
        pressureDropFilePtr_() <<endl;
    }

    boundaryFaces_.setSize(models.size());

    forAll(porousZoneNames_, zoneI)
    {
        const porosityModel& pm = *models[porousZoneNames_[zoneI]];
        boundaryFaces_[zoneI] = pm.zoneBoundaryFaces();
    }
    massFlux_.setSize(models.size(), 0.0);
    pressureDrop_.setSize(models.size(), 0.0);
    force_.setSize(models.size(), vector::zero);

    return true;
}


// ************************************************************************* //
