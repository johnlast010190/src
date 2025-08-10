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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "acousticAnalogy/DavidsonCurle.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(DavidsonCurle, 0);
    addToRunTimeSelectionTable(functionObject, DavidsonCurle, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::DavidsonCurle::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "DavidsonCurle Acoustic Analogy");
    writeCommented(os, "Time");
    forAll(observers_, obsI)
    {
        writeDelimited(os,observers_[obsI].name());
    }
    os << endl;
}


void Foam::functionObjects::DavidsonCurle::writeDavidsonCurle()
{
    forAll(observers_, obsI)
    {
        Log << observers_[obsI].name() << "   "
             << observers_[obsI].pPrime() << endl;
    }

    // File output
    file() << obr_.time().timeName() << tab << setw(1) << "   ";
    forAll(observers_, obsI)
    {
        file() << observers_[obsI].pPrime() << "   ";
    }
    file() << endl;
}


void Foam::functionObjects::DavidsonCurle::initialise()
{
    if (initialised_)
    {
        return;
    }

    if
    (
        !obr_.foundObject<volVectorField>(UName_)
     || !obr_.foundObject<volScalarField>(pName_)
    )
    {
        FatalErrorInFunction
            << "Could not find " << UName_ << ", " << pName_
            << exit(FatalError);
    }

    if
    (
        cellZoneID_ == -1
     && cellZoneName_ != word::null
    )
    {
        FatalErrorInFunction
            << "Could not find cellZone '" << cellZoneName_
            << "' in database." << nl
            << exit(FatalError);
    }

    initialised_ = true;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::DavidsonCurle::p() const
{
    return
    (
        rhoRef_ * obr_.lookupObject<volScalarField>(pName_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::DavidsonCurle::dpdt() const
{
    return
    (
        rhoRef_ * Foam::fvc::ddt(obr_.lookupObject<volScalarField>(pName_))
    );
}


Foam::tmp<Foam::volTensorField> Foam::functionObjects::DavidsonCurle::Tij() const
{
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

    return
    (
        rhoRef_*(U*U)
    );
}


Foam::tmp<Foam::volTensorField> Foam::functionObjects::DavidsonCurle::dTijdt() const
{
    // Get velocity field of current and last two time steps
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volVectorField& U0 = lookupObject<volVectorField>(UName_).oldTime();
    const volVectorField& U00 =
        lookupObject<volVectorField>(UName_).oldTime().oldTime();

    // Get time stepping size
    scalar deltaT = mesh_.time().deltaTValue();
    scalar deltaT0;
    if (mesh_.time().timeIndex() < 2)
    {
        // Revert to Euler
        deltaT0 = GREAT;
    }
    else
    {
        deltaT0 = mesh_.time().deltaT0Value();
    }

    // Calculate coefficients
    scalar coefft   = 1.0 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    // Second-order backward-differencing time derivative
    // assuming constant time stepping
    // see: backwardDdtScheme.C
    return
    (
        rhoRef_*
        (
            (
                coefft*U*U
              - coefft0*U0*U0
              + coefft00*U00*U00
            )
          / mesh_.time().deltaT()
        )
    );
}


Foam::tmp<Foam::volTensorField> Foam::functionObjects::DavidsonCurle::d2Tijdt2() const
{
    // Get velocity field of current and last two time steps
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volVectorField& U0 = lookupObject<volVectorField>(UName_).oldTime();
    const volVectorField& U00 =
        lookupObject<volVectorField>(UName_).oldTime().oldTime();

    // Get time stepping size
    scalar deltaT = mesh_.time().deltaTValue();
    scalar deltaT0;
    if (mesh_.time().timeIndex() < 2)
    {
        // For start-up, assume deltaT0 = deltaT
        deltaT0 = deltaT;
    }
    else
    {
        deltaT0 = mesh_.time().deltaT0Value();
    }
    scalar dt = 4.0 / sqr(deltaT + deltaT0);

    // Calculate coefficients
    scalar coefft   = (deltaT + deltaT0)/(2.*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2.*deltaT0);
    scalar coefft0  = coefft + coefft00;

    // First-order Euler second time derivative
    // see: EulerD2dt2Scheme.H
    return
    (
        rhoRef_*dt*
        (
            (
                coefft*U*U
              - coefft0*U0*U0
              + coefft00*U00*U00
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::DavidsonCurle::DavidsonCurle
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    initialised_(false),
    patches_(0),
    cellZoneName_(word::null),
    cellZoneID_(-1),
    pName_(word::null),
    UName_(word::null),
    rhoRef_(-1),
    cRef_(-1),
    observers_(0)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::functionObjects::DavidsonCurle::~DavidsonCurle()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::DavidsonCurle::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    patches_ =
        mesh_.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));

    cellZoneName_ = dict.lookupOrDefault<word>("cellZone", word::null);

    cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

    pName_ = dict.lookupOrDefault<word>("pName", "p");

    UName_ = dict.lookupOrDefault<word>("UName", "U");

    rhoRef_ = readScalar(dict.lookup("rhoRef"));

    cRef_ = readScalar(dict.lookup("cRef"));

    // Read observers
    {
        const dictionary& obsDict = dict.subDict("observers");
        wordList obsNames = obsDict.toc();

        forAll(obsNames, obsI)
        {
            word oName = obsNames[obsI];
            vector oPos (vector::zero);
            obsDict.subDict(oName).lookup("position") >> oPos;

            observers_.append
            (
                SoundObserver
                (
                    oName,
                    oPos
                )
            );
        }
    }

    return true;
}


bool Foam::functionObjects::DavidsonCurle::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    calculate();

    writeDavidsonCurle();

    return true;
}


bool Foam::functionObjects::DavidsonCurle::write()
{
    return true;
}


void Foam::functionObjects::DavidsonCurle::calculate()
{
    initialise();

    // Pressure field and its time derivative
    volScalarField p( this->p() );
    p.correctBoundaryConditions();
    volScalarField dpdt( this->dpdt() );
    dpdt.correctBoundaryConditions();

    // Lighthill tensor and its time derivatives
    volTensorField Tij( this->Tij() );
    volTensorField dTijdt( this->dTijdt() );
    volTensorField d2Tijdt2( this->d2Tijdt2() );

    // Calculate constant
    scalar coeff = 1.0 / ( 4.0 * Foam::constant::mathematical::pi );

    // Identity tensor
    tensor I(1,0,0,0,1,0,0,0,1);

    // Loop over all observer
    forAll(observers_, obsI)
    {
        SoundObserver& obs = observers_[obsI];
        scalar pPrime = 0.0;

        // Volume integral
        if (cellZoneID_ != -1)
        {
            // List of cells in cellZoneID
            const labelList& cells = mesh_.cellZones()[cellZoneID_];

            // Cell volume and cell center
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();

            // Loop over all cells
            forAll(cells, i)
            {
                label cellI = cells[i];

                // Distance to observer
                scalar r = mag(obs.position() - C[cellI]);
                vector l = (obs.position() - C[cellI]) / r;

                // Calculate pressure fluctuation
                pPrime += coeff * V[cellI] *
                        (
                            ((l*l) && d2Tijdt2[cellI]) / (cRef_ * cRef_ * r)
                          + ((3.0 * l*l - I) && dTijdt[cellI]) / (cRef_ * r * r)
                          + ((3.0 * l*l - I) && Tij[cellI]) / (r * r * r)
                        );
            }
            reduce(pPrime, sumOp<scalar>());
        }



        // Surface integral - loop over all patches
        forAllConstIter(labelHashSet, patches_, iter)
        {
            // Get patch ID
            label patchI = iter.key();

            // Surface area vector and face center at patch
            vectorField Sf = mesh_.Sf().boundaryField()[patchI];
            vectorField Cf = mesh_.Cf().boundaryField()[patchI];

            // Normal vector pointing towards fluid
            vectorField n( -Sf/mag(Sf) );

            // Pressure field and time derivative at patch
            scalarField pp = p.boundaryField()[patchI];
            scalarField dpdtp = dpdt.boundaryField()[patchI];

            // Lighthill tensor on patch
            tensorField Tijp = Tij.boundaryField()[patchI];
            tensorField dTijdtp = dTijdt.boundaryField()[patchI];

            // Distance surface-observer
            scalarField r( mag(obs.position() - Cf) );
            vectorField l( (obs.position() - Cf) / r );

            // Calculate pressure fluctuations
            pPrime += coeff * gSum
            (
                (
                    (l*n)
                  &&
                    (
                        (dpdtp*I - dTijdtp) / (cRef_*r)
                      + (pp*I - Tijp) / sqr(r)
                    )
                )
              * mag(Sf)
            );
        }
        obs.pPrime(pPrime);
    }
}


// ************************************************************************* //
