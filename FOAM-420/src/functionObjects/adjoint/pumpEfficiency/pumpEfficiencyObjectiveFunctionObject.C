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

#include "pumpEfficiency/pumpEfficiencyObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/derived/resistivePressure/resistivePressureFvPatchScalarField.H"
#include "sources/derived/MRFSource/MRFSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pumpEfficiencyObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        pumpEfficiencyObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void
Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::costDp() const
{
    scalar powerGain = 0;
    scalar mass = 0;

    tmp<volScalarField> rhoInf = rho();

    forAll(mesh_.boundary(), patchI)
    {
        if (objectivePatch_[patchI])
        {
            const fvPatch& p = mesh_.boundary()[patchI];

            const scalarField& rhop(rhoInf().boundaryField()[patchI]);
            const vectorField& Vp(U().boundaryField()[patchI]);

            scalarField pp
            (
                objectiveFunctionObject::p().boundaryField()[patchI]
            );
            scalarField phip(phi().boundaryField()[patchI]);

            //grab outside efficiency from resistive outlet
            if
            (
                isA<resistivePressureFvPatchScalarField>
                (
                    objectiveFunctionObject::p().boundaryField()[patchI]
                )
            )
            {
                pp =
                    dynamic_cast<const resistivePressureFvPatchScalarField&>
                    (
                        objectiveFunctionObject::p().boundaryField()[patchI]
                    ).p0();
            }

            // Different formulation based on flow solver
            if (phi().dimensions() == dimensionSet(0,3,-1,0,0)) //incompressible
            {
                // Calculate inlet and outlet power
                if (p.type() == "inlet" || p.patch().physicalType() == "inlet")
                {
                    powerGain += sum
                    (
                        rhop*phip*(pp + 0.5*magSqr(Vp))
                    );
                }
                else if (p.type() == "outlet" || p.patch().physicalType() == "outlet")
                {
                    powerGain += sum
                    (
                        rhop*phip*(pp + 0.5*magSqr(Vp))
                    );
                }

                mass += sum(phip);

            }
            else //compressible
            {
                //different treatment based on type of thermo used
                scalarField h(p.size(), 0.0);
                mass += sum(phip);

                if (adiabatic_)
                {
                    h = pp/rhop;
                }
                else if (mesh_.foundObject<volScalarField>("h"))
                {
                    h = p.lookupPatchField<volScalarField, scalar>("h");
                }
                else if (mesh_.foundObject<volScalarField>("e"))
                {
                    h = p.lookupPatchField<volScalarField, scalar>("e")
                        + pp/rhop;
                }
                else
                {
                    FatalErrorInFunction
                        << "Could not find sensible enthalpy or internal "
                        << "energy field in database." << exit(FatalError);
                }

                if
                (
                    p.type() == "inlet"
                 || p.patch().physicalType() == "inlet"
                )
                {
                    powerGain -= sum
                    (
                        phip * (h + 0.5*magSqr(Vp))
                    );
                }
                else if
                (
                    p.type() == "outlet"
                 || p.patch().physicalType() == "outlet"
                )
                {
                    powerGain -= sum
                    (
                        phip * (h + 0.5*magSqr(Vp))
                    );
                }
            }
        }
    }

    reduce(powerGain, sumOp<scalar>());
    reduce(mass, sumOp<scalar>());

    autoPtr<scalar>& cpGPtr =  const_cast<autoPtr<scalar>&>(pressureGain_);
    if (!cpGPtr.valid())
    {
       cpGPtr.reset(new scalar(powerGain));
    }
    else
    {
        cpGPtr() = powerGain;
    }
}


void
Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::costTw() const
{
    tmp<volScalarField> P = objectiveFunctionObject::P();
    const volScalarField::Boundary& Pw
        = P().boundaryField();
    tmp<volSymmTensorField> tau = devRhoReff();
    const volSymmTensorField::Boundary& tauw
        = tau().boundaryField();
    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    scalar power = 0.0;

    fv::options& fvOpt(fv::options::New(this->mesh_));

    for (int fvI=0; fvI < (fvOpt.optionList::size()); fvI++)
    {
        fv::option& fvOptI = fvOpt.optionList::operator[](fvI);
        if (fvOptI.isActive())
        {
            if (fvOptI.type() == "MRFSource")
            {
                fv::MRFSource& mrfS = refCast<fv::MRFSource>(fvOptI);

                label zoneNamesN = zoneNames_.size();
                bool zoneOpt(false);
                if (zoneNamesN == 0)
                {
                    zoneOpt = true;
                }
                forAll(zoneNames_, zNI)
                {
                    if (zoneNames_[zNI] == mrfS.cellSetName())
                    {
                        zoneOpt = true;
                    }
                }
                if (zoneOpt)
                {
                    const vector& axisPatch = mrfS.axis();
                    const scalar& omegaPatch = mag(mrfS.Omega());
                    const vector& centerOfRotation = mrfS.CofR();
                    const labelListList& incFaces = mrfS.getFrameSourceFaces().includedFaces();
                    vector omegaV = axisPatch*omegaPatch;
                    forAll(mesh_.boundary(), pI)
                    {
                        if (objectivePatch_[pI])
                        {
                            tmp<vectorField> normal = mesh_.boundary()[pI].nf();
                            vectorField pf(Sfb[pI]*Pw[pI]);
                            vectorField vf(Sfb[pI]&tauw[pI]);
                            forAll(incFaces[pI],i)
                            {
                                label pF = incFaces[pI][i];
                                vector Cp_n = mesh_.C().boundaryField()[pI][pF];
                                scalar Cp_naxis = Cp_n & axisPatch;
                                vector localVector =
                                (
                                    (Cp_n - (Cp_naxis * axisPatch)- centerOfRotation)^omegaV
                                );
                                power += ((pf[pF]+vf[pF]) & localVector);
                            }
                        }
                    }
                }
            }
        }
    }

    reduce(power, sumOp<scalar>());
    autoPtr<scalar>& cTOPtr =  const_cast<autoPtr<scalar>&>(TOmega_);
    if (!cTOPtr.valid())
    {
        cTOPtr.reset(new scalar(power));
    }
    else
    {
        cTOPtr() = power;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::
pumpEfficiencyObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    adiabatic_
    (
        objectiveDict.lookupOrDefault<Switch>("adiabatic", false)
    ),
    pressureGain_(nullptr),
    TOmega_(nullptr)
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::
pumpEfficiencyObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    pumpEfficiencyObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::
~pumpEfficiencyObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    adiabatic_ = dict.lookupOrDefault<Switch>("adiabatic", false);

    return true;
}


bool Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::execute()
{
    costDp();
    costTw();
    if (TOmega_()!=0.)
    {
        objectiveValue_ = pressureGain_()/TOmega_() * scalar(100);
    }
    else
    {
        WarningInFunction<< "Pump power equals to zero " << nl
        <<"Please check you MRF definition "<<endl;
        objectiveValue_ = 0;
    }

    Info<< type() << " " << name() << " execute:" << nl
        << "Pump efficiency = " << objectiveValue_ << " [%]" << nl << endl;

    return true;
}


bool Foam::functionObjects::pumpEfficiencyObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
