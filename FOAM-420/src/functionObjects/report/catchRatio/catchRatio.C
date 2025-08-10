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
    (c) 2011 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "catchRatio.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(catchRatio, 0);
    addToRunTimeSelectionTable(functionObject, catchRatio, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volVectorField&
Foam::functionObjects::catchRatio::Ud
(
    const label idx
) const
{
    return lookupObject<volVectorField>("U." + phaseNames_[idx]);
}


const Foam::volScalarField&
Foam::functionObjects::catchRatio::alphad
(
    const label idx
) const
{
    return lookupObject<volScalarField>("alpha." + phaseNames_[idx]);
}

void Foam::functionObjects::catchRatio::createAndStoreFields()
{
    word fieldname("gcr");
    tmp<volScalarField> gcr
    (
        new volScalarField
        (
            IOobject
            (
                fieldname,
                mesh_.time().timeName(),
                obr(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("gcr", dimless, 0.0)
        )
    );
    store(fieldname, gcr);

    forAll(phaseNames_, phaseI)
    {
        fieldname = "scr" + phaseNames_[phaseI];
        tmp<volScalarField> scr
        (
            new volScalarField
            (
                IOobject
                (
                    fieldname,
                    mesh_.time().timeName(),
                    obr(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("scr", dimless, 0.0)
            )
        );
        store(fieldname, scr);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::catchRatio::catchRatio
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseNames_(),
    psd_(),
    Rh_()
{
    read(dict);
    createAndStoreFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::catchRatio::~catchRatio()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::catchRatio::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    dict.lookup("psd") >> psd_;
    Rh_ = dict.lookupOrDefault<scalar>("Rh", 1.0);

    IOdictionary disperseEulerianProperties
    (
        IOobject
        (
            "disperseEulerianProperties",
            mesh_.time().constant(),
            obr(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    PtrList<entry> phaseEntries(disperseEulerianProperties.lookup("phases"));

    if (psd_.size() != phaseEntries.size())
    {
        FatalErrorInFunction
            << "Size of specied scalarList psd: "
            << psd_.size() << nl
            << " does not match the number of disperse phases: "
            << phaseEntries.size() << nl
            << "Please check your functionObject setup." << endl
            << exit(FatalError);
    }

    phaseNames_.setSize(phaseEntries.size());

    forAll(phaseNames_, phaseI)
    {
        phaseNames_[phaseI] = phaseEntries[phaseI].keyword();
    }

    return true;
}


bool Foam::functionObjects::catchRatio::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;
    Log << tab << "Computing local and global catch ratio." << nl << endl;

// compute local (phasic) and global catch ratio
    volScalarField& gcr =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>("gcr")
        );
    gcr.forceAssign(0.0); // re-set to zero

    // for all spray fields
    forAll(phaseNames_, phaseI)
    {
        surfaceScalarField normalvel
        (
            mag((mesh_.Sf()/mesh_.magSf()) & fvc::interpolate(Ud(phaseI)))
        );
        surfaceScalarField surfaceScr
        (
            (normalvel * fvc::interpolate(alphad(phaseI))) * ((3600*1E3)/(Rh_*psd_[phaseI]))
        );
        const surfaceScalarField::Boundary& patchSurfaceScr = surfaceScr.boundaryField();

        word fieldname("scr" + phaseNames_[phaseI]);
        volScalarField& scr =
            const_cast<volScalarField&>
            (
                lookupObject<volScalarField>(fieldname)
            );
        volScalarField::Boundary& scrBf = scr.boundaryFieldRef();

        forAll(scrBf,patchi)
        {
            scrBf[patchi] = patchSurfaceScr[patchi];
        }

        scr.write();
        gcr += scr*psd_[phaseI];
        gcr.write();
    }

    return true;
}


bool Foam::functionObjects::catchRatio::write()
{
    return true;
}


// ************************************************************************* //
