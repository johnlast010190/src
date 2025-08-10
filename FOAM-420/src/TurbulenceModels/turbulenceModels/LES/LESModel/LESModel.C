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
    (c) 2010-2016 Esi Ltd.
    (c) 2013-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "LESModel.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::LESModel<BasicTurbulenceModel>::LESModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    LESDict_(this->subOrEmptyDict("LES")),
    turbulence_(LESDict_.lookup("turbulence")),
    printCoeffs_(LESDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(LESDict_.optionalSubDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kMin",
            LESDict_,
            sqr(dimVelocity),
            SMALL
        )
    ),

    epsilonMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilonMin",
            LESDict_,
            kMin_.dimensions()/dimTime,
            SMALL
        )
    ),

    omegaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaMin",
            LESDict_,
            dimless/dimTime,
            SMALL
        )
    ),

    nutSmall_("nutSmall", dimLength*dimLength/dimTime, SMALL),
    ySmall_("ySmall", dimLength, SMALL),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", U.group()),
            *this,
            LESDict_
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::LESModel<BasicTurbulenceModel>>
Foam::LESModel<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                IOobject::groupName(propertiesName, U.group()),
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("LES").lookup("LESModel")
    );

    Info<< "Selecting LES turbulence model " << modelType << endl;

    const auto ctor = ctorTableLookup("LESModel type", dictionaryConstructorTable_(), modelType);
    return autoPtr<LESModel>
    (
        ctor(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::LESModel<BasicTurbulenceModel>::yPlus
(
    const label patchNo,
    const scalar Cmu
) const
{
    //Cmu not used, needed for RASModel compatibility

    const fvPatch& curPatch = this->mesh_.boundary()[patchNo];

    tmp<scalarField> tYp(new scalarField(curPatch.size()));
    scalarField& Yp = tYp.ref();

    if (isA<wallFvPatch>(curPatch))
    {
        Yp = this->rho_.boundaryField()[patchNo] * sqrt
            (
                this->muEff(patchNo)
               *mag(this->U_.boundaryField()[patchNo].snGrad())
               /this->rho_.boundaryField()[patchNo]
            )
            /this->mu(patchNo)
            * this->y()[patchNo];
    }
    else
    {
        WarningInFunction
            << "Patch " << patchNo << " is not a wall. Returning null field"
            << nl << endl;

        Yp.setSize(0);
    }

    return tYp;
}

template<class BasicTurbulenceModel>
bool Foam::LESModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        LESDict_ <<= this->subDict("LES");
        LESDict_.lookup("turbulence") >> turbulence_;

        coeffDict_ <<= LESDict_.optionalSubDict(type() + "Coeffs");

        delta_().read(LESDict_);

        kMin_.readIfPresent(LESDict_);

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void Foam::LESModel<BasicTurbulenceModel>::correct()
{
    delta_().correct();
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
