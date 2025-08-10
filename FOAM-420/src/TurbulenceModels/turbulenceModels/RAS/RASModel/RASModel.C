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
    (c) 2010-2012 Esi Ltd.
    (c) 2013-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "RASModel.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::RASModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::RASModel<BasicTurbulenceModel>::RASModel
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

    RASDict_(this->subOrEmptyDict("RAS")),
    turbulence_(RASDict_.lookup("turbulence")),
    printCoeffs_(RASDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(RASDict_.optionalSubDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kMin",
            RASDict_,
            sqr(dimVelocity),
            SMALL
        )
    ),

    epsilonMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilonMin",
            RASDict_,
            kMin_.dimensions()/dimTime,
            SMALL
        )
    ),

    omegaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaMin",
            RASDict_,
            dimless/dimTime,
            SMALL
        )
    ),

    nutSmall_("nutSmall", dimLength*dimLength/dimTime, SMALL),
    ySmall_("ySmall", dimLength, SMALL)
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::RASModel<BasicTurbulenceModel>>
Foam::RASModel<BasicTurbulenceModel>::New
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
        ).subDict("RAS").lookup("RASModel")
    );

    Info<< "Selecting RAS turbulence model " << modelType << endl;

    const auto ctor = ctorTableLookup("RASModel type", dictionaryConstructorTable_(), modelType);
    return autoPtr<RASModel>
    (
        ctor(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::scalar Foam::RASModel<BasicTurbulenceModel>::yPlusLam
(
    const scalar kappa,
    const scalar E
) const
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField> Foam::RASModel<BasicTurbulenceModel>::yPlus
(
    const label patchNo,
    const scalar Cmu
) const
{
    const fvPatch& curPatch = this->mesh_.boundary()[patchNo];

    tmp<scalarField> tYp(new scalarField(curPatch.size()));
    scalarField& Yp = tYp.ref();

    if (isA<wallFvPatch>(curPatch))
    {
         Yp =
         sqrt
         (
             this->nuEff()().boundaryField()[patchNo]
            *mag(this->U_.boundaryField()[patchNo].snGrad())
         )
         /this->nu()().boundaryField()[patchNo]
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
Foam::tmp<Foam::volScalarField> Foam::RASModel<BasicTurbulenceModel>::L() const
{
    FatalError << "Length scale function, L(), not implemented for: "
               << IOdictionary::lookup("RASModel")
               << exit(FatalError);

    return tmp<volScalarField>();
}



template<class BasicTurbulenceModel>
bool Foam::RASModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        RASDict_ <<= this->subDict("RAS");
        RASDict_.lookup("turbulence") >> turbulence_;

        coeffDict_ <<= RASDict_.optionalSubDict(type() + "Coeffs");

        kMin_.readIfPresent(RASDict_);
        epsilonMin_.readIfPresent(RASDict_);
        omegaMin_.readIfPresent(RASDict_);
        nutSmall_.readIfPresent(*this);
        ySmall_.readIfPresent(*this);

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void Foam::RASModel<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
