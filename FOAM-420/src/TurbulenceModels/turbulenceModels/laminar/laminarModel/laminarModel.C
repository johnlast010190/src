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
    (c) 2016-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "laminarModel.H"
#include "laminar/Stokes/Stokes.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::laminarModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::laminarModel<BasicTurbulenceModel>::laminarModel
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

    laminarDict_(this->subOrEmptyDict("laminar")),
    printCoeffs_(laminarDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(laminarDict_.optionalSubDict(type + "Coeffs"))
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::laminarModel<BasicTurbulenceModel>>
Foam::laminarModel<BasicTurbulenceModel>::New
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
    IOdictionary modelDict
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
    );

    if (modelDict.found("laminar"))
    {
        // get model name, but do not register the dictionary
        // otherwise it is registered in the database twice
        const word modelType
        (
            modelDict.subDict("laminar").lookup("laminarModel")
        );

        Info<< "Selecting laminar stress model " << modelType << endl;

        const auto ctor = ctorTableLookup("laminar stress model", dictionaryConstructorTable_(), modelType);
        return autoPtr<laminarModel>
        (
            ctor
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport, propertiesName)
        );
    }
    else
    {
        Info<< "Selecting laminar stress model "
            << laminarModels::Stokes<BasicTurbulenceModel>::typeName << endl;

        return autoPtr<laminarModel>
        (
            new laminarModels::Stokes<BasicTurbulenceModel>
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport,
                propertiesName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::laminarModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        laminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= laminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::nut() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("nut", this->U_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("nut", dimViscosity, 0.0)
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::nuEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("nuEff", this->U_.group()), this->nu()
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->U_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("k", sqr(this->U_.dimensions()), 0.0)
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->U_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar
            (
                "epsilon", sqr(this->U_.dimensions())/dimTime, 0.0
            )
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::laminarModel<BasicTurbulenceModel>::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", this->U_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedSymmTensor
            (
                "R", sqr(this->U_.dimensions()), Zero
            )
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField> Foam::laminarModel<BasicTurbulenceModel>::yPlus
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
        Yp = sqrt
            (
                nuEff()().boundaryField()[patchNo]
               *mag(this->U_.boundaryField()[patchNo].snGrad())
            )
            /this->nu()().boundaryField()[patchNo]/curPatch.deltaCoeffs();
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
void Foam::laminarModel<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
