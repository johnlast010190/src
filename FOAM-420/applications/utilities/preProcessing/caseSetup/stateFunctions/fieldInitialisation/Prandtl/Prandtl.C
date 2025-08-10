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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "Prandtl.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "stateFunction/stateFunction.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(Prandtl, 0);
    addToRunTimeSelectionTable(fieldInit, Prandtl, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::Prandtl::Prandtl
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
    // requires velocity field to be initialised first
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::Prandtl::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    //always scalar field
    volScalarField& f
    (
        const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()))
    );

    //wall dist
    volScalarField d = wallDist(mesh()).y();

    //velocity field
    word UName = word
    (
        initDict().lookupOrDefault<word>("U","U")
    );

    if (localDb().foundObject<volVectorField>(UName))
    {
        const volVectorField& U
        (
            localDb().lookupObject<volVectorField>(UName)
        );

        //calculate reasonable boundary layer thickness
        //mean wall dist * constant

        scalar cBL = 0.5;

        dimensionedScalar blThickness
        (
            "boundaryLayerThickness",
            dimLength,
            gAverage(d)*cBL
        );

        blThickness.value() = initDict().lookupOrDefault<scalar>
        (
            "boundaryLayerThickness",
            blThickness.value()
        );


        Info<< "   Assumed boundary thickness = " << blThickness.value()
             << "   m." << endl;

        d = min(d, blThickness);

        //read turbulenceProperties dict
        IOdictionary  turbProps
        (
            IOobject
            (
                "turbulenceProperties",
                mesh().time().constant(),
                localDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );


        if
        (

            stateFunction
            ::turbulenceTypeNames_[word(turbProps.lookup("simulationType"))]
         == stateFunction::tuLES
        )
        {
            forAll(d, cI)
            {
                d[cI] = min(d[cI], ::pow(mesh().V()[cI], (1.0/3.0)));
            }
        }

        dimensionedScalar kappa("kappa", dimless, 0.4187);
        dimensionedScalar Cmu("Cmu", dimless, 0.09);

        //needs U, use UPtr for now
        //insert grad scheme for U
        fvMesh& refMesh = const_cast<fvMesh&>(mesh());
        dictionary fvSchemes = refMesh.schemes().localSchemeDict();
        {
            if (!fvSchemes.found("gradSchemes"))
            {
                fvSchemes.add(word("gradSchemes"), dictionary(), false);
            }
            if (!fvSchemes.subDict("gradSchemes").found("grad(U)"))
            {
                char UGrad[] = "Gauss linear";
                fvSchemes.subDict("gradSchemes").add
                (
                    word("grad(U)"), UGrad, false
                );
            }

            if (!fvSchemes.found("interpolationSchemes"))
            {
                fvSchemes.add
                (
                    word("interpolationSchemes"),
                    dictionary(),
                    false
                );
            }
            if
            (
                !fvSchemes.subDict("interpolationSchemes")
                .found("interpolate(U)")
            )
            {
                fvSchemes.subDict("interpolationSchemes").add
                (
                    word("interpolate(U)"), "linear", false
                );
            }

            refMesh.schemes().setLocalSchemeDict(fvSchemes);
        }


        volScalarField nut
        (
            sqr(kappa*d)
            *::sqrt(2.0)*mag(symm(fvc::grad(U)))

        );

        //limit nut field
        {
            scalar Imax = 0.2;
            scalar Imin = 0.001;
            scalar ck0 = ::pow(Cmu.value(), 0.25) * kappa.value();

            scalarField sqrtkmax( sqrt(0.5*magSqr(Imax*U.primitiveField())) );
            scalarField sqrtkmin( sqrtkmax*Imin/Imax );

            scalar nuMin = initDict().lookupOrDefault<scalar>("nu", 1e-5);

            nut.primitiveFieldRef() = min
            (
                max(ck0*d*sqrtkmin, nut.primitiveField()),
                ck0*d*sqrtkmax
            ) + 0.5 * nuMin;
        }

        if (name() == "k" || name() == "kl")
        {
            scalar ck0 = ::pow(Cmu.value(), 0.25) * kappa.value();

            f.primitiveFieldRef()
                = ::sqr(nut.primitiveField()/(ck0*d));
        }
        else if (name() == "epsilon")
        {
            scalar ck0 = ::pow(Cmu.value(), 0.25) * kappa.value();
            scalar ce0 = ::pow(Cmu.value(), 0.75)/kappa.value();

            f.primitiveFieldRef()
                = ce0*Foam::pow(::sqr(nut.primitiveField()/(ck0*d)), 1.5)/d;
        }
        else if (name() == "nuTilda")
        {
            f.primitiveFieldRef() = nut.primitiveField();
        }
        else if (name() == "nut")
        {
            f.primitiveFieldRef() = nut.primitiveField();
        }
        else if (name() == "mut")
        {
            scalar rho = readScalar(initDict().lookup("rho"));

            f.primitiveFieldRef() = rho*nut.primitiveField();
        }
        else if (name() == "omega")
        {
            scalar ck0 = ::pow(Cmu.value(), 0.25) * kappa.value();

            f.primitiveFieldRef()
                = nut.primitiveField()/::sqr(ck0*d.primitiveField());
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find velocity field "<< UName
            <<" to initialize turbulent values"
            << exit(FatalError);
    }


    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
    initTurbBoundaries();
}



// ************************************************************************* //
