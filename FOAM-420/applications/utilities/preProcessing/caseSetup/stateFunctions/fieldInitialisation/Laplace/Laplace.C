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

#include "Laplace.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(Laplace, 0);
    addToRunTimeSelectionTable(fieldInit, Laplace, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::Laplace::Laplace
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::Laplace::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    label nNonOrthCorr
        = initDict().lookupOrDefault<label>("nNonOrthogonalCorrectors", 10);


    dimensionedScalar alphaL
    (
        "alpha",
        dimensionSet(0,2,-1,0,0,0,0),
        scalar(initDict().lookupOrDefault<scalar>("alpha", 1.0))
    );

    //generate unit 1/AU field for boundary condition compatibility
    volScalarField rUA
    (
        IOobject
        (
            initDict().lookupOrDefault<word>("coefficientFieldName", "(1|A(U))"),
            mesh().time().timeName(),
            localDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("(1|A(U))", dimTime, alphaL.value()),
        zeroGradientFvPatchScalarField::typeName
    );

    fvMesh& refMesh = const_cast<fvMesh&>(mesh());

    dictionary fvSchemes = refMesh.schemes().localSchemeDict();
    {
        if (!fvSchemes.found("gradSchemes"))
        {
            fvSchemes.add
            (
                word("gradSchemes"), dictionary(), false
            );
        }
        char fGrad[] = "Gauss linear";

        fvSchemes.subDict("gradSchemes").add
        (
            word("grad("+name()+")"), fGrad, false
        );

        fvSchemes.subDict("gradSchemes").add
        (
            word("snGradCorr("+name()+")"), fGrad, false
        );


        if (!fvSchemes.found("laplacianSchemes"))
        {
            fvSchemes.add(word("laplacianSchemes"), dictionary(), false);
        }

        char fLapl[]
            = "Gauss linear limited 0.333";

        fvSchemes.subDict("laplacianSchemes").add
        (
            word("laplacian(alpha,"+name()+")"), fLapl, false
        );

        if (!fvSchemes.found("fluxRequired"))
        {
            fvSchemes.add(word("fluxRequired"), dictionary(), false);
        }
        fvSchemes.subDict("fluxRequired").add(name(), "");


        if (!fvSchemes.found("interpolationSchemes"))
        {
            fvSchemes.add(word("interpolationSchemes"), dictionary(), false);
        }
        fvSchemes.subDict("interpolationSchemes").add
        (
            word("interpolate("+name()+")"),
            "linear",
            false
        );

        refMesh.schemes().setLocalSchemeDict(fvSchemes);
    }

    //solve steady Laplace's equation based on field type

    Info<< "Solving the Steady Laplace Equation" << endl;

    if (localDb().foundObject<volVectorField>(name()))
    {
        volVectorField& f
            = const_cast<volVectorField&>
            (localDb().lookupObject<volVectorField>(name()));

        for (int i=0; i<=nNonOrthCorr; i++)
        {
            solve
            (
                fvm::laplacian(alphaL, f)
            );
        }

    }
    else if (localDb().foundObject<volScalarField>(name()))
    {
        volScalarField& f
            = const_cast<volScalarField&>
            (localDb().lookupObject<volScalarField>(name()));

        for (int i=0; i<=nNonOrthCorr; i++)
        {
            solve
            (
                fvm::laplacian(alphaL, f)
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Trying to initialise unsupported field type."
            << exit(FatalError);
    }

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
