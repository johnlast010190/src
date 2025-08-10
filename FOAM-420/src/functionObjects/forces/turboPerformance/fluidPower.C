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
    (c) D. Boger, B. Lewis M. Auvinen, H. Nilsson

\*---------------------------------------------------------------------------*/

#include "turboPerformance/fluidPower.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"    // added
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fluidPower, 0);
    addToRunTimeSelectionTable(functionObject, fluidPower, dictionary);
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Same constructor layout as in 'forces'.
// The differences lie only in the private data --mikko
Foam::functionObjects::fluidPower::fluidPower
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool baseClass
)
:
    forces(name, runTime, dict, true),
    turbine_(false),          // Set to true if turbine  (Bryan)
    inletPatchSet_(),         // new labelHashSet
    outletPatchSet_(),        //  -- " --    -- mikko
    phiName_("phi"),
    inletSectors_(dict.lookupOrDefault<label>("inletSectors", 1)),
    outletSectors_(dict.lookupOrDefault<label>("outletSectors", 1)),
    volflux_(0)
{
    if (!baseClass)
    {
        read(dict);
        resetFile(fluidPower::typeName);
        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fluidPower::~fluidPower()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::functionObjects::fluidPower::writeFileHeader(Ostream& os) const
{
    writeCommented(os,"Time");
    writeDelimited(os,"dEm (W)");
    writeDelimited(os,"Head (m)");
    os  << endl;
}


// The computation of dEm and Head
Foam::functionObjects::fluidPower::dEmHead
Foam::functionObjects::fluidPower::calcDEmHead() const
{
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);
    tmp<volScalarField> trho = rho();

    dEmHead dEmH( scalar(0) , scalar(0) );

    const volVectorField::Boundary& Ub = U.boundaryField();
    const volScalarField::Boundary& pb = p.boundaryField();
    const surfaceScalarField::Boundary&  phib = phi.boundaryField();

    scalar mflowInlet  = 0.0;        // mass flow into the domain
    scalar mflowOutlet = 0.0;        // mass flow out of the domain
                                     //(currently not needed)
    scalar EmInlet     = 0.0;        // Mechanical energy flow into
                                     //the domain (in Watts!)
    scalar EmOutlet    = 0.0;        // Mechanical energy flow out of the domain

    // inflow of mechanical energy
    forAllConstIter(labelHashSet, inletPatchSet_, iter)
    {
        label patchi = iter.key();

        scalarField ptInlet( 0.5*magSqr(Ub[patchi]));
        if (p.dimensions() == dimPressure)
        {
            ptInlet += pb[patchi]/trho().boundaryField()[patchi];
        }
        else
        {
            ptInlet += pb[patchi];
        }

        mflowInlet  += rho(p)* sum( phib[patchi]) + VSMALL;  // Sign: Inflow (-)
        EmInlet     += rho(p)* sum( phib[patchi] * ptInlet );
    }

    // outflow of mechanical energy
    forAllConstIter(labelHashSet, outletPatchSet_, iter)
    {
        label patchi = iter.key();

        scalarField ptOutlet( 0.5*magSqr(Ub[patchi]));
        if (p.dimensions() == dimPressure)
        {
            ptOutlet += pb[patchi]/trho().boundaryField()[patchi];
        }
        else
        {
            ptOutlet += pb[patchi];
        }

        mflowOutlet += rho(p)*sum( phib[patchi])+VSMALL; // - Signs: Outflow(+)
        EmOutlet    += rho(p)* sum( phib[patchi] * ptOutlet );
    }


    reduce(EmInlet , sumOp<scalar>());
    reduce(EmOutlet , sumOp<scalar>());
    reduce(mflowInlet , sumOp<scalar>());
    reduce(mflowOutlet , sumOp<scalar>());

    volflux_ = (Foam::mag(scalar(inletSectors_*mflowInlet))
        +Foam::mag(scalar(outletSectors_*mflowOutlet)))/(2*rho(p));

    //dEmH.first() = (EmOutlet + EmInlet ); // The rate of work output
                                            // from the system (W)
    //dEmH.second() = (EmOutlet + EmInlet )/(-mflowInlet*scalar(9.81));
                                           // Hydrodynamic (total) head

    if (turbine_) // (Bryan)
    {
        scalar dE =
            scalar(-1.0)*(outletSectors_*EmOutlet + inletSectors_*EmInlet );
        dEmH.first() = dE;
        // The rate of work output from the system (W)
        dEmH.second() = dE/(-inletSectors_*mflowInlet*scalar(9.81));
        // Hydrodynamic (total) head
    }
    else
    {
        scalar dE = outletSectors_*EmOutlet + inletSectors_*EmInlet ;
        dEmH.first() = dE;
        // The rate of work output from the system (W)
        dEmH.second() = dE/(-inletSectors_*mflowInlet*scalar(9.81));
        // Hydrodynamic (total) head
    }

    return dEmH;
}


bool Foam::functionObjects::fluidPower::read(const dictionary& dict)
{
    // duplicate neccessary force read functions
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;
    nBin_ = 1;

    Info<< type() << " " << name() <<  " read:" << nl;

    // Field names
    wordList pNames(2, "p"); pNames[1] = "pName";
    pName_ = dict.lookupOrDefault<word>(pNames, "p");

    wordList UNames(2, "U"); UNames[1] = "UName";
    UName_ = dict.lookupOrDefault<word>(UNames, "U");

    wordList rhoNames(2, "rho"); rhoNames[1] = "rhoName";
    rhoName_ = dict.lookupOrDefault<word>(rhoNames, "rho");

    wordList phiNames(2, "phi"); rhoNames[1] = "phiName";
    phiName_ = dict.lookupOrDefault<word>(phiNames, "phi");

    // Reference density needed for incompressible calculations
    if (rhoName_ == "rhoInf")
    {
        dict.lookup("rhoInf") >> rhoRef_;
    }

    // Reference pressure, 0 by default
    pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
    if (dict.isDict("referenceFields"))
    {
        const dictionary& refDict = dict.subDict("referenceFields");
        pRef_ =
            refDict.found("p")
          ? refDict.lookup<dimensionedScalar>("p").value()
          : 0;
    }
    coordSys_.origin() = dict.lookupOrDefault<vector>("CofR", Zero);

    // Fluid power specific read
    // Read if pump or turbine (Bryan)
    turbine_ = dict.lookupOrDefault<Switch>("turbine", false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    inletPatchSet_ = pbm.patchSet(wordReList(dict.lookup("inletPatches")));
    outletPatchSet_ = pbm.patchSet(wordReList(dict.lookup("outletPatches")));

    return true;
}


bool Foam::functionObjects::fluidPower::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    // Obtain the difference dEm = (Em_output - Em_input) and
    //the hydrodynamic head
    dEmHead  dEmH = calcDEmHead();

    file() << obr_.time().timeName()
           << tab << dEmH.first() << tab
           << dEmH.second()  << endl;

    Log<< " Fluid power output:" << nl
       << "  dEm (W) = " << dEmH.first() << nl
       << "  Head (m) = " << dEmH.second() << nl
       << endl;

    return true;
}


bool Foam::functionObjects::fluidPower::write()
{
    return true;
}

// ************************************************************************* //
