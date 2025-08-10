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
    (c) 2011-2017 OpenFOAM Foundation

Application
    equilibriumCO

Group
    grpThermophysicalUtilities

Description
    Calculates the equilibrium level of carbon monoxide.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "include/OSspecific.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "thermo/thermo/thermo.H"
#include "thermo/janaf/janafThermo.H"
#include "thermo/absoluteEnthalpy/absoluteEnthalpy.H"

#include "containers/LinkedLists/user/SLPtrList.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

#include "cfdTools/general/include/fvCFD.H"

using namespace Foam;

typedef species::thermo<janafThermo<perfectGas<specie>>, absoluteEnthalpy>
    thermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates the equilibrium level of carbon monoxide."
    );
    argList::noParallel();
    argList::noFunctionObjects();

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    word instance = runTime.findInstance(
        polyMesh::meshSubDir, "points", IOobject::READ_IF_PRESENT
    );

    bool meshExists(isFile(runTime.path()/instance/polyMesh::meshSubDir/"faces"));

    fvMesh* meshPtr;

    if (meshExists)
    {
        Info<< "Creating mesh ... " << endl;
        #include "include/createMesh.H"
        meshPtr = &mesh;
    }
    else
    {
        Info<< "mesh not found" << endl;
        Info<< "Creating createSingleCellMesh on the fly" << endl;
        #include "chemFoam/createSingleCellMesh.H"
        meshPtr = &mesh;
    }

    const fvMesh& mesh = *meshPtr;

    Info<< nl << "Reading thermodynamic data IOdictionary" << endl;

    IOdictionary thermoData
    (
        IOobject
        (
            "thermoData",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );



    const scalar P = 1e5;
    const scalar T = 3000.0;

    // Oxidant (mole-based)
    thermo O2(mesh, thermoData.subDict("O2")); O2 *= O2.W();
    thermo N2(mesh, thermoData.subDict("N2")); N2 *= N2.W();

    // Intermediates (mole-based)
    thermo H2(mesh, thermoData.subDict("H2")); H2 *= H2.W();
    thermo OH(mesh, thermoData.subDict("OH")); OH *= OH.W();
    thermo H(mesh, thermoData.subDict("H")); H *= H.W();
    thermo O(mesh, thermoData.subDict("O")); O *= O.W();

    // Products (mole-based)
    thermo CO2(mesh, thermoData.subDict("CO2")); CO2 *= CO2.W();
    thermo H2O(mesh, thermoData.subDict("H2O")); H2O *= H2O.W();
    thermo CO(mesh, thermoData.subDict("CO")); CO *= CO.W();

    SLPtrList<thermo> EQreactions;

    EQreactions.append
    (
        new thermo(CO2 == CO + 0.5*O2)
    );

    EQreactions.append
    (
        new thermo(O2 == 2*O)
    );

    EQreactions.append
    (
        new thermo(H2O == H2 + 0.5*O2)
    );

    EQreactions.append
    (
        new thermo(H2O == H + OH)
    );


    forAllConstIter(SLPtrList<thermo>, EQreactions, iter)
    {
        Info<< "Kc(EQreactions) = " << iter().Kc(P, T) << endl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
