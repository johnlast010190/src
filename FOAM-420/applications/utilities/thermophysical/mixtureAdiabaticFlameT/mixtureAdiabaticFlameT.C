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
    (c) 2011-2016 OpenFOAM Foundation

Application
    mixtureAdiabaticFlameT

Group
    grpThermophysicalUtilities

Description
    Calculates the adiabatic flame temperature for a given mixture
    at a given temperature.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "include/OSspecific.H"
#include "global/etcFiles/etcFiles.H"

#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "thermo/thermo/thermo.H"
#include "thermo/janaf/janafThermo.H"
#include "thermo/absoluteEnthalpy/absoluteEnthalpy.H"
#include "mixture.H"

#include "cfdTools/general/include/fvCFD.H"

using namespace Foam;

typedef species::thermo<janafThermo<perfectGas<specie>>, absoluteEnthalpy>
    thermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates the adiabatic flame temperature for a given mixture\n"
        "at a given temperature."
    );
    argList::noParallel();
    argList::noFunctionObjects();
    argList::validArgs.append("controlFile");
    argList args(argc, argv);

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

    const fileName controlFileName(args[1]);

    // Construct control dictionary
    IFstream controlFile(controlFileName);

    // Check controlFile stream is OK
    if (!controlFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << controlFileName
            << abort(FatalError);
    }

    dictionary control(controlFile);


    scalar P(readScalar(control.lookup("P")));
    scalar T0(readScalar(control.lookup("T0")));
    mixture rMix(control.lookup("reactants"));
    mixture pMix(control.lookup("products"));


    Info<< nl << "Reading thermodynamic data dictionary" << endl;

    fileName thermoDataFileName(findEtcFile("thermoData/thermoData"));

    // Construct control dictionary
    IFstream thermoDataFile(thermoDataFileName);

    // Check thermoData stream is OK
    if (!thermoDataFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << thermoDataFileName
            << abort(FatalError);
    }

    dictionary thermoData(thermoDataFile);


    thermo reactants
    (
        rMix[0].volFrac()*thermo(mesh, thermoData.subDict(rMix[0].name()))
    );

    for (label i = 1; i < rMix.size(); i++)
    {
        reactants = reactants
            + rMix[i].volFrac()*thermo(mesh, thermoData.subDict(rMix[i].name()));
    }


    thermo products
    (
        2*pMix[0].volFrac()*thermo(mesh, thermoData.subDict(pMix[0].name()))
    );

    for (label i = 1; i < pMix.size(); i++)
    {
        products = products
            + 2*pMix[i].volFrac()*thermo(mesh, thermoData.subDict(pMix[i].name()));
    }

    Info<< "Adiabatic flame temperature of mixture " << rMix.name() << " = "
         << products.THa(reactants.Ha(P, T0), P, 1000.0) << " K" << endl;

    return 0;
}


// ************************************************************************* //
