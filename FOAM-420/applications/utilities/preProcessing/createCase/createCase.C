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
    (c) 2008 Icon CG Ltd.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "primitives/strings/fileName/fileName.H"
#include "primitives/strings/lists/stringList.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("case name");
    argList::validOptions.insert("nProcs", "number of processors");
    argList::validOptions.insert("force", "");

    Foam::argList::addOption
    (
        "templateDir",
        "file",
        "read case set-up templates from specified location"
    );

    argList args(argc, argv);
    fileName caseName(args[1]);

    string newCase = args.path()/caseName;

    Info<< endl;

    if (args.optionFound("force"))
    {
        Info<< "Overwriting " << newCase << endl;
        #if defined( WIN32 ) || defined( WIN64 )
        Foam::system(string("rmdir \"" + newCase +"\" \/"+"s"+" \/"+"q").c_str());
        #else
        Foam::system(string("rm -rf " + newCase).c_str());
        #endif

    }
    else
    {
        Info<< "Creating new case : " << newCase << endl;
    }


    if (Foam::isDir(newCase))
    {
        FatalError << "Case directory " << newCase << " already exists."
             << exit(FatalError);
    }
    else
    {
        Foam::mkDir(newCase);
    }

    Foam::mkDir(newCase/"system");
    Foam::mkDir(newCase/"log");
    Foam::mkDir(newCase/"constant");
    Foam::mkDir(newCase/"constant"/"triSurface");

    // create case.foam file for Paraview
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "type nul >> \""
            + newCase/caseName.name()+".foam\" "
        ).c_str()
        #else
        string
        (
            "touch "
            + newCase/caseName.name()+".foam"
        ).c_str()
        #endif

    );

    #if defined( WIN32 ) || defined( WIN64 )
    fileName baseDir
    (
//        "%FOAM_PROJECT_DIR%\\etc\\caseDicts\\preProcessing\\createCase"
        "${FOAM_PROJECT_DIR}/etc/caseDicts/preProcessing/createCase"
    );
    #else
    fileName baseDir
    (
        "${FOAM_PROJECT_DIR}/etc/caseDicts/preProcessing/createCase"
    );
    #endif

    if (args.optionFound("templateDir"))
    {
        baseDir = args["templateDir"];
    }

    baseDir.expand();
    baseDir.toAbsolute();

    if (!isDir(baseDir))
    {
        FatalErrorInFunction
            << "templateDir " << baseDir
            << " should point to the folder containing the "
            << "case set-up templates" << exit(FatalError);
    }


    //controlDict
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.controlDict\" \""
            + newCase + "\\system\\controlDict\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.controlDict "
            + newCase/"system"/"controlDict"
        ).c_str()
        #endif
    );
    //foamHexMeshDict
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.foamHexMeshDict\" \""
            + newCase + "\\system\\foamHexMeshDict\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.foamHexMeshDict "
            + newCase/"system"/"foamHexMeshDict"
        ).c_str()
        #endif
    );
    //caseSetupDict
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.caseSetupDict\" \""
            + newCase + "\\system\\caseSetupDict\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.caseSetupDict "
            + newCase/"system"/"caseSetupDict"
        ).c_str()
        #endif
    );
    //decomposeParDict
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.decomposeParDict\" \""
            + newCase + "\\system\\decomposeParDict\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.decomposeParDict "
            + newCase/"system"/"decomposeParDict"
        ).c_str()
        #endif
    );
    //fvSchemes
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.fvSchemes\" \""
            + newCase + "\\system\\fvSchemes\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.fvSchemes "
            + newCase/"system"/"fvSchemes"
        ).c_str()
        #endif
    );
    //fvSolution
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.fvSolution\" \""
            + newCase + "\\system\\fvSolution\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.fvSolution "
            + newCase/"system"/"fvSolution"
        ).c_str()
        #endif
    );
    //blockMeshDict
    Foam::system
    (
        #if defined( WIN32 ) || defined( WIN64 )
        string
        (
            "copy  \""+ baseDir + "\\createCase.blockMeshDict\" \""
            + newCase + "\\system\\blockMeshDict\" /B"
        ).c_str()
        #else
        string
        (
            "cp " + baseDir + "/createCase.blockMeshDict "
            + newCase/"system"/"blockMeshDict"
        ).c_str()
        #endif
    );

    //check if number of processors is larger than 1
    //if so do factoring to get optmal hierarchical decomposition
    //load decomposeParDict, modify and write back to file.

    if (args.optionFound("nProcs"))
    {
        Time runTime
        (
            Foam::Time::controlDictName,
            args.path(),
            caseName
        );

        IOdictionary decomposePar
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );

        label np = -1;
        args.optionReadIfPresent("nProcs", np);

        decomposePar.add("numberOfSubdomains", np, true);

        //factor nprocs for optmial hierarchical decomposition
        FixedList<label, 3> hdecomp(-1);
        scalar sfactor= Foam::cbrt(scalar(np));
        label factor = sfactor;
        if (factor == 1 && sfactor > 1)
        {
            factor++;
        }

        for (; factor < np+1; factor++)
        {
            if (np%factor == 0)
            {
                break;
            }
        }
        hdecomp[0] = factor;
        label np2 = np/factor;
        sfactor= Foam::sqrt(scalar(np2));
        factor = sfactor;
        if (factor == 1 && sfactor > 1)
        {
            factor++;
        }

        for (; factor < np2+1; factor++)
        {
            if (np2%factor == 0)
            {
                break;
            }
        }
        hdecomp[1] = factor;
        hdecomp[2] = np2/factor;

        decomposePar.subDict("hierarchicalCoeffs").add("n", hdecomp, true);


        decomposePar.regIOobject::write();

    }


    Info<< nl << "End." << endl;

    return(0);
}


// ************************************************************************* //
