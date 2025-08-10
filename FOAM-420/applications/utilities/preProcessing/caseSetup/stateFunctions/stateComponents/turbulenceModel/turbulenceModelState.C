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

#include "turbulenceModelState.H"
#include "db/dictionary/functionEntries/includeEtcEntry/includeEtcEntry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(turbulenceModelState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

primitiveEntry turbulenceModelState::divScheme
(
    const word& fieldName, const dictionary& dict
) const
{
    return
    {
        word("div(phi," + fieldName + ")"),
        dict.lookup("divSchemes")
    };
}

void turbulenceModelState::addTurbulenceSpecificSchemes
(
    const dictionary& turbDict, dictionary& fvSchemes
)
{
    const dictionary& specificDivSchemes
    (
        turbDict.subDict("settings").subDict("divSchemes")
    );
    fvSchemes.subDict("divSchemes").merge(specificDivSchemes, true);
    const dictionary& specificLaplacianSchemes
    (
        turbDict.subDict("settings").subDict("laplacianSchemes")
    );
    fvSchemes.subDict("laplacianSchemes").merge(specificLaplacianSchemes, true);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModelState::turbulenceModelState
(
    const word& region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    const word& meshName
)
:
    regionState(region, input, defaults, master, index, meshName)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void turbulenceModelState::fvSolutionInject
(
    const dictionary& fieldDict,
    const dictionary& modelDict,
    dictionary& fvSolDict,
    const word& tfname
)
{

    fvSolDict.subDict("relaxationFactors").subDict("equations").add
    (
        primitiveEntry
        (
            tfname,
            modelDict.lookup("relaxationFactor")
        ),
        false
    );

    //store existing entry (could be user defined state)
    dictionary prev;
    if (fvSolDict.subDict("solvers").found(tfname))
    {
        prev = fvSolDict.subDict("solvers").subDict(tfname);
        fvSolDict.subDict("solvers").remove(tfname);
    }

    primitiveEntry tol("vTol", modelDict.lookup("vTol"));
    primitiveEntry reltol
    (
        "vRelTol",
        modelDict.lookup("vRelTol")
    );

    fvSolDict.subDict("solvers").add(tfname, dictionary());
    dictionary& solverDict
    (
        fvSolDict.subDict("solvers").subDict(tfname)
    );
    solverDict.add(tol);
    solverDict.add(reltol);

    //cfg file
    fileName solverCfg
    (
        "caseDicts/preProcessing/caseSetup/settings/matrixSolvers/"
        + word(fieldDict.lookup("solver")) + ".cfg"
    );

    functionEntries::includeEtcEntry::mergeFile(solverDict, solverCfg);

    //overlay potential content from custom state
    fvSolDict.subDict("solvers").subDict(tfname).merge(prev);

    if (matrix() == maCoupled)
    {
        if (modelDict.lookup("maxIter"))
        {
            primitiveEntry maxIter
            (
                "maxIter",
                modelDict.lookup("maxIter")
            );
            solverDict.add(maxIter);
        }
    }

    //remove variables
    fvSolDict.subDict("solvers").subDict(tfname).remove("vTol");
    fvSolDict.subDict("solvers").subDict(tfname).remove("vRelTol");
}

void turbulenceModelState::initialise()
{
    //set turbulentProperties dictionary
    constant().add("turbulenceProperties", dictionary());
    dictionary& turbProp
    (
        constant().subDict("turbulenceProperties")
    );

    word turbCat(turbulenceTypeNames_[turbulence()]);

    turbProp.add("simulationType", turbCat);


    //get turbulence model defaults
    const dictionary& turbDB(defaults().subDict("turbulenceProperties"));

    word turbModelType
    (
        input().lookupOrDefault<word>("turbulenceModel", "laminar")
    );

    //workaround for laminar that doesnt require redoing everything
    //for singular model category
    if (turbModelType == "laminar")
    {
        turbProp.add("simulationType", turbModelType, true);
    }


    if (turbulence() != tuLam && turbModelType != "laminar")
    {
        //add to fieldMaps and inject turbulenceProperties dictionary

        word turbModelCategory
        (
            (usf() ? "USF" : word(compressibilityTypeNames_[compressibility()]))
            + word(turbulenceTypeNames_[turbulence()])
        );

        dictionary turbModelDefaults
        (
            turbDB.subDict(turbModelCategory).subDict(turbModelType)
        );

        //update fieldMaps
        if (turbModelDefaults.found("fieldMaps"))
        {
            const dictionary& turbFieldMaps
            (
                turbModelDefaults.subDict("fieldMaps")
            );

            stateDict_.subDict("fieldMaps").merge(turbFieldMaps);

            turbulenceFields_.setSize(turbFieldMaps.size());

            label i = 0;
            forAllConstIter(dictionary, turbFieldMaps, iter)
            {
                turbulenceFields_[i++] = iter().keyword();
            }

            turbModelDefaults.remove("fieldMaps");
        }

        //inject default settings into the specific turb model dictionary
        turbProp.add
        (
            turbCat, turbDB.subDict("baseTurbTypeDict")
        );

        dictionary& turbPropType(turbProp.subDict(turbCat));

        turbPropType.merge(turbModelDefaults);


        //do not overwrite
        turbPropType.add(word(turbCat + word("Model")), turbModelType, false);

        // add delta coefs
        if (turbulence() == tuLES && turbPropType.found("delta"))
        {
            word deltaType(turbModelDefaults.lookup("delta"));
            turbPropType.add
            (
                word(deltaType + "Coeffs"),
                turbDB.subDict(deltaType + "Coeffs")
            );
        }

    }

    stateFunction::initialise();
}

void turbulenceModelState::correct()
{
    //skip turbulent field settings injection for laminar
    word turbModelType
    (
        input().lookupOrDefault<word>("turbulenceModel", "laminar")
    );

    if (turbModelType != "laminar")
    {

        //inject fvSchemes & fvSolution settings for turbulence fields
        word turbCat(turbulenceTypeNames_[turbulence()]);
        dictionary& turbCatDict
        (
            constant().subDict("turbulenceProperties").subDict(turbCat)
        );

        forAll(turbulenceFields_, i)
        {
            //
            word tfname(turbulenceFields_[i]);

            //get field definition
            const dictionary& turbFieldDict
            (
                fieldDefinitions().subDict(tfname)
            );

            //check if field is solved for (solver)
            if (turbFieldDict.found("solver"))
            {
                //get settings
                word timeTN(timeTypeNames_[time()]);
                word matrixTN(matrixTypeNames_[matrix()]);

                const dictionary& turbSettings
                (
                    turbCatDict.subDict("settings").
                        subDict(matrixTN).subDict(timeTN)
                );

                //fvSchemes
                dictionary& fvSchemes(system().subDict("fvSchemes"));

                fvSchemes.subDict("divSchemes").add
                (
                    divScheme(tfname, turbSettings), false
                );

                addTurbulenceSpecificSchemes(turbCatDict, fvSchemes);

                word lapScheme(word::null);

                if (compressibility() == ctIncomp)
                {
                    lapScheme = "laplacian(D"+tfname+"Eff,"+tfname+")";
                }
                else if (compressibility() == ctComp)
                {
                    lapScheme = "laplacian((rho*D"+tfname+"Eff),"+tfname+")";
                }
                fvSchemes.subDict("laplacianSchemes").add
                (
                    primitiveEntry
                    (
                        lapScheme,
                        turbSettings.lookup("laplacianSchemes")
                    ),
                    false
                );

                if (turbSettings.found("convectionGradient"))
                {
                    fvSchemes.subDict("gradSchemes").add
                    (
                        primitiveEntry
                        (
                            "turbulence",
                            turbSettings.lookup("convectionGradient")
                        ),
                        false
                    );
                }

                //field gradients
                if (turbSettings.found("fieldGradient"))
                {
                    fvSchemes.subDict("gradSchemes").add
                    (
                        primitiveEntry
                        (
                            word("grad("+tfname+")"),
                            turbSettings.lookup("fieldGradient")
                        ),
                        false
                    );
                }
                if (turbSettings.found("nRequired"))
                {
                    if (turbSettings.lookup<bool>("nRequired"))
                    {
                        fvSchemes.subDict("wallDist").add
                        (
                            "nRequired", turbSettings.lookup("nRequired")
                        );
                    }
                }
                //fvSolution
                dictionary& fvSolution(system().subDict("fvSolution"));

                fvSolutionInject
                (
                    turbFieldDict,
                    turbSettings,
                    fvSolution,
                    tfname
                );

                if (time() == ttTrans)
                {
                    fvSolutionInject
                    (
                        turbFieldDict,
                        turbSettings.subDict("final"),
                        fvSolution,
                        tfname+"Final"
                    );
                }
            }
        }

        //remove settings
        if (turbCatDict.found("settings"))
        {
            turbCatDict.remove("settings");
        }
    }
    stateFunction::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
