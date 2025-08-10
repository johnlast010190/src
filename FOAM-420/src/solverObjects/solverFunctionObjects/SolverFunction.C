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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "SolverFunction.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::functionObjects::SolverFunction<Type>::SolverFunction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    fvMeshFunctionObject(name, runTime, dict),
    solverObject_(name, const_cast<fvMesh&>(this->mesh_), dict),
    postProcess_(dict.lookupOrDefault<Switch>("postProcess", false))
{
    if (debug)
    {
        Info<< "call to SolverFunction Constructor" << endl;
    }
    read(dict);
    solverObject_.initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::functionObjects::SolverFunction<Type>::~SolverFunction()
{
    if (debug)
    {
        Info<< "call to SolverFunction Destructor" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::SolverFunction<Type>::read(const dictionary& dict)
{
    solverObject_.readBase(dict);
    return true;
}


template<class Type>
bool Foam::functionObjects::SolverFunction<Type>::execute()
{
    if (!postProcess_)
    {
        Info<< "functionObject:"  << type() << ":" << name() << endl;
        // Run solver
        solverObject_.solve();
    }
    return true;
}


template<class Type>
bool Foam::functionObjects::SolverFunction<Type>::write()
{
    if (!postProcess_)
    {
        solverObject_.write();
    }
    return true;
}


template<class Type>
bool Foam::functionObjects::SolverFunction<Type>::end()
{
    const dictionary& dict = solverObject_.dict();
    fvMesh& mesh = const_cast<fvMesh&>(mesh_);

    if (dict.found("solutionControl"))
    {
        Info<< "functionObject:"  << type() << ":" << name() << endl;
        Info<< "Executing post-processing" << endl;

        const word solutionControl = word(dict.lookup("solutionControl"));
        Info<< "solutionControl: " << solutionControl << endl;

        // reset time and time index to start values
        Time& time = const_cast<Time&>(mesh.time());
        time.setTime(solverObject_.timeStart(), solverObject_.timeStart());

        // call write function of the solver
        solverObject_.write();

        if (solutionControl == "SIMPLE" || solutionControl == "simple" || solutionControl == "Simple")
        {
            simpleControl simpleC(mesh);
            scalar nIter(solverObject_.timeEnd() - solverObject_.timeStart());
            label nonOrthCorr(simpleC.nNonOrthCorr());

            // hard-coded time loop
            // cannot use simle.loop() as it inlcudes functionObjects
            // and thus creates a loop in a loop ...
            for (label i=0; i<nIter; i++)
            {
                time++;
                Info<< "Time = " << time.timeName() << nl << endl;

                for (label j=0; j<=nonOrthCorr; j++)
                {
                    // run the solver
                    solverObject_.solve();
                }
            }
        }
        else
        {
            WarningInFunction
                << "solutionControl method not supported, doing nothing!"
                << endl;
        }
    }
    return true;
}

// ************************************************************************* //
