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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "ccnvReport/ccnvReport.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ccnvReport, 0);
    addToRunTimeSelectionTable(functionObject, ccnvReport, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ccnvReport::ccnvReport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    ccnv_pr_(8),
    UName_("U"),
    PhiName_("phi"),
    PName_("p"),
    stepNumber_(0),
    logToFile_(false),
    logFileName_("ccnv.log"),
    logFilePtr_(nullptr),
    maxDivergence_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ccnvReport::~ccnvReport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ccnvReport::execute()
{
    Info<< type() << " " << name() <<  " write:" << nl;

    //calculate
    calculate();
    stepNumber_++;

    const dictionary& solverDict = mesh_.solverPerformanceDict();

    double pressureRes = 1;
    double pressureIters = 1;
    double momentumRes = 1;
    double momentumIters = 1;

    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();

        if (variableName == UName_ || variableName ==PName_)
        {
            const List<solverPerformance> sp(iter().stream());

            if (variableName == PName_)
            {
                pressureRes = sp.first().initialResidual();
                pressureIters = sp.first().nIterations();
            } else if (variableName == UName_)
            {
                momentumRes = sp.first().initialResidual();
                momentumIters = sp.first().nIterations();
            }
        }
    }

    logFilePtr_.open
    (
        logFileName_.c_str(),
        std::ios::out |
        std::ios::app |
        std::ios::binary
     );

    double actualTime = obr_.time().value();

    //write to file
    if (logToFile_ && (Pstream::master() || !Pstream::parRun()))
    {
        double snd = stepNumber_;
        Info<<"    Step Number: "<<snd<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&snd), ccnv_pr_);

        Info<<"    Time: "<<actualTime<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&actualTime), ccnv_pr_);

        Info<<"    Maximum Divergence: "<<maxDivergence_<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&maxDivergence_), ccnv_pr_);

        Info<<"    Residual of Pressure: "<<pressureRes<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&pressureRes), ccnv_pr_);

        Info<<"    Iterations of Pressure Eq: "<<pressureIters<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&pressureIters), ccnv_pr_);

        Info<<"    Residual of Momentum: "<<momentumRes<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&momentumRes), ccnv_pr_);

        Info<<"    Iterations of Momentum Eq: "<<momentumIters<<endl;
        logFilePtr_.write( reinterpret_cast<char*>(&momentumIters), ccnv_pr_);
    }

    logFilePtr_.close();

    Info<< endl;

    return true;
}


void Foam::functionObjects::ccnvReport::calculate()
{
    //- find the max divergence
    if (foundObject<surfaceScalarField>(PhiName_))
    {
        const surfaceScalarField& phi
            (lookupObject<surfaceScalarField>(PhiName_));
        tmp<volScalarField> divPhi(fvc::div(phi));
        maxDivergence_ = gMax(divPhi().internalField());
    }
    else
    {
        WarningInFunction << "Field " << PhiName_
                << " could not be found in the object registry."
                << " No ccnvReport compiled." << endl;
    }

}

bool Foam::functionObjects::ccnvReport::write()
{
    return true;
}


bool Foam::functionObjects::ccnvReport::read(const dictionary& dict)
{
    //- CCNV ID
    long ccnv_id = 64408111;

    //- CCNV Version
    char ccnv_ver[16] = "1.0.0";

    //- account ID (future versions)
    char ccnv_actid[128]="";

    //- job ID (future versions)
    char ccnv_jobid[128]="";

    //- password (future versions)
    char ccnv_pass[128]="";

    //- CCNV Solver ID
    char ccnvSolverID[128]="OS-0044-0004-000000-0000-0000-0000";

    //- number of x-axis titles
    long num_x=2;

    //- x-axis titles
    char x_title[2][128]=
    {
        {"STEP NO."},
        {"TIME"}
    };

    //- number of y-axis titles
    long num_y=5;

    //- y-axis titles
    char y_title[5][128]=
    {
        {"MAXIMUM DIVERGENCE"},
        {"RESIDUAL OF PRESSURE EQ"},
        {"ITERATION OF PRESSURE EQ"},
        {"RESIDUAL OF MOMENTUM EQ"},
        {"ITERATION OF MOMENTUM EQ"}
    };

    ccnv_pr_ = dict.lookupOrDefault<int>("writePrecision",8);

    //get field names
    UName_ = dict.lookupOrDefault<word>("UName","U");
    PhiName_ = dict.lookupOrDefault<word>("PhiName","phi");
    PName_ = dict.lookupOrDefault<word>("PName","p");

    logToFile_= dict.lookupOrDefault<bool>("logToFile", true);

    //make sure file directories exist and create OFstream
    if (logToFile_ && (Pstream::master() || !Pstream::parRun()))
    {
        if (Pstream::parRun())
        {
            logFileName_
                =  obr_.rootPath()/obr_.caseName()/".."/"postProcessing"
                /name()/logFileName_;
        }
        else
        {
            logFileName_
                = obr_.rootPath()/obr_.caseName()/"postProcessing"
                /name()/logFileName_;
        }

        Info<< "    Logging CCNV information to file: " << logFileName_
             << endl;

        //create file and dir if necessary
        if (Pstream::parRun())
        {
            if
            (
                !Foam::isDir
                (obr_.rootPath()/obr_.caseName()/".."/"postProcessing"/name())
            )
            {
                Foam::mkDir(obr_.rootPath()/obr_.caseName()/".."
                            /"postProcessing"/name());
            }
        }
        else
        {
            if
            (
                !Foam::isDir
                (obr_.rootPath()/obr_.caseName()/"postProcessing"/name())
            )
            {
                Foam::mkDir(obr_.rootPath()/obr_.caseName()
                            /"postProcessing"/name());
            }
        }

        logFilePtr_.open(logFileName_.c_str(), std::ios::out | std::ios::binary);

        //- write header information
        logFilePtr_.write( reinterpret_cast<char*>(&ccnv_id), 4);
        logFilePtr_.write( reinterpret_cast<char*>(&ccnv_ver), 16);
        logFilePtr_.write( reinterpret_cast<char*>(&ccnv_actid), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&ccnv_jobid), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&ccnv_pass), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&ccnvSolverID), 128);

        logFilePtr_.write( reinterpret_cast<char*>(&num_x), 4);
        logFilePtr_.write( reinterpret_cast<char*>(&x_title[0]), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&x_title[1]), 128);

        logFilePtr_.write( reinterpret_cast<char*>(&num_y), 4);
        logFilePtr_.write( reinterpret_cast<char*>(&y_title[0]), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&y_title[1]), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&y_title[2]), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&y_title[3]), 128);
        logFilePtr_.write( reinterpret_cast<char*>(&y_title[4]), 128);

        logFilePtr_.write( reinterpret_cast<char*>(&ccnv_pr_), 4);

        logFilePtr_.close();
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
