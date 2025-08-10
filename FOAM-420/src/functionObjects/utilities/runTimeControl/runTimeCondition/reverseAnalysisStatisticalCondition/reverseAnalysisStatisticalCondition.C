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
    (c) 2018, 2020 Esi Ltd
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "runTimeControl/runTimeCondition/reverseAnalysisStatisticalCondition/reverseAnalysisStatisticalCondition.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/functionObjects/stateFunctionObject/stateFunctionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(reverseAnalysisStatisticalCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        reverseAnalysisStatisticalCondition,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::functionObjects::runTimeControls
::reverseAnalysisStatisticalCondition::calc
(
    const word& entryName,
    const scalar& currentValue,
    bool& satisfied,
    bool& processed
)
{

    //Get latest value and append to data set
    DLList<scalar> boxValues;
    dictionary& dict = this->conditionDict().subDict(entryName);
    dict.readIfPresent("boxValues", boxValues);
    boxValues.append(currentValue);


    //Prepare necessary arrays
    const label pts = boxValues.size();
    scalarField avg(pts,0.0);
    scalarField fit(pts,0.0);
    scalarField std(pts,0.0);
    scalarField stdx(pts,0.0);
    scalarField stdInt(pts,0.0);
    scalarField solQual(pts,0.0);


    //Populate cumulative average, gradient and standard deviation arrays
    DLList<scalar>::const_reverse_iterator valueIter = boxValues.crbegin();
    avg[0] = valueIter();
    ++valueIter;

    if (pts > 1)
    {
        avg[1] = (avg[0] + valueIter())/2;
        std[1] = fabs(valueIter() - avg[0])/sqrt(2.0);
        stdx[1] = 1.0/sqrt(2.0);
        fit[1] = valueIter() - avg[0];
        ++valueIter;
    }

    for (label i = 2; i < pts; ++i, ++valueIter)
    {
        const label in = i - 1;
        const label ip = i + 1;

        avg[i] = (i * avg[in] + valueIter())/ip;

        scalar a = valueIter() - avg[in];
        std[i] = sqrt(in*std[in]*std[in]/i + a*a/ip);

        stdx[i] = sqrt(in*stdx[in]*stdx[in]/i + 0.25*ip);
    }


    //Populate fit array
    for (label i = 2; i < pts; ++i)
    {
        scalar sum = 0.0;
        valueIter = boxValues.crbegin();
        for (label j = 0; j <= i; ++j, ++valueIter)
        {
            sum += (j - 0.5*i)*(valueIter() - avg[i]);
        }
        fit[i] = sum/stdx[i]/stdx[i]/i;
    }


    //Populate std 'Integral' array
    for (label i = 1; i < pts; ++i)
    {
        stdInt[i]=0;
        for (label j = 1 ; j <= i; ++j)
        {
            stdInt[i] += min(std[i], std[j]);
        }
        if (std[i] != 0)
        {
            stdInt[i] /= (i + 1)*std[i];
        }
        else
        {
            stdInt[i] = 1;
        }
    }

    //Overwrite based on std tolerance
    if (pts > minBox_)
    {
        if (std[minBox_-1] < stdTol_)
        {
            for (label i = 0; i < pts; ++i)
            {
                if (std[i] < stdTol_)
                {
                    stdInt[i] = 1.0;
                }
                else
                {
                    break;
                }
            }
        }
    }


    //Populate solution Quality array
    scalar maxSolQual = 0.0;
    label count = 1;
    label ind = 0;
    for (label i = 1; i < pts; ++i)
    {
        solQual[i] = stdInt[i] - (i + 1)*fabs(fit[i]);
        //Find maximum solution Quality and store index
        if (solQual[i] > maxSolQual)
        {
            maxSolQual = solQual[i];
            ind = i;
        }
        //Count instances of solution Quality over user tolerance
        if (solQual[i] > solQualTol_)
        {
            count++;
        }
    }


    //Output
    Info<< "  Total data points: " << pts << nl;
    Info<< "  Average over " << ind + 1 << " pts: " << avg[ind] << nl;
    Info<< "  Best fit gradient: " << -(ind + 1)*fit[ind] << nl;
    Info<< "  Max Sol Quality: " << maxSolQual << nl;

    if (Pstream::master() || !Pstream::parRun())
    {
        file() << obr_.time().timeName() << tab
               << pts << tab
               << ind + 1 << tab
               <<  avg[ind] << tab
               <<  -(ind + 1)*fit[ind] << tab
               <<  maxSolQual
               << endl;
    }


    //Discarding points from the data set based on max solution quality achieved
    label i = pts - 1;

    while (maxSolQual - solQual[i] > discardRate_)
    {
        boxValues.removeHead();
        --i;
    }

    if (pts - i > 1)
    {
    	Info<< "  discarded " << pts - i -1 << " pts..." << nl;
    }


    //Judging convergence
    if (count < minBox_ || fabs((ind + 1)*fit[ind]) > fitTol_)
    {
        satisfied = false;
    }

    processed = true;

    // Store the state information for the next step
    dict.set("boxValues", boxValues);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls
::reverseAnalysisStatisticalCondition::reverseAnalysisStatisticalCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    writeFile(obr_, name, typeName, dict),
    functionObjectName_(dict.lookup("functionObject")),
    entryNames_(dict.lookup("fields")),
    minBox_(max(2,dict.lookupOrDefault<scalar>("minBox",2))),
    stdTol_(dict.lookupOrDefault<scalar>("stdTolerance", VSMALL)),
    solQualTol_(dict.lookupOrDefault<scalar>("solQualTolerance", 0.95)),
    fitTol_(dict.lookupOrDefault<scalar>("fitTolerance", VGREAT)),
    discardRate_(dict.lookupOrDefault<scalar>("discardRate", 0.5)),
    vectorComp_(dict.lookupOrDefault<bool>("vectorComp", "false")),
    vectorMag_(dict.lookupOrDefault<bool>("vectorMag", "false")),
    vectorSum_(dict.lookupOrDefault<bool>("vectorSum", "false")),
    starting_(true),
    resetOnRestart_(dict.lookupOrDefault<bool>("resetOnRestart","false"))
{
    //set up file IO
    if (Pstream::master() || !Pstream::parRun())
    {
        //add headers
        writeCommented(file(), "Time");
        writeDelimited(file(), "TotalDataPoints");
        writeDelimited(file(), "AveragePoints");
        writeDelimited(file(), "AverageValue");
        writeDelimited(file(), "BestFitGradient");
        writeDelimited(file(), "MaxSolQuality");
        file() <<endl;
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls
::reverseAnalysisStatisticalCondition::~reverseAnalysisStatisticalCondition()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls
::reverseAnalysisStatisticalCondition::apply()
{
    bool satisfied = true;

    if (!active_)
    {
        return satisfied;
    }

    DynamicList<label> unprocessedFields(entryNames_.size());

    dictionary& conditionDict = this->conditionDict();

    if (starting_)
    {
        starting_=false;

        forAll(entryNames_, fieldi)
        {
            const word& entryName = entryNames_[fieldi];

            const word valueType =
                state_.objectResultType(functionObjectName_, entryName);

            if (resetOnRestart_)
            {
                if (pTraits<scalar>::typeName == valueType)
                {
                    conditionDict.set(entryName, dictionary());
                }

                if (pTraits<vector>::typeName == valueType && vectorComp_)
                {
                    const vector tempVector =
                        state_.getObjectResult<vector>(functionObjectName_, entryName);

                    for (label i = 0; i<tempVector.size(); ++i)
                    {
                        word entryName2 = entryName + "_" + Foam::name(i);
                        conditionDict.set(entryName2, dictionary());
                    }
                }

                if (pTraits<vector>::typeName == valueType && vectorMag_)
                {
                    word entryName2 = entryName + "_Mag";
                    conditionDict.set(entryName2, dictionary());
                }

                if (pTraits<vector>::typeName == valueType && vectorSum_)
                {
                    word entryName2 = entryName + "_Sum";
                    conditionDict.set(entryName2, dictionary());
                }
            }
            else
            {
                if (pTraits<scalar>::typeName == valueType)
                {
                    conditionDict.set(entryName, dictionary());
                }

                if (pTraits<vector>::typeName == valueType && vectorComp_)
                {
                    const vector tempVector =
                        state_.getObjectResult<vector>(functionObjectName_, entryName);

                    for (label i = 0; i<tempVector.size(); ++i)
                    {
                        word entryName2 = entryName + "_" + Foam::name(i);
                        conditionDict.set(entryName2, dictionary());
                    }
                }

                if (pTraits<vector>::typeName == valueType && vectorMag_)
                {
                    word entryName2 = entryName + "_Mag";
                    conditionDict.set(entryName2, dictionary());
                }

                if (pTraits<vector>::typeName == valueType && vectorSum_)
                {
                    word entryName2 = entryName + "_Sum";
                    conditionDict.set(entryName2, dictionary());
                }
            }
        }
    }

    forAll(entryNames_, entryi)
    {
        const word& entryName(entryNames_[entryi]);

        bool processed = false;

        Log << nl << " " << type() << ": " << entryNames_[entryi] << nl;

        const word valueType =
            state_.objectResultType(functionObjectName_, entryName);

        scalar currentValue = 0.0;

        if (pTraits<scalar>::typeName == valueType)
        {
            currentValue +=
                state_.getObjectResult<scalar>(functionObjectName_, entryName);
            calc(entryName, currentValue, satisfied, processed);
        }

        if (pTraits<vector>::typeName == valueType && vectorComp_)
        {
            const vector tempVector =
                state_.getObjectResult<vector>(functionObjectName_, entryName);

            for (label i = 0; i<tempVector.size(); ++i)
            {

                word entryName2 = entryName + "_" + Foam::name(i);
                calc(entryName2, tempVector[i], satisfied, processed);
                Info<< nl;
            }
        }

        if (pTraits<vector>::typeName == valueType && vectorMag_)
        {
            const vector tempVector =
                state_.getObjectResult<vector>(functionObjectName_, entryName);

            for (label i = 0; i<tempVector.size(); ++i)
            {
                currentValue += tempVector[i]*tempVector[i];
            }

            calc(entryName + "_Mag", sqrt(currentValue), satisfied, processed);
        }

        if (pTraits<vector>::typeName == valueType && vectorSum_)
        {
            const vector tempVector =
                state_.getObjectResult<vector>(functionObjectName_, entryName);

            for (label i = 0; i<tempVector.size(); ++i)
            {
                currentValue += tempVector[i];
            }

            if (vectorMag_)
            {
                Info<< nl;
            }

            calc(entryName + "_Sum", currentValue, satisfied, processed);
        }

        if (!processed)
        {
            unprocessedFields.append(entryi);
        }
    }

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        forAll(unprocessedFields, i)
        {
            label entryi = unprocessedFields[i];
            Info<< " " << entryNames_[entryi] << nl;
        }
    }

    Log << endl;

    return satisfied;
}


void Foam::functionObjects::runTimeControls
::reverseAnalysisStatisticalCondition::write()
{
}


// ************************************************************************* //
