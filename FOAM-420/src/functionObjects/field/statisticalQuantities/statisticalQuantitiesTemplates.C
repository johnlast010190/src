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
    (c) 2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::statisticalQuantities::pop
(
    const DLList<Type>& runFieldList
)
{
    bool pop = false;
    if (sAver2)
    {
        pop = runFieldList.size()>2*sampleSize_;
    }
    else
    {
        pop = runFieldList.size()>sampleSize_;
    }
    return pop;
}


template<class Type>
void Foam::functionObjects::statisticalQuantities::calc
(
    List<DLList<Type>>& runFieldLists,
    const labelList& map
)
{
    forAll(runFieldLists, fieldi)
    {
        /*
        if (runTimeLists_.size()<runFieldLists[fieldi].size())
        {
            Info<< "Incompatible value and dt sizes for: "
                 << fieldNames_[map[fieldi]]
                 << " . Skipping moving average for that field."
                 << endl;
            return;
        }
        */

        const word& fieldName(fieldNames_[map[fieldi]]);

        const word valueType = objectResultType(functionObjectName_, fieldName);

        if (pTraits<Type>::typeName != valueType)
        {
            return;
        }

        Type currentValue =
            getObjectResult<Type>(functionObjectName_, fieldName);

        runFieldLists[fieldi].append(currentValue);
        if (pop(runFieldLists[fieldi]))
        {
            runFieldLists[fieldi].removeHead();
        }

        const word meanName(fieldName + "MovMean");
        Type meanValue = Zero;
        if (sAver)
        {
            typename DLList<Type>::const_reverse_iterator fieldIter =
                runFieldLists[fieldi].rbegin();

            //DLList<scalar>::const_reverse_iterator dtIter =
            //    runTimeLists_.rbegin();

            label i=0;
            for
            (
                fieldIter = runFieldLists[fieldi].rbegin();
                (fieldIter != runFieldLists[fieldi].rend()) && i< sampleSize_;
                ++fieldIter, ++i
            )
            {
                meanValue += (*fieldIter);
            }
            meanValue /= i;

            setResult(meanName, meanValue);

            file() << tab << meanValue;

            Log<< "    " << meanName << ": " << meanValue << nl;

            if (sDevAver)
            {
                runFieldMeanScalarLists_[fieldi].append(meanValue);
                if (pop(runFieldMeanScalarLists_[fieldi]))
                {
                    runFieldMeanScalarLists_[fieldi].removeHead();
                }
                this->dict_.set(meanName, runFieldMeanScalarLists_[fieldi]);
            }
        }
        const word mean2Name(fieldName + "MovMean2");
        Type mean2Value = Zero;
        if (sAver2)
        {
            typename DLList<Type>::const_reverse_iterator fieldIter =
                runFieldLists[fieldi].rbegin();

            for
            (
                fieldIter = runFieldLists[fieldi].rbegin();
                (fieldIter != runFieldLists[fieldi].rend());
                ++fieldIter
            )
            {
                mean2Value += (*fieldIter);

            }
            mean2Value /= runFieldLists[fieldi].size();

            setResult(mean2Name, mean2Value);

            file() << tab << mean2Value;

            Log<< "    " << mean2Name << ": " << mean2Value << nl;
        }
        Type sigmaValue = Zero;
        if (sDev)
        {
            const word movStdDevName(fieldName + "MovStdDev");
            typename DLList<Type>::const_reverse_iterator fieldIter =
                runFieldLists[fieldi].rbegin();

            label i=0;
            for
            (
                fieldIter = runFieldLists[fieldi].rbegin();
                (fieldIter != runFieldLists[fieldi].rend()) && i< sampleSize_;
                ++fieldIter, ++i
            )
            {
                sigmaValue += sqr((*fieldIter)-meanValue);
            }
            sigmaValue = sqrt(sigmaValue/i);

            setResult(movStdDevName, sigmaValue);

            file() << tab << sigmaValue;

            Log<< "    " << movStdDevName << ": " << sigmaValue << nl;
        }
        Type sigmaMeanValue = Zero;
        Type meanMeanValue = Zero;
        if (sDevAver)
        {
            const word movStdDevAverName(fieldName + "MovStdDevAver");
            {
                typename DLList<Type>::const_reverse_iterator fieldIter =
                    runFieldMeanScalarLists_[fieldi].rbegin();

                label i=0;
                for
                (
                    fieldIter = runFieldMeanScalarLists_[fieldi].rbegin();
                    (fieldIter != runFieldMeanScalarLists_[fieldi].rend()) &&
                    i< sampleSize_;
                    ++fieldIter, ++i
                )
                {
                    meanMeanValue += (*fieldIter);
                }
                meanMeanValue /= i;
            }
            {
                typename DLList<Type>::const_reverse_iterator fieldIter =
                    runFieldMeanScalarLists_[fieldi].rbegin();

                label i=0;
                for
                (
                    fieldIter = runFieldMeanScalarLists_[fieldi].rbegin();
                    (fieldIter != runFieldMeanScalarLists_[fieldi].rend()) &&
                    i< sampleSize_;
                    ++fieldIter, ++i
                )
                {
                    sigmaMeanValue += sqr((*fieldIter)-meanMeanValue);
                }
                sigmaMeanValue = sqrt
                    (
                        sigmaMeanValue/i
                    );
            }

            setResult(movStdDevAverName, sigmaMeanValue);

            file() << tab << sigmaMeanValue;

            Log<< "    " << movStdDevAverName << ": " << sigmaMeanValue << nl;
        }
        Type covValue = Zero;
        if (sCov)
        {
            const word movCovName(fieldName + "MovCoV");
            if (meanValue==pTraits<Type>::zero)
            {
                 Warning << "meanValue is zero for the CoV computation. "
                         << "Setting CoV to zero."
                         << endl;
            }
            else
            {
                covValue = sigmaValue/meanValue;
            }

            setResult(movCovName, covValue);

            file() << tab << covValue;

            Log<< "    " << movCovName << ": " << covValue << nl;
        }

        if (sCoDav)
        {
            const word movCoDaVName(fieldName + "MovCoDaV");
            Type coDaVValue = Zero;

            if (meanValue==pTraits<Type>::zero)
            {
                 Warning << "meanValue is zero for the CoDaV computation. "
                         << "Setting CoDaV to zero."
                         << endl;
            }
            else
            {
                coDaVValue = Foam::mag(scalar(meanValue-mean2Value))/meanValue;
            }

            setResult(movCoDaVName, coDaVValue);

            file() << tab << coDaVValue;

            Log<< "    " << movCoDaVName << ": " << coDaVValue << nl;
        }

        if (sCoVav)
        {
            const word movCoVaVName(fieldName + "MovCoVaV");
            Type coVaVValue = Zero;

            if (meanValue==pTraits<Type>::zero)
            {
                 Warning << "meanValue is zero for the CoVaV computation. "
                         << "Setting CoVaV to zero."
                         << endl;
            }
            else
            {
                coVaVValue = sigmaMeanValue/meanValue;
            }

            setResult(movCoVaVName, coVaVValue);

            file() << tab << coVaVValue;

            Log<< "    " << movCoVaVName << ": " << coVaVValue << nl;
        }

        this->dict_.set(fieldName, runFieldLists[fieldi]);
    }
}


template<class Type>
void Foam::functionObjects::statisticalQuantities::populate
(
    List<DLList<Type>>& fieldTypeList,
    labelList& map
)
{
    label nType(0);
    DynamicList<label> dynList(fieldNames_.size());

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName(fieldNames_[fieldi]);
        const word valueType =
            this->objectResultType(functionObjectName_, fieldName);
        if (pTraits<Type>::typeName == valueType)
        {
            dynList.append(fieldi);
            nType++;
        }
        else
        {
             Warning << "Quantity " << fieldNames_[fieldi]
                     << " was not found in " <<  functionObjectName_
                     << " FO data."
                     << endl;
        }
    }
    dynList.shrink();
    map = dynList;

    fieldTypeList = List<DLList<Type>>(nType);

    forAll(map, i)
    {
        const label mapI = map[i];
        const word& fieldName(fieldNames_[mapI]);
        dict_.readIfPresent(fieldName, fieldTypeList[i]);
    }
}


// ************************************************************************* //
