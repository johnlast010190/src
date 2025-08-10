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
    (c) 2017-2021 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "TimeDependentBoundarySource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "primitives/strings/lists/wordReList.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::fv::TimeDependentBoundarySource<Type>::surfaceModeTypeNames_
(
    IStringStream("(absolute specific)")()
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<class Type>
void Foam::fv::TimeDependentBoundarySource<Type>::setTimeSeries(const dictionary& dict)
{
    labelList patchIDs;
    wordReList patches(dict.lookup("patches"));
    wordRes patchMatch(patches);
    forAll(mesh().boundary(), i)
    {
        if (patchMatch.match(mesh().boundary()[i].name()))
        {
            patchIDs.append(i);
        }
    }

    dictionary fluxDict(dict.subDict("fluxSource"));

    wordList fieldNames(fluxDict.toc().size());
    label i = 0;

    forAllConstIter(dictionary, fluxDict, iter)
    {
        fieldNames[i] = iter().keyword();
        i++;
    }

    Info<<"    Reading boundary source time series data"<<endl;

    forAll(fieldNames, i)
    {
        timeSeriesList_.set
        (
            fieldNames[i],
            Function1<Type>::New(fieldNames[i], fluxDict)
        );
    }

    forAll(fieldNames, i)
    {
        boundarySourcePatchIDs_.set(fieldNames[i], patchIDs);
        boundaryApplied_.set(fieldNames[i], boolList(patchIDs.size(), false));
    }

    // Set volume normalisation
    if (surfaceMode_ == smAbsolute)
    {
        areaNorm_ = scalar(0);
        forAll(patchIDs, i)
        {
            areaNorm_ += sum(mesh().boundary()[patchIDs[i]].magSf());
        }
        reduce(areaNorm_, sumOp<scalar>());
        Info<< "   Boundary source spread over area: " << areaNorm_ << endl;
    }
}

template<class Type>
typename Foam::fv::TimeDependentBoundarySource<Type>::surfaceModeType
Foam::fv::TimeDependentBoundarySource<Type>::wordToSurfaceModeType
(
    const word& smtName
) const
{
    forAll(surfaceModeTypeNames_, i)
    {
        if (smtName == surfaceModeTypeNames_[i])
        {
            return surfaceModeType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown surfaceMode type " << smtName
        << ". Valid surfaceMode types are:" << nl << surfaceModeTypeNames_
        << exit(FatalError);

    return surfaceModeType(0);
}


template<class Type>
Foam::word Foam::fv::TimeDependentBoundarySource<Type>::surfaceModeTypeToWord
(
    const surfaceModeType& vmtType
) const
{
    if (vmtType > surfaceModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return surfaceModeTypeNames_[vmtType];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::TimeDependentBoundarySource<Type>::TimeDependentBoundarySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(name, modelType, dict, obr),
    surfaceMode_(smSpecific),
    areaNorm_(1.0),
    timeSeriesList_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::TimeDependentBoundarySource<Type>::addBoundarySource
(
    const word& fieldName,
    const label patchID,
    const Field<Type>& pf,
    Field<Type>& f,
    Field<Type>& df
)
{
    if (debug)
    {
        Info<< "TimeDependentBoundarySource<" << pTraits<Type>::typeName
            << ">::addBoundarySource for source " << name_ << endl;
    }

    f +=
        timeSeriesList_[fieldName]().value
        (
            mesh_.time().timeOutputValue()
        )/areaNorm_;
}


// ************************************************************************* //
