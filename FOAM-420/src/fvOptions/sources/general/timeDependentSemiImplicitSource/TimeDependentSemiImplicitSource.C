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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "TimeDependentSemiImplicitSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "finiteVolume/fvm/fvmSup.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::fv::TimeDependentSemiImplicitSource<Type>::volumeModeTypeNames_
(
    IStringStream("(absolute specific)")()
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fv::TimeDependentSemiImplicitSource<Type>::sourceMatrix
(
    const dimensionSet& dSet,
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const label fieldi
) const
{
    if (debug)
    {
        Info<< "TimeDependentSemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    DimensionedField<Type, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            dSet/dimVolume,
            Zero
        ),
        false
    );

    if (SuTimeSeriesList_(fieldi) != nullptr)
    {
        UIndirectList<Type>(Su, cells_) =
            SuTimeSeriesList_[fieldi].value
            (mesh_.time().timeOutputValue())/VDash_;
    }

    DimensionedField<scalar, volMesh> Sp
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>
        (
            "zero",
            Su.dimensions()/psi.dimensions(),
            0.0
        ),
        false
    );

    if (SpTimeSeriesList_(fieldi) != nullptr)
    {
        UIndirectList<scalar>(Sp, cells_) =
            SpTimeSeriesList_[fieldi].value
            (mesh_.time().timeOutputValue())/VDash_;
    }

    return Su + fvm::SuSp(Sp, psi);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<class Type>
bool Foam::fv::TimeDependentSemiImplicitSource<Type>::wordListFind(const wordList& wl, const word& w)
{
    forAll(wl, i)
    {
        if (wl[i] == w)
        {
            return true;
        }
    }

    return false;
}

template<class Type>
void Foam::fv::TimeDependentSemiImplicitSource<Type>::setTimeSeries(const dictionary& dict)
{
    dictionary SuRateDict(dict.subDict("injectionRateSu"));
    dictionary SpRateDict(dict.subDict("injectionRateSp"));

    // check Su and Sp lists for non-unique field names
    {
        // add all unique field names to filedNames
        // toc and iterators seem smart enough not
        // to use duplicate
        fieldNames_.setSize(SuRateDict.toc().size());
        label i = 0;

        forAllConstIter(dictionary, SuRateDict, iter)
        {
            fieldNames_[i] = iter().keyword();
            i++;
        }

        wordList SpFields;
        SpFields.setSize(SpRateDict.toc().size());
        label j = 0;

        // grab all unique field names in SpRateDcit
        forAllConstIter(dictionary, SpRateDict, iter)
        {
            SpFields[j] = iter().keyword();
            j++;
        }

        // merge the two to make one unique list
        forAll(SpFields,fI)
        {
            if (!wordListFind(fieldNames_, SpFields[fI]))
            {
                fieldNames_.append(SpFields[fI]);
            }
        }
    }

    SuTimeSeriesList_.setSize(fieldNames_.size());
    SpTimeSeriesList_.setSize(fieldNames_.size());

    Info<<"    Reading Su and Sp semi-implicit source time series data"<<endl;

    forAll(fieldNames_,i)
    {
        fileName lookupFileName;

        if (SuRateDict.found(fieldNames_[i]))
        {
            SuTimeSeriesList_.set
            (
                i,
                Function1<Type>::New(fieldNames_[i], SuRateDict)
            );
            Info<< "      Found Field "<< fieldNames_[i]
                 << " Su time series data"
                 << endl;

        } else
        {
            SuTimeSeriesList_.set
            (
                i,
                nullptr
            );

            Info<< "      Field "<< fieldNames_[i]
                 << " does not have Su time series data, assumed as zero. "
                 << endl;

        }

        if (SpRateDict.found(fieldNames_[i]))
        {

            SpTimeSeriesList_.set
            (
                i,
                Function1<scalar>::New(fieldNames_[i], SpRateDict)
            );

            Info<< "      Found Field "<< fieldNames_[i]
                 << " Sp time series data"
                 << endl;

        } else
        {
            SpTimeSeriesList_.set
            (
                i,
                nullptr
            );

            Info<< "      Field "<< fieldNames_[i]
                 << " does not have Sp time series data, assumed as zero. "
                 << endl;
        }
    }

    applied_.setSize(fieldNames_.size(), false);

    // Set volume normalisation
    if (volumeMode_ == vmAbsolute)
    {
        VDash_ = V_;
    }
}

template<class Type>
typename Foam::fv::TimeDependentSemiImplicitSource<Type>::volumeModeType
Foam::fv::TimeDependentSemiImplicitSource<Type>::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


template<class Type>
Foam::word Foam::fv::TimeDependentSemiImplicitSource<Type>::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::TimeDependentSemiImplicitSource<Type>::TimeDependentSemiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
    SuTimeSeriesList_(),
    SpTimeSeriesList_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::TimeDependentSemiImplicitSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    eqn += this->sourceMatrix(eqn.dimensions(), eqn.psi(), fieldi);
}


template<class Type>
void Foam::fv::TimeDependentSemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "TimeDependentSemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldi);
}


template<class Type>
void Foam::fv::TimeDependentSemiImplicitSource<Type>::addSup
(
    fvBlockMatrix<Type>& eqn,
    const label fieldi
)
{
    eqn += this->sourceMatrix(eqn.dimensionSets()[0], eqn.psi(), fieldi);
}


template<class Type>
void Foam::fv::TimeDependentSemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldi);
}


// ************************************************************************* //
