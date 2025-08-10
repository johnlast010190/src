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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::word Foam::referenceFrameFvPatch<Type>::loadUName
(
    const DimensionedField<Type, volMesh>& iF
) const
{
    if (isNull(iF) || !isA<DimensionedField<vector, volMesh>>(iF))
    {
        return "U";
    }
    return iF.name();
}


template<class Type>
const Foam::objectRegistry& Foam::referenceFrameFvPatch<Type>::loadDb
(
    const DimensionedField<Type, volMesh>& iF
) const
{
    if (notNull(iF))
    {
        return iF.db();
    }
    else
    {
        return patch_.boundaryMesh().mesh();
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::referenceFrameFvPatch<Type>::setCoorFramePtr() const
{
    if (!coorFramePtr_ && frameName_ != word::null)
    {
        const fvMesh& mesh = patch_.boundaryMesh().mesh();
        coorFramePtr_ = &coordinateFrame::New(mesh, frameName_);
        if (!inletFlux_)
        {
            coorFramePtr_->attachPatch(patch_.index());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::referenceFrameFvPatch<Type>::referenceFrameFvPatch
(
    const fvPatch& patch,
    const DimensionedField<Type, volMesh>& iF
)
:
    coorFramePtr_(nullptr),
    patch_(patch),
    inputValue_(patch.size(), Zero),
    inletFlux_(false),
    frameName_(),
    refPatchUName_(loadUName(iF)),
    obr_(loadDb(iF)),
    definedInFrame_(false)
{}


template<class Type>
Foam::referenceFrameFvPatch<Type>::referenceFrameFvPatch
(
    const fvPatch& patch,
    coordinateFrame* coorFramePtr,
    const Field<Type>& inputValue,
    Switch inletFlux,
    const word& frameName,
    const DimensionedField<Type, volMesh>& iF
)
:
    coorFramePtr_(coorFramePtr),
    patch_(patch),
    inputValue_(inputValue),
    inletFlux_(inletFlux),
    frameName_(frameName),
    refPatchUName_(loadUName(iF)),
    obr_(loadDb(iF)),
    definedInFrame_(false)
{}


template<class Type>
Foam::referenceFrameFvPatch<Type>::referenceFrameFvPatch
(
    const fvPatch& patch,
    coordinateFrame* coorFramePtr,
    const Field<Type>& inputValue,
    Switch inletFlux,
    const word& frameName,
    const fvPatchFieldMapper& mapper,
    const DimensionedField<Type, volMesh>& iF
)
:
    coorFramePtr_(coorFramePtr),
    patch_(patch),
    inputValue_(patch.size(), Zero),
    inletFlux_(inletFlux),
    frameName_(frameName),
    refPatchUName_(loadUName(iF)),
    obr_(loadDb(iF)),
    definedInFrame_(false)
{
    mapper(inputValue_, inputValue);
}


template<class Type>
Foam::referenceFrameFvPatch<Type>::referenceFrameFvPatch
(
    const dictionary& dict,
    const fvPatch& patch,
    const DimensionedField<Type, volMesh>& iF
)
:
    coorFramePtr_(nullptr),
    patch_(patch),
    inputValue_(patch.size(), Zero),
    inletFlux_(dict.lookupOrDefault<Switch>("inletFlux", false)),
    frameName_(dict.lookupOrDefault<word>("referenceFrame", "")),
    refPatchUName_(loadUName(iF)),
    obr_(loadDb(iF)),
    definedInFrame_(dict.lookupOrDefault<Switch>("definedInFrame", false))
{
    if (dict.found("inputValue"))
    {
        inputValue_ = Field<Type>("inputValue", dict, patch.size());
    }
    if
    (
        isDefinedInFrame()
     && coorFramePtr()
    && !coorFramePtr()->coorSys().uniform()
    )
    {
        FatalErrorInFunction
            << "Option \"definedInFrame\" isn't supported for non-uniform"
            << " coordinate system yet."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::referenceFrameFvPatch<Type>::autoMap(const fvPatchFieldMapper& m)
{
    m(inputValue_, inputValue_);
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::rmap
(
    const referenceFrameFvPatch<Type>& ptf,
    const labelList& addr
)
{
    const referenceFrameFvPatch<Type>& trptf =
        dynamic_cast<const referenceFrameFvPatch<Type>&>(ptf);

    inputValue_.rmap(trptf.inputValue_, addr);
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mapper.map(inputValue_, pTraits<Type>::zero);
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::addFrameVelocity
(
    Field<Type>& Up,
    bool setInValue,
    bool setToGlobal
)
{}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::addFrameGradient
(
    Field<Type>& gradient
) const
{}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::makeVectorGlobal
(
    Field<Type>& value
) const
{}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::makeVectorGlobal
(
    Type& value
) const
{}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::makeRelative(Field<scalar>& phiP)
{
    if (coorFramePtr())
    {
        forAll(phiP, facei)
        {
            phiP[facei] -=
                coorFramePtr()->frameVelocity(patch_.Cf()[facei])
               &patch_.Sf()[facei];
        }
    }
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::makeRelative(vectorField& Up)
{
    if (coorFramePtr())
    {
        Up -= coorFramePtr()->frameVelocity(patch_.Cf());
    }
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::updateCordinateFrameRegistry() const
{
    if (!inletFlux_ && coorFramePtr())
    {
        coorFramePtr()->attachPatch(patch_.index());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::referenceFrameFvPatch<Type>::getFrameVelocity() const
{
    return tmp<Field<Type>>(new Field<Type>(patch_.size(), Zero));
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::read(const dictionary& dict)
{
    inletFlux_ = dict.lookupOrDefault<Switch>("inletFlux", false);
    frameName_ = dict.lookupOrDefault<word>("referenceFrame", "");
    definedInFrame_ = dict.lookupOrDefault<Switch>("definedInFrame", false);
    if
    (
        isDefinedInFrame()
     && coorFramePtr()
     &&!coorFramePtr()->coorSys().uniform()
    )
    {
        FatalErrorInFunction
            << "Option \"definedInFrame\" isn't supported for non-uniform"
            << " coordinate system yet."
            << abort(FatalError);
    }
    if (dict.found("inputValue"))
    {
        inputValue_ = Field<Type>("inputValue", dict, patch_.size());
    }
    else
    {
        inputValue_ = Zero;
    }
}


template<class Type>
void Foam::referenceFrameFvPatch<Type>::write(Ostream& os) const
{
    inputValue_.writeEntry("inputValue", os);
    if (inletFlux_)
    {
        os.writeEntry("inletFlux", inletFlux_);
    }
    if (frameName_ != word::null)
    {
        os.writeEntry("referenceFrame", frameName_);
        if (definedInFrame_)
        {
            os.writeEntry("definedInFrame", definedInFrame_);
        }
    }
}


// ************************************************************************* //
