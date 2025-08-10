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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //
namespace Foam
{

    template<class Type>
    List<fileName> mappedProfileFvPatchField<Type>::inputFileNameList_(0);

    template<class Type>
    PtrList<pointProfileDataList>
        mappedProfileFvPatchField<Type>::inputData_(0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::mappedProfileFvPatchField<Type>::updateValues()
{
    this->forceAssign(mappedField_);
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::initialiseFluent
(
    const dictionary& dict
)
{
    pointField scaledBoundaryPoints(this->patch().Cf()/convertToMeters_);

    tmp<Field<Type>> tmappedField
    (
        mapProfileToField
        (
            scaledBoundaryPoints,
            fieldname_,
            regionName_
        )
    );

    if (tmappedField.valid())
    {
        mappedField_ = tmappedField();
    }

    //check availability of alternative field data
    if (dict.found("altDirection") && dict.found("altField"))
    {
        negSurfDotDir_.reset(new vector(dict.lookup("altDirection")));
        altFieldName_.reset(new word(dict.lookup("altField")));
    }

    if (negSurfDotDir_.valid() && altFieldName_.valid())
    {
        DynamicList<label> negSurfDotFaceIndices(this->patch().size());

        scalarField dotP(this->patch().Sf() & negSurfDotDir_());
        forAll(dotP, fI)
        {
            if (dotP[fI] < 0)
            {
                negSurfDotFaceIndices.append(fI);
            }
        }
        negSurfDotFaceIndices.shrink();

        //labelList altMap.reset(new labelList(negSurfDotFaceIndices));
        pointField altPoints(negSurfDotFaceIndices.size());

        forAll(negSurfDotFaceIndices, ii)
        {
            altPoints[ii]
                = this->patch().Cf()[negSurfDotFaceIndices[ii]]
                /convertToMeters_;
        }

        tmp<Field<Type>> altMappedField
        (
            mapProfileToField
            (
                altPoints,
                altFieldName_,
                regionName_
            )
        );

        if (altMappedField.valid())
        {
            forAll(altMappedField(), ii)
            {
                mappedField_[negSurfDotFaceIndices[ii]] = altMappedField()[ii];
            }
        }
    }
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::initialiseCSV
(
    const dictionary& dict
)
{
    inputFileNameList_.setSize(1);
    dict.lookup("file") >> inputFileNameList_[0];

    if (!this->size())
    {
        return;
    }

    pointField scaledBPoints(this->patch().Cf()/convertToMeters_);

    fileName fName = proFile_;
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fName.expand()));
    ISstream& is = isPtr();

    if (!is.good())
    {
        FatalIOErrorInFunction(is) << "Problem with CSV file."
            << exit(FatalIOError);
    }


    Switch mergeSeparators =
        dict.lookupOrDefault<Switch>("mergeSeparators", false);
    char separator =
        dict.lookupOrDefault<string>("separator", string(","))[0];

    //- Header read
    string headerLine;
    is.getLine(headerLine);

    List<string> header
    (
        readLineElements(headerLine, separator, mergeSeparators)
    );

    //- find column maps
    labelList locB(3, label(-1));
    labelList locField(pTraits<Type>::nComponents, label(-1));

    if (!coorFramePtr_)
    {
        forAll(header, sI)
        {
            if (header[sI] == "x") locB[0] = sI;
            if (header[sI] == "y") locB[1] = sI;
            if (header[sI] == "z") locB[2] = sI;
        }
        forAll(header, sI)
        {
            if (locField.size()!=1)
            {
                if (header[sI] == fieldname_+"_x") locField[0] = sI;
                if (header[sI] == fieldname_+"_y") locField[1] = sI;
                if (header[sI] == fieldname_+"_z") locField[2] = sI;
            }
            else
            {
                if (header[sI] == fieldname_) locField[0] = sI;
            }
        }
    }
    else
    {
        forAll(header, sI)
        {
            if (coorFramePtr_->coorSys().type()=="cartesian")
            {
                if (header[sI] == "x") locB[0] = sI;
                if (header[sI] == "y") locB[1] = sI;
                if (header[sI] == "z") locB[2] = sI;
            }
            else if (coorFramePtr_->coorSys().type()=="cylindrical")
            {
                if (header[sI] == "r") locB[0] = sI;
                if (header[sI] == "theta") locB[1] = sI;
                if (header[sI] == "z") locB[2] = sI;
            }
        }
        forAll(header, sI)
        {
            if (locField.size()!=1)
            {
                if (coorFramePtr_->coorSys().type()=="cartesian")
                {
                    if (header[sI] == fieldname_+"_x") locField[0] = sI;
                    if (header[sI] == fieldname_+"_y") locField[1] = sI;
                    if (header[sI] == fieldname_+"_z") locField[2] = sI;
                }
                else if (coorFramePtr_->coorSys().type()=="cylindrical")
                {
                    if (header[sI] == fieldname_+"_r") locField[0] = sI;
                    if (header[sI] == fieldname_+"_theta") locField[1] = sI;
                    if (header[sI] == fieldname_+"_z") locField[2] = sI;
                }
            }
            else
            {
                if (header[sI] == fieldname_) locField[0] = sI;
            }
        }
    }

    forAll(locB, cI)
    {
        if (locB[cI]==-1)
        {
            FatalErrorInFunction
                << "Error in defining the columns of the point coordiates"
                << this->internalField().name() << tab << locB
                << " of the csv file."
                << exit(FatalError);
        }
    }
    forAll(locField, cI)
    {
        if (locField[cI]==-1)
        {
            FatalErrorInFunction
                << "Error in defining the columns of the field "
                << this->internalField().name() << tab << locField
                << " of the csv file."
                << exit(FatalError);
        }
    }


    //- read data
    DynamicList<point> boundaryPoints;
    DynamicList<Type> field;

    while (is.good())
    {
        string line;
        is.getLine(line);

        List<string> splitted
        (
            readLineElements(line, separator, mergeSeparators)
        );
        if (splitted.size() <= 1) break;

        //- coordinate
        {
            scalar p0 = readScalar(IStringStream(splitted[locB[0]])());
            scalar p1 = readScalar(IStringStream(splitted[locB[1]])());
            scalar p2 = readScalar(IStringStream(splitted[locB[2]])());
            boundaryPoints.append(point(p0, p1, p2));
        }

        //- field
        setTypeValue(field, locField, splitted);
    }
    boundaryPoints.shrink();
    field.shrink();


    //- transformations
    pointField transformedBP(boundaryPoints);

    if (coorFramePtr_)
    {
        transformedBP = coorFramePtr_->coorSys().globalPosition(transformedBP);
        forAll(field, pI)
        {
            transformTypeValue(field[pI], transformedBP[pI]);
        }
    }

    //- calculate mapping
    labelList map(scaledBPoints.size(), label(-1));
    calculateMap(map, transformedBP, scaledBPoints);

    //- map and assign field
    label nFailedMaps(0);

    mappedField_ = Field<Type>(this->size(), Zero);
    forAll(map, pI)
    {
        if (map[pI]==-1) nFailedMaps++;
        mappedField_[pI] = field[map[pI]];
    }
    if (nFailedMaps>0)
    {
        WarningInFunction << "Find nearest map failed. "
            << "Check source points again."
            << endl;
    }
}


template<class Type>
Foam::List<Foam::string> Foam::mappedProfileFvPatchField<Type>::readLineElements
(
    const string& line,
    const char& separator,
    const Switch& mergeSeparators
)
{
    //- code copied from CSV.C and modified
    label n = 0;
    std::size_t pos = 0;
    const label nEntries = 100; // max for sanity reading

    DynamicList<string> splitted;

    if (mergeSeparators)
    {
        std::size_t nPos = 0;

        while ((pos != std::string::npos) && (n <= nEntries))
        {
            bool found = false;
            while (!found)
            {
                nPos = line.find(separator, pos);

                if ((nPos != std::string::npos) && (nPos - pos == 0))
                {
                    pos = nPos + 1;
                }
                else
                {
                    found = true;
                }
            }

            nPos = line.find(separator, pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos = nPos;
                n++;
            }
            else
            {
                splitted.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
                n++;
            }
        }
    }
    else
    {
        while ((pos != std::string::npos) && (n <= nEntries))
        {
            std::size_t nPos = line.find(separator, pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos = nPos;
                n++;
            }
            else
            {
                splitted.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
                n++;
            }
        }
    }
    return List<string>(splitted, true);
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::setTypeValue
(
    DynamicList<Type>& field,
    const List<label>& locField,
    const List<string>& splitted
)
{
    Type fieldElement(Zero);
    forAll(locField, cI)
    {
        scalar p0 = readScalar(IStringStream(splitted[locField[cI]])());
        scalar& gieldElementI = fieldElement.component(cI);
        gieldElementI = p0;
    }
    field.append(fieldElement);
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::transformTypeValue
(
    Type& data,
    const vector& global
)
{
      data = coorFramePtr_->coorSys().R(global)&data;
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::calculateMap
(
    labelList& map,
    const Field<point>& sPoints,
    const Field<point>& tPoints
)
{
    //construct octree
    treeBoundBox overallBb(tPoints);
    treeBoundBox overallBbs(sPoints);
    overallBb.min() = min(overallBb.min(), overallBbs.min());
    overallBb.max() = max(overallBb.max(), overallBbs.max());

    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    indexedOctree<treeDataPoint> pointTree
    (
        treeDataPoint(sPoints),
        overallBb,  // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    //loop over points, find nearest and assign mapped data
    scalar spanSqr = Foam::sqr(pointTree.bb().mag());

    forAll(tPoints, pI)
    {
        const point& p = tPoints[pI];

        pointIndexHit info = pointTree.findNearest
        (
            p,
            spanSqr
        );

        if (!info.hit())
        {
            info = pointTree.findNearest
            (
                p,
                GREAT*spanSqr
            );
        }

        map[pI] = info.index();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedProfileFvPatchField<Type>::mapProfileToField
(
    const pointField& targetPoints,
    word fieldName,
    word regionName
)
{
    FatalErrorInFunction
        << "Function only defined for scalar and vector types."
        << exit(FatalError);

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            this->patch().size(), pTraits<Type>::zero
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::mappedProfileFvPatchField<Type>::mappedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    proFile_(""),
    regionName_(""),
    fieldname_(""),
    convertToMeters_(0),
    negSurfDotDir_(),
    altFieldName_(),
    mappedField_(p.size(), Zero),
    frameName_(""),
    coorFramePtr_(nullptr)
{
}


template<class Type>
Foam::mappedProfileFvPatchField<Type>::mappedProfileFvPatchField
(
    const mappedProfileFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    proFile_(ptf.proFile_),
    regionName_(ptf.regionName_),
    fieldname_(ptf.fieldname_),
    convertToMeters_(ptf.convertToMeters_),
    negSurfDotDir_(),
    altFieldName_(),
    mappedField_(mapper(ptf.mappedField_)),
    frameName_(ptf.frameName_),
    coorFramePtr_(ptf.coorFramePtr_)
{
    if (ptf.negSurfDotDir_.valid() && ptf.altFieldName_.valid())
    {
        negSurfDotDir_.reset(new vector(ptf.negSurfDotDir_()));
        altFieldName_.reset(new word(ptf.altFieldName_()));
    }
}


template<class Type>
Foam::mappedProfileFvPatchField<Type>::mappedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    proFile_(dict.lookup("file")),
    regionName_(""),
    fieldname_(""),
    convertToMeters_(dict.lookupOrDefault<scalar>("convertToMeters", 1)),
    negSurfDotDir_(),
    altFieldName_(),
    mappedField_(p.size(), Zero),
    frameName_(dict.lookupOrDefault<word>("referenceFrame", "")),
    coorFramePtr_(nullptr)
{
    if (!coorFramePtr_ && frameName_ != word::null)
    {
        coorFramePtr_ = &coordinateFrame::New(iF.mesh(), frameName_);
    }
    if (proFile_.ext()!= "csv")
    {
        regionName_  = word((dict.lookup("region")));
        fieldname_  = word(dict.lookup("fieldname"));
        initialiseFluent(dict);
    }
    else if (proFile_.ext()== "csv")
    {
        fieldname_  =
        (
            dict.lookupOrDefault<word>
            (
                "field",
                this->internalField().name()
            )
        );
        initialiseCSV(dict);
    }
    else
    {
        FatalErrorInFunction
            << "Only prof and CSV based inputs are supported."
            << exit(FatalError);
    }
}


template<class Type>
Foam::mappedProfileFvPatchField<Type>::mappedProfileFvPatchField
(
    const mappedProfileFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    proFile_(ptf.proFile_),
    regionName_(ptf.regionName_),
    fieldname_(ptf.fieldname_),
    convertToMeters_(ptf.convertToMeters_),
    negSurfDotDir_(),
    altFieldName_(),
    mappedField_(ptf.mappedField_),
    frameName_(ptf.frameName_),
    coorFramePtr_(ptf.coorFramePtr_)
{
    if (ptf.negSurfDotDir_.valid() && ptf.altFieldName_.valid())
    {
        negSurfDotDir_.reset(new vector(ptf.negSurfDotDir_()));
        altFieldName_.reset(new word(ptf.altFieldName_()));
    }
}


template<class Type>
Foam::mappedProfileFvPatchField<Type>::mappedProfileFvPatchField
(
    const mappedProfileFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    proFile_(ptf.proFile_),
    regionName_(ptf.regionName_),
    fieldname_(ptf.fieldname_),
    convertToMeters_(ptf.convertToMeters_),
    negSurfDotDir_(),
    altFieldName_(),
    mappedField_(ptf.mappedField_),
    frameName_(ptf.frameName_),
    coorFramePtr_(ptf.coorFramePtr_)
{
    if (ptf.negSurfDotDir_.valid() && ptf.altFieldName_.valid())
    {
        negSurfDotDir_.reset(new vector(ptf.negSurfDotDir_()));
        altFieldName_.reset(new word(ptf.altFieldName_()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    m(mappedField_, mappedField_);
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
    const mappedProfileFvPatchField<Type>& dmptf =
        refCast<const mappedProfileFvPatchField<Type>>(ptf);
    mappedField_.rmap(dmptf.mappedField_, addr);
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::autoMapGIB(mapper);
    mapper.map(mappedField_, pTraits<Type>::zero);
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->updateValues();

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedProfileFvPatchField<Type>::write(Ostream& os) const
{
    this->fvPatchField<Type>::write(os);
    os.writeEntry("file", proFile_);
    this->template writeEntryIfDifferent<word>
    (
        os, "region", word::null, regionName_
    );
    this->template writeEntryIfDifferent<word>
    (
        os, "field", this->internalField().name(), fieldname_
    );

    this->template writeEntryIfDifferent<scalar>
    (
        os, "convertToMeters", 1.0, convertToMeters_
    );

    this->template writeEntryIfDifferent<word>
    (
        os, "referenceFrame", word::null, frameName_
    );

    if (negSurfDotDir_.valid() && altFieldName_.valid())
    {
        os.writeEntry("altDirection", negSurfDotDir_());
        os.writeEntry("altField", altFieldName_());

    }

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
