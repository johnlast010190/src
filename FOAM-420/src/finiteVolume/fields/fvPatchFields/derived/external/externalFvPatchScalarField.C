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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/external/externalFvPatchScalarField.H"
#include "db/Time/Time.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "db/IOobjects/IOField/IOField.H"
#include "volMesh/volMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::externalFvPatchScalarField::IFormats, 3>::names[] =
{
    "IOFIELD",
    "RAW",
    "MONORAW"
};

//template<>
const Foam::NamedEnum<Foam::externalFvPatchScalarField::IFormats, 3>
    Foam::externalFvPatchScalarField::IFormatNames_;

template<>
const char* Foam::NamedEnum<Foam::externalFvPatchScalarField::ufFormats, 2>
::names[] =
{
    "FIXED",
    "FILE"
};

//template<>
const Foam::NamedEnum<Foam::externalFvPatchScalarField::ufFormats, 2>
    Foam::externalFvPatchScalarField::ufFormatNames_;


label externalFvPatchScalarField::monorawCounter = 0;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void externalFvPatchScalarField::fileState(const fileName& f, bool shouldExist) const
{
    label tooLong = 1;
    label iters = 0;

    if (shouldExist)
    {
        while (!isFile(f))
        {
            sleep(1);
            iters++;
            if (iters > tooLong)
            {
                iters = 0;
                Pout<< "Waiting for file: "
                     << f << " to be created." << endl;
            }
        }
    }
    else
    {
        while (isFile(f))
        {
            sleep(1);
            iters++;
            if (iters > tooLong)
            {
                iters = 0;
                Pout<< "Waiting for file: "
                     << f << " to be deleted." << endl;
            }
        }

    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void externalFvPatchScalarField::getMonoProperties()
{

    const volScalarField& field
        = this->db().lookupObject<volScalarField>(fieldName_);

    const fvPatchList& patches = this->patch().boundaryMesh();

    forAll(field.boundaryField(), patchI)
    {

        if
        (
            field.boundaryField()[patchI].type()
            == externalFvPatchScalarField::typeName
        )
        {
            IFormats patchIFormat
                = (refCast<const externalFvPatchScalarField>
                (field.boundaryField()[patchI])).inputFormat();



            if (patchIFormat == eMONORAW)
            {
                label pSize = patches[patchI].size();

                if (&(patches[patchI].patch()) == &(this->patch().patch()))
                {
                    monoStart_ = monoSize_;
                }

                if (Pstream::parRun())
                {
                    reduce(pSize, sumOp<label>());
                }
                monoSize_ += pSize;
            }
        }
    }
}

void externalFvPatchScalarField::updateUpdateFreq()
{

    switch(ufFormat_)
    {
        case eFIXED:
            //
            break;

        case eFILE:

            fileState(fileName(IOroot_/IODir_/ufFile_),true);

            fileState
            (
                fileName(IOroot_/IODir_/ufFile_+".busy"),false
            );

            {
                IFstream infile(IOroot_/IODir_/ufFile_);
                infile.read(updateFreq_);
                if (updateFreq_ < 1 || updateFreq_%1 != 0)
                {
                    FatalError << "Invalid or non-existant update frequency: "
                               << updateFreq_
                               << " specified in control file: "
                               << ufFile_ << endl
                               << "Update frequency must be a positive "
                               << "integer number."
                               << exit(FatalError);
                }
            }

            break;

        default:

            FatalError << ufFormatNames_[ufFormat_] << " is not a valid "
                       << "update frequency input format." << endl
                       << "Valid formats are: FIXED and FILE."
                       << exit(FatalError);

    }


}

Field<scalar> externalFvPatchScalarField::readExternalField()
{
    label pSize = this->size();

    if (Pstream::parRun())
    {
        reduce(pSize, sumOp<label>());
    }

    Field<scalar> newField(0, 0.0);


    if (Pstream::master() || !Pstream::parRun())
    {
        newField.setSize(pSize);
        newField = 0;
        fileState(fileName(IOroot_/IODir_/IFile_),true);

        switch(IFormat_)
        {
            case eIOFIELD:
            {
                fileState
                (
                    fileName(IOroot_/IODir_/IFile_+ string(".done")),
                    true
                );
                newField = IOField<scalar>
                (
                    IOobject
                    (
                        IFile_,
                        IODir_,
                        this->patch().boundaryMesh().mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
                rm(IOroot_/IODir_/IFile_+ string(".done"));
                break;
            }

            case eRAW:
            {
                fileState
                (
                    fileName(IOroot_/IODir_/IFile_+ string(".done")),
                    true
                );
                IFstream infile
                (
                    IOroot_/IODir_/IFile_
                );

                forAll(newField, fI)
                {
                    if (infile.good())
                    {
                        scalar var = -GREAT;
                        infile>>var;
                        if (var != -GREAT)
                        {
                            newField[fI] = var;
                        }
                        else
                        {
                            FatalError << "Blank line in file: "
                                       << IFile_
                                       << " at position : " << fI + 1
                                       << exit(FatalError);

                        }

                    }
                    else
                    {
                        FatalError << "Read beyond end of file: "
                                   << IFile_ << exit(FatalError);
                    }
                }
                rm(IOroot_/IODir_/IFile_+ string(".done"));
            }
            break;

            case eMONORAW:
            {
                if (pSize > 0)
                {
                    fileState
                    (
                        fileName(IOroot_/IODir_/IFile_+ string(".done")),
                        true
                    );
                    IFstream infile
                    (
                        IOroot_/IODir_/IFile_
                    );
                    //scroll through entries to start of this boundary
                    for (label i = 0; i < monoStart_; i++)
                    {
                        infile.read(newField[0]);
                    }

                    //read this boundary
                    forAll(newField, fI)
                    {
                        if (infile.good())
                        {
                            scalar var = -GREAT;
                            infile>>var;
                            if (var != -GREAT)
                            {
                                newField[fI] = var;
                            }
                            else
                            {
                                FatalError << "Blank line in file: "
                                           << IFile_
                                           << " at position : "
                                           << monoStart_ + fI + 1
                                           << exit(FatalError);
                            }
                        }
                        else
                        {
                            FatalError << "Read beyond end of file: "
                                       << IFile_ << exit(FatalError);
                        }
                    }
                    externalFvPatchScalarField::monorawCounter += pSize;

                    if (externalFvPatchScalarField::monorawCounter == monoSize_)
                    {
                        rm(IOroot_/IODir_/IFile_+ string(".done"));
                        externalFvPatchScalarField::monorawCounter = 0;
                    }
                }
            }
            break;

            default:
                FatalError << IFormatNames_[IFormat_] << " is not a valid external "
                           << "boundary input file format." << endl
                           << "Valid formats are: IOFIELD, RAW, MONORAW."
                           << exit(FatalError);
        }
    }

    //synchronise
    Pstream::barrier();

    return newField;
}

void externalFvPatchScalarField::assignField(const Field<scalar> nf)
{
    if (Pstream::parRun())
    {
        List<scalarList> psib(Pstream::nProcs());

        psib[Pstream::myProcNo()].setSize(size());
        psib[Pstream::myProcNo()] = -1;

        Pstream::gatherList(psib);

        if (Pstream::master())
        {
            forAll(psib, procI)
            {
                forAll(psib[procI], bfI)
                {
                    label fromFace = parMap_[procI][bfI];
                    psib[procI][bfI] = nf[fromFace];
                }
            }
        }

        forAll(psib, procI)
        {
            Pstream::scatter(psib[procI]);
        }

        Field<scalar>& patchField = *this;

        forAll(psib[Pstream::myProcNo()], bfI)
        {
            patchField[bfI] = psib[Pstream::myProcNo()][bfI];
        }
    }
    else
    {
        Field<scalar>& patchField = *this;
        patchField = nf;
    }
}

void externalFvPatchScalarField::buildParMap()
{
    if (Pstream::parRun())
    {
        parMap_.setSize(Pstream::nProcs());

        label mpn = Pstream::myProcNo();

        parMap_[mpn].setSize(this->size());

        const fileName& inst
            = this->patch().boundaryMesh().mesh().polyMesh::instance();

        //read faceProc addressing (check whether faceProcAddressing is
        // available in "inst". If not set "inst" to "constant")

        autoPtr<IOobject> fPA_header
        (
            new IOobject
            (
                "faceProcAddressing",
                inst,
                this->patch().boundaryMesh().mesh().meshSubDir,
                this->patch().boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (!fPA_header->typeHeaderOk<labelIOList>(true))
        {
            fPA_header.reset
            (
                new IOobject
                (
                    "faceProcAddressing",
                    this->patch().boundaryMesh().mesh().time().constant(),
                    this->patch().boundaryMesh().mesh().meshSubDir,
                    this->patch().boundaryMesh().mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
        }

        labelIOList faceProcAddress(fPA_header());

        //map global face index to local patch index
        forAll(parMap_[mpn], mfI)
        {
            label currentFace = this->patch().patch().start() + mfI;

            parMap_[mpn][mfI] = mag(faceProcAddress[currentFace]) -1;
        }


        //gather maps to master
        Pstream::gatherList(parMap_);

        //make global face indexing local to current boundary
        //by subtracting the smallest face index from all entries
        //i.e. startFace = 0

        label minFace = this->patch().boundaryMesh().mesh().nFaces();
        reduce(minFace, sumOp<label>());

        if (Pstream::master())
        {
            forAll(parMap_, procI)
            {
                if (parMap_[procI].size() != 0)
                {
                    minFace = min(minFace, min(parMap_[procI]));
                }
            }

            forAll(parMap_, procI)
            {
                if (parMap_[procI].size() != 0)
                {
                    parMap_[procI] = parMap_[procI] - minFace;
                }
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

externalFvPatchScalarField::externalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    curTimeIndex_(-1),
    ufFormat_(eFIXED),
    ufFile_(fileName::null),
    updateFreq_(1),
    IODir_(fileName::null),
    IOroot_(fileName::null),
    fieldName_(word::null),
    IFormat_(eIOFIELD),
    IFile_(fileName::null),
    monoStart_(0),
    monoSize_(0),
    firstUpdate_(true),
    parMap_(0)
{
}


externalFvPatchScalarField::externalFvPatchScalarField
(
    const externalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    curTimeIndex_(-1),
    ufFormat_(ptf.ufFormat_),
    ufFile_(ptf.ufFile_),
    updateFreq_(ptf.updateFreq_),
    IODir_(ptf.IODir_),
    IOroot_(ptf.IOroot_),
    fieldName_(ptf.fieldName_),
    IFormat_(ptf.IFormat_),
    IFile_(ptf.IFile_),
    monoStart_(ptf.monoStart_),
    monoSize_(ptf.monoSize_),
    firstUpdate_(ptf.firstUpdate_),
    parMap_(ptf.parMap_)
{
}


externalFvPatchScalarField::externalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    curTimeIndex_(-1),
    ufFormat_(ufFormatNames_.read(IStringStream(string("FILE"))())),
    ufFile_("OpenFOAM.ctrl"),
    updateFreq_(-1),
    IODir_(fileName("AURA")),
    IOroot_(fileName::null),
    fieldName_(this->internalField().name()),
    IFormat_(IFormatNames_.read(IStringStream(string("MONORAW"))())),
    IFile_(fileName::null),
    monoStart_(0),
    monoSize_(0),
    firstUpdate_(true),
    parMap_(0)
{

    if (Pstream::parRun())
    {
        IOroot(this->db().time().path()/"..");
    }
    else
    {
        IOroot(this->db().time().path());
    }
    if (dict.found("commsRoot"))
    {
        IOroot(fileName(dict.lookup("commsRoot")));
    }

    if (dict.found("fileFormat"))
    {
        IFormat_ = (IFormatNames_.read(dict.lookup("fileFormat")));
    }
    if (dict.found("commsDir"))
    {
        IODir_ = (fileName(dict.lookup("commsDir")));
    }
    if (dict.found("freqFormat"))
    {
        ufFormat_ = ufFormatNames_.read(dict.lookup("freqFormat"));
    }
    switch(ufFormat_)
    {
        case eFIXED:
            updateFreq_ = readLabel(dict.lookup("freqValue"));
            break;

        case eFILE:
            if (dict.found("freqFile"))
            {
                ufFile_ = fileName(dict.lookup("freqFile"));
            }
            //updateUpdateFreq();
            break;

        default:
            FatalError << ufFormatNames_[ufFormat_] << " is not a valid update "
                       << "frequency input format." << endl
                       << "Valid formats are: FIXED and FILE."
                       << exit(FatalError);
    }


    setIFile();
}


externalFvPatchScalarField::externalFvPatchScalarField
(
    const externalFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    curTimeIndex_(-1),
    ufFormat_(ptf.ufFormat_),
    ufFile_(ptf.ufFile_),
    updateFreq_(ptf.updateFreq_),
    IODir_(ptf.IODir_),
    IOroot_(ptf.IOroot_),
    fieldName_(ptf.fieldName_),
    IFormat_(ptf.IFormat_),
    IFile_(ptf.IFile_),
    monoStart_(ptf.monoStart_),
    monoSize_(ptf.monoSize_),
    firstUpdate_(ptf.firstUpdate_),
    parMap_(ptf.parMap_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Edit
void externalFvPatchScalarField::setIFile()
{
    switch(IFormat_)
    {
        case eIOFIELD:
            IFile_ = fileName(fieldName_+"."+this->patch().name());
            break;

        case eRAW:
            IFile_ = fileName(fieldName_+"."+this->patch().name());
            break;

        case eMONORAW:
            IFile_ = fileName(fieldName_+".boundaries");
            break;

        default:
            FatalError << IFormatNames_[IFormat_] << " is not a valid external "
                       << "boundary input file format." << endl
                       << "Valid formats are: IOFIELD, RAW, MONORAW."
                       << exit(FatalError);
    }
}


void externalFvPatchScalarField::setIFile(const fileName& ifile)
{
    IFile_ = ifile;
}

void externalFvPatchScalarField::setIOFormat(const word& IOf)
{
    IFormat_ = IFormatNames_[IOf];
}

void externalFvPatchScalarField::setUpdateFrequency(const label f)
{
    ufFormat_ = eFIXED;
    updateFreq_ = f;
}
void externalFvPatchScalarField::setUpdateFrequency(const fileName& uff)
{
    ufFile_ = uff;
    ufFormat_ = eFILE;
    updateUpdateFreq();
}

void externalFvPatchScalarField::IOroot(const fileName& ir)
{
    if (IFormat_ != eIOFIELD)
    {
        IOroot_ = ir;
    }
    else
    {
        IOroot_ = this->db().time().path();

        Warning << "Non-const access to IOroot_ is not allowed "
                << "for IOFIELD file format." << endl
                << "IOroot_ : " << IOroot_ << endl;
    }
}

// Map from self
void externalFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void externalFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<scalar>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void externalFvPatchScalarField::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    if (firstUpdate_)
    {
        updateUpdateFreq();

        if (IFile_ == fileName::null)
        {
            setIFile();
        }

        getMonoProperties();
        buildParMap();

        curTimeIndex_ = 0;
        updateUpdateFreq();
        assignField(readExternalField());

        firstUpdate_ = false;
    }
    else
    {
        if (curTimeIndex_ != this->db().time().timeIndex())
        {
            curTimeIndex_ = this->db().time().timeIndex();

            //first change update frequency if required
            updateUpdateFreq();

            if (curTimeIndex_%updateFreq_ == 0)
            {
                assignField(readExternalField());
            }
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


// Write
void externalFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("commsDir", IOdirectory());

    bool patchMatch = false;

    fileName eIOroot = IOroot();
    eIOroot.expand();

    if (Pstream::parRun())
    {
        fileName cPath = this->db().time().path()/"..";
        cPath.expand();

        if (cPath == eIOroot)
        {
            patchMatch = true;
        }
    }
    else
    {
        fileName cPath = this->db().time().path();
        cPath.expand();
        if (cPath == eIOroot)
        {
            patchMatch = true;
        }
    }

    if (!patchMatch)
    {
        os.writeEntry("commsRoot", IOroot());
    }

    os.writeEntry
    (
        "fileFormat",
        externalFvPatchScalarField::IFormatNames_[inputFormat()]
    );

    os.writeEntry
    (
        "freqFormat",
        externalFvPatchScalarField::ufFormatNames_[updateFormat()]
    );

    switch(updateFormat())
    {
        case eFIXED:
        {
            os.writeEntry("freqValue", updateFrequency());
        }
        break;

        case eFILE:
        {
            os.writeEntry("freqFile", frequencyFile());
            break;
        }
        default:
        {
            FatalError << ufFormatNames_[updateFormat()] << " is not a valid update "
                       << "frequency input format." << endl
                       << "Valid formats are: FIXED and FILE."
                       << exit(FatalError);
        }
    }

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, externalFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
