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

Class
    externalFvPatchScalarField

Description
    Reads boundary field values from a specified location

    Example of the boundary condition specification:
    @verbatim
    externalBoundary
    {
        type            external;
        value           uniform 0;  // placeholder

        //entries below are optional, default values provided
        freqFormat        FILE; // FIXED|FILE
        freqFile        "OpenFOAM.ctrl";
        freqValue         10; //will not be read by default
        commsDir        "constant/boundaryData"; //directory to read from
        commsRoot        "/users/home/case"; //default case directory
        fileFormat        MONORAW; //IOFIELD|RAW|MONORAW
    }
    @endverbatim


    Formats supported are:
    IOFIELD - OpenFOAM <scalar>IOField
            - default input filename is <fileName>.<boundaryName>
    RAW - text file, space/newline delimited, uses boundary face ordering
            - default input filename is <fileName>.<boundaryName>
    MONORAW - As RAW, but single text file for all MONORAW external boundaries
              Boundary region ordering uses foam patch ordering
            - default input filename is <fileName>.boundaries

    IO conventions:
    - Read is delayed while <inputfile>.done does not exists in comms directory
    - <inputfile>.done is deleted by externalBoundary upon read completion

SourceFiles
    externalFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef externalFvPatchScalarField_H
#define externalFvPatchScalarField_H

#include "primitives/random/Random/Random.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class externalFvPatch Declaration
\*---------------------------------------------------------------------------*/

class externalFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
public:

    // Public data types
        enum IFormats
        {
            eIOFIELD,
            eRAW,
            eMONORAW
        };

        enum ufFormats
        {
            eFIXED,
            eFILE
        };

protected:

        static const NamedEnum<IFormats, 3> IFormatNames_;
        static const NamedEnum<ufFormats, 2> ufFormatNames_;

        static label monorawCounter;

    // Protected member functions

        //- Waits for file to reach specified state
        void fileState(const fileName&, bool exists) const;

private:

    // Private data

        //- Current time index
        label curTimeIndex_;

        //- How to determine update frequency, legal options are: FIXED and FILE
        ufFormats ufFormat_;

        //- Update frequency file name
        fileName ufFile_;

        //- Update frequency
        label updateFreq_;

        //- Communication directory
        fileName IODir_;

        //- Root of communications directory
        fileName IOroot_;

        //- this fields name
        word fieldName_;

        //- Input file format, legal options are: IOFIELD, RAW, MONORAW
        IFormats IFormat_;

        //- Input fileName
        fileName IFile_;


protected:
        //- monolithic input file start index
        label monoStart_;

        //- monolithic input file size
        label monoSize_;

        //- needed to set mono properties before first update
        bool firstUpdate_;

        //- parallel boundary face map
        List<labelList> parMap_;

private:

    // Private member functions

        //- update monolithic input file settings
        void getMonoProperties();

        //- update frequency if necessary
        void updateUpdateFreq();

        //- read external field file
        Field<scalar> readExternalField();

        //- assign field to boundary region
        void assignField(const Field<scalar>);

        //- build local to global boundary face maps for parallel runs
        void buildParMap();

public:

    //- Runtime type information
    TypeName("external");


    // Constructors

        //- Construct from patch and internal field
        externalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct by mapping given externalFvPatchScalarField
        //  onto a new patch
        externalFvPatchScalarField
        (
            const externalFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct from patch, internal field and dictionary
        externalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new externalFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        externalFvPatchScalarField
        (
            const externalFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (const DimensionedField<scalar, volMesh>& iF) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new externalFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return const access to input file format
            const IFormats& inputFormat() const
            {
                return IFormat_;
            }

            //return const access to IFile
            const fileName& IFile() const
            {
                return IFile_;
            }

            //- return const access to frequency format
            const ufFormats& updateFormat() const
            {
                return ufFormat_;
            }

            //- Return const access to update frequency
            const label& updateFrequency() const
            {
                return updateFreq_;
            }

            //- return const access to frequency file
            const fileName& frequencyFile() const
            {
                return ufFile_;
            }

            //- Return const access to IODir
            const fileName& IOdirectory() const
            {
                return IODir_;
            }

            //- Return non-const access to IODir
            fileName& IOdirectory()
            {
                return IODir_;
            }

            //- Return const access to IOroot
            const fileName& IOroot() const
            {
                return IOroot_;
            }

            //- Return const access to fieldName_
            const word& fieldName() const
            {
                return fieldName_;
            }

        // Edit
            //- Set input file name
            void setIFile();
            void setIFile(const fileName&);

            //- Change input format
            void setIOFormat(const word&);

            //- Modify update frequency
            void setUpdateFrequency(const label);
            void setUpdateFrequency(const fileName&);

            //- Change IOroot
            void IOroot(const fileName&);


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
