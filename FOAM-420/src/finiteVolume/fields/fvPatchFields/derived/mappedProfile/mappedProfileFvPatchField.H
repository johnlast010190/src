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
    (c) 2017 Esi Ltd.

Class
    Foam::mappedProfileFvPatchField

Description
    Foam::mappedProfileFvPatchField

     A fixed value boundary mapped from either
     i. Fluent compatible .prof input file
     using a 0th order proximity method.
     If a region or field is not present
     in the prof file, a warning will be printed and the provided value field
     will be used. The altDirection and altField entries are optional. If they
     are both present and the dot product of the altDirection and the surface
     normal (witch points out of the fluid domain) is NEGATIVE, then the faces
     so marked receive their value from the altField. This allows distinction
     between front and back of simple baffles even if they are in the same
     patch.

     ii. csv file. The coordinates and field from cloud points are read
     automatically based on the header names. For locations, the header must
     contain either (x y z) or (r theta z) for cartesian and cylindrical
     coordanates respectivelly.
     The field is read based on the name in the header. In case of vectors, the
     header must contain _x _y _z (or _r _theta _z)
     extensions to identify the column correctly.

    Example of the boundary condition specification:

    //- Fluent based
    @verbatim
    hotWall
    {
        type        mappedProfile;
        file        "profiles/temperature.prof";//relative to the case directory
        region      heat_shield;
        field       front_temp;
        convertToMeters 0.001;
        value uniform 300; //required entry

        //optional entries
        altDirection (0 0 1);
        altField     back_temp;

    }
    @endverbatim

    //- csv based
    inlet
    {
        type        mappedProfile;
        file        "0/data.csv";//relative to the case directory
        coordinateSystem cartesian;
        value        uniform ( 10 0 0 );
    }



SourceFiles
    mappedProfileFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedProfileFvPatchField_H
#define mappedProfileFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/mappedProfile/pointProfileDataList.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "memory/autoPtr/autoPtr.H"
#include "coordinate/systems/coordinateSystem.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mappedProfileFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedProfileFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- input file name
        fileName proFile_;

        //- input boundary region name
        string regionName_;

        //- input field name
        word fieldname_;

        //- convert to meters
        scalar convertToMeters_;

        //- alternative input negative surface dot product direction
        autoPtr<vector> negSurfDotDir_;

        //- alt field name
        autoPtr<word> altFieldName_;

        //- field after mapping on the patch
        Field<Type> mappedField_;

        //- coordinateSystem
        const word frameName_;
        mutable coordinateFrame* coorFramePtr_;

    // Private member functions

        //- updateValues
        virtual void updateValues();

        //- Based on fluent .prof data structure
        virtual void initialiseFluent(const dictionary&);

        //- Based on CSV data structure
        virtual void initialiseCSV(const dictionary&);

        virtual List<string> readLineElements
        (
            const string&,
            const char&,
            const Switch&
        );

        virtual void setTypeValue
        (
            DynamicList<Type>&,
            const List<label>&,
            const List<string>&
        );

        virtual void transformTypeValue
        (
            Type&,
            const vector&
        );


        virtual void calculateMap
        (
            labelList&,
            const Field<point>&,
            const Field<point>&
        );

        virtual tmp<Field<Type>> mapProfileToField
        (
            const pointField&,
            word, word
        );

protected:

    static List<fileName> inputFileNameList_;
    static PtrList<pointProfileDataList> inputData_;

public:

    //- Runtime type information
    TypeName("mappedProfile");


    // Constructors

        //- Construct from patch and internal field
        mappedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct by mapping given mappedProfileFvPatchField
        //  onto a new patch
        mappedProfileFvPatchField
        (
            const mappedProfileFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct from patch, internal field and dictionary
        mappedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct as copy
        mappedProfileFvPatchField
        (
            const mappedProfileFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedProfileFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedProfileFvPatchField
        (
            const mappedProfileFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedProfileFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access
            virtual const Field<Type>& mappedField() const
            {
                return mappedField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


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
