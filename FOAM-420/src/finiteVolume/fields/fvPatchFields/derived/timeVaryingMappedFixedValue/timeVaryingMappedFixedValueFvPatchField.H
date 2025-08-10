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

Class
    Foam::timeVaryingMappedFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space and time.

    Supplied data should be specified in constant/boundaryData/\<patchname\>/
      - points             : pointField of locations
      - \<time\>/\<field\> : field of values at time \<time\>

    The default mode of operation (mapMethod planarInterpolation) is to project
    the points onto a plane (constructed from the first threee points) and
    construct a 2D triangulation and finds for the face centres the triangle it
    is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and triangulation
    and just use the value at the nearest vertex.

    Values are interpolated linearly between times.

Usage
    \table
        Property     | Description                | Required    | Default value
        setAverage   | Switch to activate setting of average value | no | false
        cylindricalCoords | Switch to convert values from cylindrical to cartesian coordinate system
        perturb      | Perturb points for regular geometries | no | 1e-5
        delabella      | Whether to use delabella triangulation | no | true
        fieldTableName | Alternative field name to sample | no| this field name
        mapMethod    | Type of mapping            | no | planarInterpolation
        offset       | Offset to mapped values    | no | Zero
        timeOffset   | Time offset to mapped values. The value is added to the execution time.| no | Zero
    \endtable

    \verbatim
    <patchName>
    {
        type            timeVaryingMappedFixedValue;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    timeVaryingMappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedFixedValueFvPatchField_H
#define timeVaryingMappedFixedValueFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "containers/Lists/FixedList/FixedList.H"
#include "db/Time/instant/instantList.H"
#include "triSurface/triSurfaceTools/pointToPointPlanarInterpolation.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class timeVaryingMappedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data
    using referenceFrameFvPatch<Type>::getFrameVelocity;

protected:

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        Switch setAverage_;

        //- If true, points are in cylindrical coordinates and must be converted to cartesian
        Switch cylindricalCoords_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Whether to use delabela triangulation method
        Switch delabella_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- Input stream
        autoPtr<IFstream> isPtr_;

        Switch f4stFormat_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Need to store sample points for use with frames
        autoPtr<pointField> samplePoints_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;

        //- Time offset value
        scalar timeOffset_;


public:

    //- Runtime type information
    TypeName("timeVaryingMappedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedFixedValueFvPatchField
        //  onto a new patch
        timeVaryingMappedFixedValueFvPatchField
        (
            const timeVaryingMappedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedFixedValueFvPatchField
        (
            const timeVaryingMappedFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMappedFixedValueFvPatchField
        (
            const timeVaryingMappedFixedValueFvPatchField<Type>&,
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
                new timeVaryingMappedFixedValueFvPatchField<Type>(*this, iF)
            );
        }

        //- Destructor
        virtual ~timeVaryingMappedFixedValueFvPatchField() = default;


    // Member functions

        // Access

            //- Return startSampledValues
            const Field<Type> startSampledValues()
            {
                 return startSampledValues_;
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

            //- Is defined in frame supported?
            virtual bool definedInFrameSupported() const
            {
                return true;
            }


        // Utility functions

            //- Transform points from local to global frame
            tmp<pointField> makePointsGlobal(pointField& samplePoints);

            //- Find boundary data inbetween current time and interpolate
            void checkTable();

            void checkTableF4st();

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
