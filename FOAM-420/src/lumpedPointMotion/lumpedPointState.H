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
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::lumpedPointState

Description
    The 'state' of lumped points corresponds to positions and rotations.

    Encapsulates the response from the external application.
    Entry place for applying relaxation and sub-stepping etc.

SourceFiles
    lumpedPointState.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointState_H
#define lumpedPointState_H

#include "db/dictionary/dictionary.H"
#include "primitives/Scalar/lists/scalarList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "fields/Fields/tensorField/tensorField.H"
#include "primitives/enums/Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                      Class lumpedPointState Declaration
\*---------------------------------------------------------------------------*/

//- Bundling of position/rotation
class lumpedPointState
{
public:

    //- Input format types
    enum class inputFormatType
    {
        PLAIN,
        DICTIONARY
    };

    // Static data

    //- Names for the input format types
    static const Enum<inputFormatType> formatNames;


private:

    // Private data

        //- Positions of lumped points
        pointField  points_;

        //- Orientation of lumped points (as Euler angles)
        vectorField angles_;

        //- Euler angles in degrees instead radians
        bool degrees_;

        //- Tensor rotation of lumped points
        mutable tensorField* rotationPtr_;


    // Private Member Functions

        void calcRotations() const;

        void readDict(const dictionary& dict);

public:

    // Constructors

        //- Construct null
        lumpedPointState();

        //- Copy constructor
        lumpedPointState(const lumpedPointState& rhs);

        //- Construct from points with zero-rotation
        lumpedPointState(const pointField& pts);

        //- Construct from points with zero-rotation
        lumpedPointState(tmp<pointField>& pts);

        //- Construct from dictionary
        lumpedPointState(const dictionary& dict);


    //- Destructor
    virtual ~lumpedPointState();


    // Member Functions

        //- Has positions and consistent number of rotations?
        inline bool valid() const;

        //- If no points were specified
        inline bool empty() const;

        //- The number of points
        inline label size() const;

        //- The points corresponding to mass centres
        inline const pointField& points() const;

        //- The orientation of the points (mass centres)
        inline const vectorField& angles() const;

        //- The local-to-global transformation for each point
        inline const tensorField& rotations() const;

        //- Relax the state
        //  alpha = 1 : no relaxation
        //  alpha < 1 : relaxation
        //  alpha = 0 : do nothing
        void relax(const scalar alpha, const lumpedPointState& prev);


        //- Read input as dictionary content
        bool readData(Istream& is);

        //- Output as dictionary content
        bool writeData(Ostream& os) const;

        //- Output as dictionary content
        void writeDict(Ostream& os) const;

        //- Read input as plain content
        bool readPlain(Istream& is);

        //- Output as plain content
        void writePlain(Ostream& os) const;

        //- Read input from file (master process only) using specified format
        bool readData(const inputFormatType& fmt, const fileName& file);

        //- Output as VTK file for debugging/visualization
        //  The points are joined as lines, the rotation is visualized
        //  by planes, write as VTK PolyData format.
        void writeVTP(const fileName& outputFile, const vector& axis) const;


    // Member Operators

        //- Assignment operator
        void operator=(const lumpedPointState& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lumpedPointStateI.H"

#endif

// ************************************************************************* //
