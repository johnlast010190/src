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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::motionSolver

Description
    Virtual base class for mesh motion solver.

SourceFiles
    motionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef motionSolver_H
#define motionSolver_H

#include "meshes/primitiveShapes/point/pointField.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class motionSolver Declaration
\*---------------------------------------------------------------------------*/

class motionSolver
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Is the motion solid body motion?
        bool isSolidBody_;

        //- The motion control function
        PtrList<solidBodyMotionFunction> SBMFs_;

        //- The reference frames with motions
        UPtrList<coordinateFrame> refFrames_;


public:

    //- Runtime type information
    TypeName("motionSolver");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionSolver,
            dictionary,
            (const polyMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Selectors

        //- Select constructed from polyMesh and dictionary
        static autoPtr<motionSolver> New(const polyMesh&, const dictionary&);

        //- Class used for the construction of PtrLists of motionSolvers
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh);

            autoPtr<motionSolver> operator()(Istream& is) const;
        };


    // Constructors

        //- Construct from polyMesh, dictionary and type
        motionSolver
        (
            const polyMesh& mesh,
            const dictionary&,
            const word& type
        );

        //- Clone function
        virtual autoPtr<motionSolver> clone() const;


    //- Destructor
    virtual ~motionSolver();


    // Member Functions

        //- Return reference to mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Const access to the coefficients dictionary
        const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return transformation septernion
        void initTransforms();

        //- Is the motion incremental?
        bool isIncrementalMotion() const;

        //- Return transformation septernion
        septernion transformation
        (
            const label transformationI = 0,
            const wordList& nestedFrames = wordList::null()
        ) const;

        //- Provide new points for motion. Solve for motion.
        virtual tmp<pointField> newPoints();

        //- Provide current points for motion. Use current motion field.
        virtual tmp<pointField> curPoints() const = 0;

        virtual void twoDCorrectPoints(pointField&) const;

        //- Solve for motion
        virtual void solve() = 0;

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&) = 0;

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&) = 0;

        //- Optionally write motion state information for restart
        virtual bool write() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
