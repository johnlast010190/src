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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::functionObjects::combineMultipointFields

Group
    grpUtilitiesFunctionObjects

Description
    Function object that combines fields and geometries from the different
    regions and combines them into region0. This function object is used
    together with the multipoint adjoint optimization to facilitate the
    post processing of the results.

Usage

    Example of function object specification:
    \verbatim
    combinedFields
    {
        type                combineMultipointFields;
        functionObjectLibs  ("libutilityFunctionObjects.so");
        writeControl        writeTime;

        fields               (p U);
        newFieldNames   (p_combined U_combined);
        operatingPoints
        (
            (geom1 p1)
            (geom2 p2)
        );
        totalFiName fi;
    };
    \endverbatim
SourceFiles
    combineMultipointFields.C

\*---------------------------------------------------------------------------*/

#ifndef combineMultipointFields_H
#define combineMultipointFields_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "primitives/Tuple2/Tuple2.H"
#include "fields/volFields/volFields.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
            Class combineMultipointFieldsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class combineMultipointFields
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        wordList fields_;

        wordList newFieldNames_;

        word totalFiName_;

        List<Tuple2<word, word>> operatingPoints_;

private:

        template<class Type>
        void calculateCombinedField
        (
            Type& field,
            const word& fieldName
        );

        void calculateAndWrite
        (
            const word& fieldName,
            const word& newFieldName
        );

        void combineAndWriteGeometry();

public:

    //- Runtime type information
    TypeName("combineMultipointFields")

    // Constructors

        combineMultipointFields
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

    virtual ~combineMultipointFields();

    // Member Functions

        //- Execute, currently does nothing
        bool execute() override;

        //- Sample and write
        bool write() override;

        //- Read the parameters
        bool read(const dictionary&) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

