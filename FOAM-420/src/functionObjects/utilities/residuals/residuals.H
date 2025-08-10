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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2019 Esi Ltd.

Class
    Foam::functionObjects::residuals

Group
    grpUtilitiesFunctionObjects

Description
    Writes out the initial residual for specified fields.

Usage
    Example of function object specification:
    \verbatim
    residuals
    {
        type            residuals;
        libs            ("libutilityFunctionObjects.so");
        ...
        fields          (U p);
        writeFields     true; // default false;
    }
    \endverbatim

    Output data is written to the dir postProcessing/residuals/\<timeDir\>/
    For vector/tensor fields, e.g. U, where an equation is solved for each
    component, the largest residual of each component is written out.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    residuals.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_residuals_H
#define functionObjects_residuals_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class residuals Declaration
\*---------------------------------------------------------------------------*/

class residuals
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Fields to write residuals
        wordList fieldSet_;

        //- Flag to write the residual as a vol field
        bool writeFields_;

        //- Initialisation flag
        bool initialised_;


    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(Ostream& os) const;

        //- Create and store a residual field on the mesh database
        void createField(const word& fieldName);

        //- Write a residual field
        void writeField(const word& fieldName) const;

        //- Output file header information per primitive type value
        template<class Type>
        void writeFileHeader(Ostream& os, const word& fileName) const;

        //- Output file header information per primitive type value
        template<class Type>
        void writeFileBlockHeader(Ostream& os, const word& fileName) const;

        //- Initialise a residual field
        template<class Type>
        void initialiseField(const word& fieldName);

        //- Initialise a residual field
        template<class Type>
        void initialiseBlockField(const word& fieldName);

        //- Calculate the field min/max
        template<class Type>
        void writeResidual(const word& fieldName);

        template<class Type>
        void writeBlockResidual(const word& fieldName);


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        residuals(const residuals&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const residuals&) = delete;


public:

    //- Runtime type information
    TypeName("residuals");


    // Constructors

        //- Construct from Time and dictionary
        residuals
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~residuals();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the residuals
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "residuals/residualsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
