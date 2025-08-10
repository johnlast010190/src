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
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::wallContactTime

Group
    grpReportFunctionObjects

Description
    Write field value on selcted patches

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    wallContactTime.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallContactTime_H
#define functionObjects_wallContactTime_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class wallContactTime Declaration
\*---------------------------------------------------------------------------*/

class wallContactTime
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        // Read from dictonary

            //- input
            wordList patchNames_;
            labelList patchIDs_;
            wordList fieldNames_;

            //- result fields
            PtrList<volScalarField> phaseContactTimes_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        wallContactTime(const wallContactTime&);
        void operator=(const wallContactTime&);


public:

    //- Runtime type information
    TypeName("wallContactTime");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        wallContactTime
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~wallContactTime();


    // Member Functions
        void createFields
        (
            PtrList<volScalarField>& flds
        ) const;

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
