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
    Foam::functionObjects::subVolume

Group
    grpReportFunctionObjects

Description
    Subsetting of volume mesh
    Call subVolume.write() to calculate&write files.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    subVolume.C
    subVolumeTemplates.C
\*---------------------------------------------------------------------------*/

#ifndef functionObjects_subVolume_H
#define functionObjects_subVolume_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/Time/Time.H"
#include "fvMeshSubset/fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class subVolume Declaration
\*---------------------------------------------------------------------------*/

class subVolume
:
    public fvMeshFunctionObject
{

protected:

    // Protected data

        //- whether to write binary
        bool binary_;

        //- flag to write the processor addressing
        bool writeProcAddressingData_;

        // Read from dictonary

            //- input
            wordList fieldNames_;

            autoPtr<fvMeshSubset> subsetMesh_;

            autoPtr<Time> subMeshTime_;

            label writeInterval_;

        // Calculated

            //- Indices (into fieldNames_) of scalar/vector/tensor fields
            wordList scalarFields_;
            wordList vectorFields_;
            wordList sphericalTensorFields_;
            wordList symmTensorFields_;
            wordList tensorFields_;

         //add write control functionality

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        subVolume(const subVolume&);
        void operator=(const subVolume&);

        //- Find the fields in the list of the given type
        template<class T>
        void findFields(wordList& typeFieldNames);

        //- Sample and write field
        template<class T>
        void writeSubset(const word& fieldName);

        //- Sample and write all the fields of the given type
        template <class T>
        void writeSubset(const wordList& typeFields);


public:

   //- Runtime type information
    TypeName("subVolume");

    // Constructors

        //- Construct for given objectRegistry and dictionary
        subVolume
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~subVolume();


    // Member Functions

        //- Setup directory structure for subCase
        void setupSubCase(const word& subName);

        //- update lift/drag values
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the subVolume
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "subVolume/subVolumeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
