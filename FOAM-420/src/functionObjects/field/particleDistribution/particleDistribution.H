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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::functionObjects::particleDistribution

Group
    grpFieldFunctionObjects

Description
    Generates a particle distribution for lagrangian data at a given time.

Usage
    \verbatim
    particleDistribution1
    {
        type        particleDistribution;
        libs        ("libfieldFunctionObjects.so");
        ...
        cloud       "myCloud";
        nameVsBinWidth
        (
            (d 0.1)
            (U 10)
        );
        setFormat   raw;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property    | Description              |  Required  | Default value
        type        | Type name: particleDistribution | yes |
        cloud       | Name of cloud to process | Yes        |
        nameVsBinWidth | List of cloud field vs bin width | Yes |
        tagField    | Name of cloud field to use to group particles | no | none
        setFormat   | Output format            | yes         |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    particleDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_particleDistribution_H
#define functionObjects_particleDistribution_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/random/Random/Random.H"
#include "primitives/Tuple2/Tuple2.H"
#include "sampledSetWriters/writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class particleDistribution Declaration
\*---------------------------------------------------------------------------*/

class particleDistribution
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Cloud name
        word cloudName_;

        //- List of field name vs. bin width
        List<Tuple2<word, scalar>> nameVsBinWidth_;

        //- Tag field name - used to filter the particles into groups
        word tagFieldName_;

        //- Random number generator - used by distribution models
        Random rndGen_;

        //- Writer
        autoPtr<writer<scalar>> writerPtr_;


    // Protected Member Functions

        //- Generate the distribution
        void generateDistribution
        (
            const word& fieldName,
            const scalarField& field,
            const scalar binWidth,
            const label tag = -1
        );

        // Helper function to retrieve the particle data
        template<class Type>
        bool processField
        (
            const objectRegistry& obr,
            const label fieldi,
            const List<DynamicList<label>>& addr
        );

        //- Disallow default bitwise copy construct
        particleDistribution(const particleDistribution&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const particleDistribution&) = delete;


public:

    //- Runtime type information
    TypeName("particleDistribution");


    // Constructors

        //- Construct from Time and dictionary
        particleDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~particleDistribution();


    // Member Functions

        //- Read the particleDistribution data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the particleDistribution
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleDistribution/particleDistributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
