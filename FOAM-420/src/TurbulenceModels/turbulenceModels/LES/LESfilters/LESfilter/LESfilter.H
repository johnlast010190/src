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
    Foam::LESfilter

Description
    Abstract class for LES filters

SourceFiles
    LESfilter.C
    newFilter.C

\*---------------------------------------------------------------------------*/

#ifndef LESfilter_H
#define LESfilter_H

#include "fields/volFields/volFields.H"
#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class LESfilter Declaration
\*---------------------------------------------------------------------------*/

class LESfilter
{
    // Private data

        const fvMesh& mesh_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LESfilter(const LESfilter&);
        void operator=(const LESfilter&);


protected:

    //- Temporary function to ensure the coupled boundary conditions of the
    //  field are correct for filtering.
    //
    //  Following the rewrite of the turbulence models to use
    //  GeometricField::InternalField for sources etc. delta() will return a
    //  GeometricField::InternalField and filters will take a
    //  tmp<GeometricField::InternalField> argument and handle the coupled BCs
    //  appropriately
    template<class GeoFieldType>
    void correctBoundaryConditions(const tmp<GeoFieldType>& tgf) const
    {
        const_cast<GeoFieldType&>(tgf()).correctBoundaryConditions();
    }


public:

    //- Runtime type information
    TypeName("LESfilter");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESfilter,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& LESfilterDict
            ),
            (mesh, LESfilterDict)
        );


    // Constructors

        //- Construct from components
        LESfilter(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a reference to the selected LES filter
        static autoPtr<LESfilter> New
        (
            const fvMesh&,
            const dictionary&,
            const word& filterDictName="filter"
        );


    //- Destructor
    virtual ~LESfilter()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&) = 0;


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const = 0;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const = 0;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const = 0;

        virtual tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
