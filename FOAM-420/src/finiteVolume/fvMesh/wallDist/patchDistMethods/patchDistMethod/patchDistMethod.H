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
    (c) 2017 OpenCFD Ltd.
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::patchDistMethod

Description
    Specialisation of patchDist for wall distance calculation

SourceFiles
    patchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistMethod_H
#define patchDistMethod_H

#include "db/dictionary/dictionary.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fields/volFields/volFieldsFwd.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class patchDistMethod Declaration
\*---------------------------------------------------------------------------*/

class patchDistMethod
{

protected:

    // Protected Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Set of patch IDs
        const labelHashSet patchIDs_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchDistMethod(const patchDistMethod&);

        //- Disallow default bitwise assignment
        void operator=(const patchDistMethod&);

public:

    //- Runtime type information
    TypeName("patchDistMethod");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            patchDistMethod,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const labelHashSet& patchIDs
            ),
            (dict, mesh, patchIDs)
        );


    // Constructors

        //- Construct from mesh and patch ID set
        patchDistMethod
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );


    // Selectors

        static autoPtr<patchDistMethod> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        static autoPtr<patchDistMethod> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs,
            const word& defaultPatchDistMethod
        );


    //- Destructor
    virtual ~patchDistMethod();


    // Static Functions

        //- Return the patch types for y and n
        //  These are fixedValue for the set provided otherwise zero-gradient
        template<class Type>
        static wordList patchTypes
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );


    // Member Functions

        //- Return the patchIDs
        const labelHashSet& patchIDs() const
        {
            return patchIDs_;
        }

        //- Update cached geometry when the mesh moves
        virtual bool movePoints()
        {
            return true;
        }

        //- Update cached topology and geometry when the mesh changes
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Correct the given distance-to-patch field
        virtual bool correct(volScalarField& y) = 0;

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(volScalarField& y, volVectorField& n) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/wallDist/patchDistMethods/patchDistMethod/patchDistMethodTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
