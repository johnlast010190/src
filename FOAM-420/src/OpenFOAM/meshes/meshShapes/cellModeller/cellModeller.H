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
    Foam::cellModeller

Description
    A static collection of cell models, and a means of looking them up.

SourceFiles
    cellModeller.C
    cellModellerIO.C
    globalCellModeller.C

\*---------------------------------------------------------------------------*/

#ifndef cellModeller_H
#define cellModeller_H

#include "meshes/meshShapes/cellModel/cellModel.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "containers/HashTables/HashTable/HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellModeller Declaration
\*---------------------------------------------------------------------------*/

class cellModeller
{
    // Private data

        //- PtrList of models
        static PtrList<cellModel> models_;

        //- List of model pointers
        static List<cellModel*> modelPtrs_;

        //- HashTable of model pointers
        static HashTable<const cellModel*> modelDictionary_;


public:

    // Constructors

        //- Construct from central "cellModels" file
        cellModeller();

    //- Destructor
    ~cellModeller();


    // Member functions

        //- Look up a model by name and return a pointer to the model or nullptr
        static const cellModel* lookup(const word&);

        //- Look up a model by index and return a pointer to the model or
        //  nullptr
        static const cellModel* lookup(const label i)
        {
            return modelPtrs_[i];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
