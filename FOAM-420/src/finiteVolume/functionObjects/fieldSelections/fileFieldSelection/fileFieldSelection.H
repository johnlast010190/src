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
    (c) 2017-19 OpenCFD Ltd.

Class
    Foam::functionObjects::fileFieldSelection

Description
    Helper class to manage file-based field selections

SourceFiles
    fieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fileFieldSelection_H
#define functionObjects_fileFieldSelection_H

#include "functionObjects/fieldSelections/fieldSelection/fieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class fileFieldSelection Declaration
\*---------------------------------------------------------------------------*/

class fileFieldSelection
:
    public fieldSelection
{
private:

    // Private Member Functions

        //- No copy construct
        fileFieldSelection(const fileFieldSelection&) = delete;


protected:

    // Protected Member Functions

        //- Add registered GeometricField types to selection
        template<template<class> class PatchType, class MeshType>
        void addGeoFieldTypes(DynamicList<fieldInfo>& set) const;

        //- Add registered Internal types to selection
        void addInternalFieldTypes(DynamicList<fieldInfo>& set) const;

        //- Add registered uniform types to selection
        void addUniformFieldTypes(DynamicList<fieldInfo>& set) const;

        //- Add objects of a given type
        template<class Type>
        void addFromFile
        (
            const IOobjectList& allFileObjects,
            DynamicList<fieldInfo>& set
        ) const;


public:

    //- Construct from object registry
    fileFieldSelection
    (
        const objectRegistry& obr,
        const bool includeComponents = false
    );


    //- Destructor
    virtual ~fileFieldSelection() = default;


    // Member Functions

        //- Update the selection
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "functionObjects/fieldSelections/fileFieldSelection/fileFieldSelectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
