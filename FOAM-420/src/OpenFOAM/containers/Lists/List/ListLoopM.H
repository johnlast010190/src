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
    (c) 2011-2015 OpenFOAM Foundation

Description
    Macros for accessing List elements

\*---------------------------------------------------------------------------*/

#ifndef ListLoop_H
#define ListLoop_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef vectorMachine

// Element access looping using [] for vector machines

#define List_FOR_ALL(f, i)                      \
        const label _n##i = (f).size();         \
        for (label i=0; i<_n##i; ++i)           \
        {

#define List_END_FOR_ALL  }

// Provide current element
#define List_CELEM(f, fp, i)  (fp[i])

// Provide current element
#define List_ELEM(f, fp, i)  (fp[i])

#define List_ACCESS(type, f, fp) \
    type* const __restrict__ fp = (f).begin()

#define List_CONST_ACCESS(type, f, fp) \
    const type* const __restrict__ fp = (f).begin()

#else

// Pointer looping for scalar machines

#define List_FOR_ALL(f, i)                      \
        label i = (f).size();                   \
        while (i--)                             \
        {                                       \

#define List_END_FOR_ALL  }

// Provide current element without incrementing pointer
#define List_CELEM(f, fp, i)  (*fp)

// Provide current element and increment pointer
#define List_ELEM(f, fp, i)  (*fp++)

#define List_ACCESS(type, f, fp) \
    type* __restrict__ fp = (f).begin()

#define List_CONST_ACCESS(type, f, fp) \
    const type* __restrict__ fp = (f).begin()

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
