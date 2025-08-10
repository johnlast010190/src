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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"
#include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::initMatrixInterfaces
(
    const bool add,
    const FieldField<Field, LUType>& interfaceCoeffs,
    const Field<Type>& psiif,
    Field<Type>& result
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces_, interfacei)
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    psiif,
                    interfaceCoeffs[interfacei],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfacei]),
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces_.size();
            interfacei++
        )
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    psiif,
                    interfaceCoeffs[interfacei],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfacei]),
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, LUType>& interfaceCoeffs,
    const Field<Type>& psiif,
    Field<Type>& result
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        // Block until all sends/receives have been finished
        if (Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking)
        {
            IPstream::waitRequests();
            OPstream::waitRequests();
        }

        forAll(interfaces_, interfacei)
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    psiif,
                    interfaceCoeffs[interfacei],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfacei]),
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll(patchSchedule, i)
        {
            label interfacei = patchSchedule[i].patch;

            if (interfaces_.set(interfacei))
            {
                if (patchSchedule[i].init)
                {
                    interfaces_[interfacei].initInterfaceMatrixUpdate
                    (
                        result,
                        add,
                        psiif,
                        interfaceCoeffs[interfacei],
                      //Amultiplier<Type, LUType>(interfaceCoeffs[interfacei]),
                        Pstream::commsTypes::scheduled
                    );
                }
                else
                {
                    interfaces_[interfacei].updateInterfaceMatrix
                    (
                        result,
                        add,
                        psiif,
                        interfaceCoeffs[interfacei],
                      //Amultiplier<Type, LUType>(interfaceCoeffs[interfacei]),
                        Pstream::commsTypes::scheduled
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces_.size();
            interfacei++
        )
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    psiif,
                    interfaceCoeffs[interfacei],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfacei]),
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
