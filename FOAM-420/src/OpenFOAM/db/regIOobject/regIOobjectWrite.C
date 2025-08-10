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

Description
    write function for regIOobjects

\*---------------------------------------------------------------------------*/

#include "db/regIOobject/regIOobject.H"
#include "db/Time/Time.H"
#include "include/OSspecific.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::regIOobject::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    if (!good())
    {
        SeriousErrorInFunction
            << "bad object " << name()
            << endl;

        return false;
    }

    if (instance().empty())
    {
        SeriousErrorInFunction
            << "instance undefined for object " << name()
            << endl;

        return false;
    }

#if defined(WIN64) || defined(WIN32)
    if (objectPath().size() > 256)
    {
        WarningInFunction << "The following file path exceeds maximum "
            << " allowable 256 characters for MS Windows : "
            << objectPath() << nl
            << endl;
    }
#endif

    //- uncomment this if you want to write global objects on master only
    //bool isGlobal = global();
    bool isGlobal = false;

    if (instance() == time().timeName())
    {
        // Mark as written to local directory
        isGlobal = false;
    }
    else if
    (
        instance() != time().system()
     && instance() != time().caseSystem()
     && instance() != time().constant()
     && instance() != time().caseConstant()
     && instance() != time().timeName()+"/uniform"
    )
    {
        // Update instance
        const_cast<regIOobject&>(*this).instance() = time().timeName();

        // Mark as written to local directory
        isGlobal = false;
    }

    if (OFstream::debug)
    {
        if (isGlobal)
        {
            Pout<< "regIOobject::write() : "
                << "writing (global) file " << objectPath();
        }
        else
        {
            Pout<< "regIOobject::write() : "
                << "writing (local) file " << objectPath();
        }
    }


    bool osGood = false;


    // Everyone check or just master
    bool masterOnly =
        isGlobal
     && (
            regIOobject::fileModificationChecking == timeStampMaster
         || regIOobject::fileModificationChecking == inotifyMaster
        );


    if (Pstream::master() || !masterOnly)
    {
        //if (mkDir(path()))
        //{
        //    // Try opening an OFstream for object
        //    OFstream os(objectPath(), fmt, ver, cmp);
        //
        //    // If any of these fail, return (leave error handling to Ostream
        //    // class)
        //    if (!os.good())
        //    {
        //        return false;
        //    }
        //
        //    if (!writeHeader(os))
        //    {
        //        return false;
        //    }
        //
        //    // Write the data to the Ostream
        //    if (!writeData(os))
        //    {
        //        return false;
        //    }
        //
        //    writeEndDivider(os);
        //
        //    osGood = os.good();
        //}
        osGood = fileHandler().writeObject(*this, fmt, ver, cmp, valid);
    }
    else
    {
        // Or scatter the master osGood?
        osGood = true;
    }

    if (OFstream::debug)
    {
        Pout<< " .... written" << endl;
    }

    // Only update the lastModified_ time if this object is re-readable,
    // i.e. lastModified_ is already set
    if (watchIndices_.size())
    {
        fileHandler().setUnmodified(watchIndices_.last());
    }

    return osGood;
}


bool Foam::regIOobject::write(const bool valid) const
{
    return writeObject
    (
        time().writeFormat(),
        IOstream::currentVersion,
        time().writeCompression(),
        valid
    );
}

// ************************************************************************* //
