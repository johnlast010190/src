#include "UPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::UPstream::commsTypes,
        3
    >::names[] =
    {
        "blocking",
        "scheduled",
        "nonBlocking"
    };
}

const Foam::NamedEnum<Foam::UPstream::commsTypes, 3>
    Foam::UPstream::commsTypeNames;


bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

Foam::LIFOStack<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<int> Foam::UPstream::myProcNo_(10);

Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(10);

Foam::DynamicList<Foam::label> Foam::UPstream::parentCommunicator_(10);

int Foam::UPstream::msgType_(1);


Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::linearCommunication_(10);

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::treeCommunication_(10);


// Allocate a serial communicator. This gets overwritten in parallel mode
// (by UPstream::setParRun())
Foam::UPstream::communicator serialComm
(
    -1,
    Foam::labelList(1, Foam::Zero),
    false
);


bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 16)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.read(Foam::debug::optimisationSwitches().lookup("commsType"))
);
