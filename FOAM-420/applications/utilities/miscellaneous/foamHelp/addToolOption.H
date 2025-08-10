argList::validArgs.append("tool");
wordList opts;
for (const auto& p : helpType::dictionaryConstructorTable_()) {
    opts.append(p.first);
}

string note = "Valid <tool> options include:";
forAll(opts, i)
{
    note = note + ' ' + opts[i];
}

argList::notes.append(note);
