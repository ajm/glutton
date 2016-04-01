
def _sanity_check_data(data, columns) :
    for d in data :
        if len(d) != columns :
            return False

    return True

def _safe_len(datum) :
    if not datum :
        return 4
    return len(str(datum))

# data is assumed to be a list of tuples
def pretty_print_table(headings, data) :
    spacer = 4
    num_columns = len(headings)

    assert _sanity_check_data(data, num_columns)
    headings = [ h.capitalize() for h in headings ]

    # generate format string using column widths
    fmt = ""
    for i in range(num_columns) :
        fmt += ("%%-%ds" % (max([len(headings[i])] + [ _safe_len(d[i]) for d in data ]) + spacer))

    # print the actual table
    s = fmt % tuple(headings)
    print ""
    print s
    print '-'*len(s)

    for d in data :
        print fmt % tuple(d)

    print ""

if __name__ == '__main__' :
    from glutton.ensembl_downloader import EnsemblDownloader
    e = EnsemblDownloader()
    pretty_print_table(['species', 'release'], e.get_all_species('metazoa'))

