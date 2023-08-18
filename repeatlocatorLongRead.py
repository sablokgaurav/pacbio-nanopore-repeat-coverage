def repeatlocatorLongRead(long_read_file, polyATGCstretch_type = None):
    """summary_line
    a long read repeat coverage calculator,given an long read file
    before assembly either direct from the sequencing runs or after 
    the cleaning, it calculates the total amount of the repeat stretches
    present in the sequencing reads and you can plot them before assembly
    to see how much repetitive your sequencing is and how to set the 
    long read assembly parameters.
    Keyword arguments:
    argument -- description
    infile = long read sequencing file
    Return: return_description
    a dataframe, with the sequences, and position of the repeats
    and also the coveage for the same.
    """
    
    if polyATGCstretch_type == "A":
        read_long_reads = list(filter(None,[x.strip() for x in open(long_read_file).readlines()]))
        long_read_conversion = {}
        for i in read_long_reads:
            if i.startswith(">"):
                path = i.strip()
                if i not in long_read_conversion:
                    long_read_conversion[i] = ""
                    continue
            long_read_conversion[path] += i.strip()
        ids = list(long_read_conversion.keys())
        sequences = list(long_read_conversion.values())
        long_read_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequences"})
        long_read_dataframe["repeat_locator"] = long_read_dataframe["sequences"]. \
                         apply(lambda n: str(n)).apply(lambda n: re.findall(r"[aA]{3,}", n))
        repeats =  list(set([j for i in (long_read_dataframe["repeat_locator"].to_list()) for j in i]))
        long_read_dataframe["fraction_length"] = long_read_dataframe["sequences"].apply(lambda n: \
                                       [(n.find(repeats[i]),(n.find(repeats[i])+len(repeats[i]))) \
                                                                         for i in range(len(repeats))])
        fractions = long_read_dataframe["fraction_length"].to_list()
        sequences = long_read_dataframe["sequences"].to_list()
        sequences_length = [len(i) for i in sequences]
        fraction_length_plot = list(map(lambda n: sum([i[1]-i[0] for i in n]),fractions))
        long_read_dataframe["fraction_length_coverage"] = pd.DataFrame([sequences_length[i]/fraction_length_plot[i] \
                                                            for i in range(len(sequences_length))]). \
                                                                        rename(columns = {0: "Repeat_coverage"})
        return long_read_dataframe
    if polyATGCstretch_type == "T":
        read_long_reads = list(filter(None,[x.strip() for x in open(long_read_file).readlines()]))
        long_read_conversion = {}
        for i in read_long_reads:
            if i.startswith(">"):
                path = i.strip()
                if i not in long_read_conversion:
                    long_read_conversion[i] = ""
                    continue
            long_read_conversion[path] += i.strip()
        ids = list(long_read_conversion.keys())
        sequences = list(long_read_conversion.values())
        long_read_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequences"})
        long_read_dataframe["repeat_locator"] = long_read_dataframe["sequences"]. \
                         apply(lambda n: str(n)).apply(lambda n: re.findall(r"[tT]{3,}", n))
        repeats =  list(set([j for i in (long_read_dataframe["repeat_locator"].to_list()) for j in i]))
        long_read_dataframe["fraction_length"] = long_read_dataframe["sequences"].apply(lambda n: \
                                       [(n.find(repeats[i]),(n.find(repeats[i])+len(repeats[i]))) \
                                                                         for i in range(len(repeats))])
        fractions = long_read_dataframe["fraction_length"].to_list()
        sequences = long_read_dataframe["sequences"].to_list()
        sequences_length = [len(i) for i in sequences]
        fraction_length_plot = list(map(lambda n: sum([i[1]-i[0] for i in n]),fractions))
        long_read_dataframe["fraction_length_coverage"] = pd.DataFrame([sequences_length[i]/fraction_length_plot[i] \
                                                            for i in range(len(sequences_length))]). \
                                                                        rename(columns = {0: "Repeat_coverage"})
        return long_read_dataframe
    if polyATGCstretch_type == "C":
        read_long_reads = list(filter(None,[x.strip() for x in open(long_read_file).readlines()]))
        long_read_conversion = {}
        for i in read_long_reads:
            if i.startswith(">"):
                path = i.strip()
                if i not in long_read_conversion:
                    long_read_conversion[i] = ""
                    continue
            long_read_conversion[path] += i.strip()
        ids = list(long_read_conversion.keys())
        sequences = list(long_read_conversion.values())
        long_read_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequences"})
        long_read_dataframe["repeat_locator"] = long_read_dataframe["sequences"]. \
                         apply(lambda n: str(n)).apply(lambda n: re.findall(r"[cC]{3,}", n))
        repeats =  list(set([j for i in (long_read_dataframe["repeat_locator"].to_list()) for j in i]))
        long_read_dataframe["fraction_length"] = long_read_dataframe["sequences"].apply(lambda n: \
                                       [(n.find(repeats[i]),(n.find(repeats[i])+len(repeats[i]))) \
                                                                         for i in range(len(repeats))])
        fractions = long_read_dataframe["fraction_length"].to_list()
        sequences = long_read_dataframe["sequences"].to_list()
        sequences_length = [len(i) for i in sequences]
        fraction_length_plot = list(map(lambda n: sum([i[1]-i[0] for i in n]),fractions))
        long_read_dataframe["fraction_length_coverage"] = pd.DataFrame([sequences_length[i]/fraction_length_plot[i] \
                                                            for i in range(len(sequences_length))]). \
                                                                        rename(columns = {0: "Repeat_coverage"})
        return long_read_dataframe
    if polyATGCstretch_type == "G":
        read_long_reads = list(filter(None,[x.strip() for x in open(long_read_file).readlines()]))
        long_read_conversion = {}
        for i in read_long_reads:
            if i.startswith(">"):
                path = i.strip()
                if i not in long_read_conversion:
                    long_read_conversion[i] = ""
                    continue
            long_read_conversion[path] += i.strip()
        ids = list(long_read_conversion.keys())
        sequences = list(long_read_conversion.values())
        long_read_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequences"})
        long_read_dataframe["repeat_locator"] = long_read_dataframe["sequences"]. \
                         apply(lambda n: str(n)).apply(lambda n: re.findall(r"[gG]{3,}", n))
        repeats =  list(set([j for i in (long_read_dataframe["repeat_locator"].to_list()) for j in i]))
        long_read_dataframe["fraction_length"] = long_read_dataframe["sequences"].apply(lambda n: \
                                       [(n.find(repeats[i]),(n.find(repeats[i])+len(repeats[i]))) \
                                                                         for i in range(len(repeats))])
        fractions = long_read_dataframe["fraction_length"].to_list()
        sequences = long_read_dataframe["sequences"].to_list()
        sequences_length = [len(i) for i in sequences]
        fraction_length_plot = list(map(lambda n: sum([i[1]-i[0] for i in n]),fractions))
        long_read_dataframe["fraction_length_coverage"] = pd.DataFrame([sequences_length[i]/fraction_length_plot[i] \
                                                            for i in range(len(sequences_length))]). \
                                                                        rename(columns = {0: "Repeat_coverage"})
        return long_read_dataframe
