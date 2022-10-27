import inspect

def PrintData(Model):
    keyword_list = list()
    label_list = list()
    DATA_PREFIX = 'Data_'

    for v in sorted(Model.__dict__.keys()):
        if v.startswith(DATA_PREFIX):
            keyword_list.append(v)

    tv = DATA_PREFIX + 'Time'
    if tv in keyword_list:
        i = keyword_list.index(tv)
        keyword_list = [tv] + keyword_list[:i] + keyword_list[i+1:]

    label_list = [v[len(DATA_PREFIX):] for v in keyword_list]

    value_list = [Model.__dict__[k] for k in keyword_list]

    # print header line
    print('\t'.join(label_list))

    # print data points
    for v in zip(*value_list):
        print('\t'.join([str(i) for i in v]))