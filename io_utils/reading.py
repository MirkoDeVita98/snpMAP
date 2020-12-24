from constants import META_DATA, HEADER, fromatdesc2code


REF_COLUMN = 3
ALT_COLUMN = 4
INFO_COLUMN = 7


class VCFReader:


    def __init__(self, fpath):
        self.fpath = fpath
        self.file = None
        self.columns = None
        self.columns_full = None
        self.formats = []
        self.infos = []
        self.flag_infos = set()


    def get_record(self):
        with open(self.fpath) as fr:
            for line in fr:
                # print(line[:5])
                if line.startswith(META_DATA):
                    if line[2:2 + len('FORMAT')] == 'FORMAT':
                        description_start = \
                            line.find('Description=') + len('Description=')
                        tmp1 = line[description_start:]
                        p1 = tmp1.find('"') + 1
                        tmp2 = tmp1[p1:]
                        p2 = tmp2.find('"')
                        format_description = tmp1[p1: p1 + p2]
                        self.formats.append(
                            fromatdesc2code[format_description])
                    if line[2:2 + len('INFO')] == 'INFO':
                        id_start = line.find('ID=') + len('ID=')
                        tmp1 = line[id_start:]
                        id_end = tmp1.find(',')
                        self.infos.append(tmp1[:id_end])
                        number_start = tmp1.find('Number=') + len('Number=')
                        tmp2 = tmp1[number_start:]
                        number_end = tmp2.find(',')
                        if tmp2[:number_end] == '0':
                            self.flag_infos.add(tmp1[:id_end])
                elif line.startswith(HEADER):
                    self.columns_full = line.split()
                    self.columns_full[0] = self.columns_full[0][1:]
                    if 'FORMAT' in line:
                        self.columns = \
                            line[:line.find('FORMAT') + len('FORMAT')].split()
                        self.columns[0] = self.columns[0][1:]
                    else:
                        self.columns = self.columns_full
                else:
                    if len(self.formats) > 0:
                        fmt_found = None
                        border = -1
                        for fmt in self.formats:
                            b = line.find(fmt)
                            if b > border:
                                border = b + len(fmt)
                                fmt_found = fmt
                        # genotype_data = line[border:]
                        line = line[:border]

                    vals = line.split()

                    # print(f'Line = {line}')
                    # print(f'Vals = {vals}')

                    refs = vals[REF_COLUMN].split(',')
                    vals[REF_COLUMN] = refs
                    alts = vals[ALT_COLUMN].split(',')
                    vals[ALT_COLUMN] = alts

                    info = vals[INFO_COLUMN]
                    # print(f'Info = {info}')
                    info_fields = info.split(';')
                    # print(f'Info fields = {info_fields}')
                    info_fields_split = [(el.split('=')) 
                        for el in info_fields if el not in self.flag_infos]
                    # print(f'Info fields split = {info_fields_split}')
                    info_dict = {el[0]: el[1] for el in info_fields_split}
                    info_dict.update({el: True 
                        for el in self.flag_infos if el in info_fields})

                    vals[INFO_COLUMN] = info_dict

                    yield {el[0]: el[1] for el in zip(self.columns, vals)}

                    # record = {el[0]: el[1] for el in zip(self.columns, vals)}
                    # record['line'] = line

                    # yield record