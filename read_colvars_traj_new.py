#!/usr/bin/env python3

class Colvars_traj:

    class NotTitleException(Exception):
        pass

    class Colvars_block:

        def __init__(self, keys):
            self._size = 0
            self._data = dict()
            for key in keys:
                self._data[key] = list()

        def __len__(self):
            if not self._data.keys():
                return 0
            else:
                for key in self._data.keys():
                    return len(self._data[key])

        def __getitem__(self, key):
            return self._data[key]

        def __contains__(self, key):
            return key in self._data.keys()

        def as_pandas(self, keys=None):
            import pandas as pd
            if keys is None:
                keys = self._data.keys()
            return pd.DataFrame({key: self._data[key] for key in keys})

    def __init__(self, filenames):
        self._data_block = list()
        self._frame = 0
        if type(filenames) == str:
            filenames = [filenames]
        if filenames:
            self._line_info = list()
            self._read_file(filenames)

    def _parse_title_line(self, line):
        new_keys = (line[1:]).split()
        if (new_keys[0] != 'step'):
            raise Colvars_traj.NotTitleException()
        result = {'step': [0]}
        # Find the boundaries of each column
        for i in range(1, len(new_keys)):
            if i == 1:
                pos = line.find(' '+new_keys[i]+' ')
                result['step'].append(pos)
            else:
                pos = line.find(' '+new_keys[i], result[new_keys[i-1]][0]+len(new_keys[i-1]))
                result[new_keys[i-1]].append(pos)
            result[new_keys[i]] = [result[new_keys[i-1]][1]]
        result[new_keys[-1]].append(len(line))
        return result

    def _parse_data_line(self, line, line_info, block):
        for key in line_info.keys():
            text = line[line_info[key][0]:line_info[key][1]].strip()
            if text[0] == '(':
                v_v = list(map(float, text[1:-1].split(',')))
            else:
                if key == 'step':
                    v_v = int(text)
                else:
                    v_v = float(text)
            block[key].append(v_v)

    def _read_file(self, filenames):
        for filename in filenames:
            with open(filename, 'r') as f_input:
                for line in f_input:
                    if (len(line) == 0):
                        continue
                    if (line[:1] == "@"):
                        continue
                    if (line[:1] == "#"):
                        try:
                            new_block_line_info = self._parse_title_line(line)
                            if self._line_info:
                                if self._line_info[-1] == new_block_line_info:
                                    continue
                            self._line_info.append(new_block_line_info)
                            self._data_block.append(Colvars_traj.Colvars_block(self._line_info[-1].keys()))
                        except Colvars_traj.NotTitleException:
                            pass
                        continue
                    else:
                        self._parse_data_line(line, self._line_info[-1], self._data_block[-1])
                        self._frame += 1

    def as_pandas(self, keys=None):
        import pandas as pd
        all_blocks = [block.as_pandas(keys) for block in self._data_block]
        return pd.concat(all_blocks, ignore_index=True)


if __name__ == '__main__':
    # test
    traj = Colvars_traj(['analyze-4.colvars.traj'])
    print(traj.as_pandas())
