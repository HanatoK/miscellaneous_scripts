#!/usr/bin/env python3

# i/o constants
it_width = 12
cv_prec = 14
cv_width = 21
en_prec = 14
en_width = 21


def wrap_string(s, nchars):
    if len(s) == 0:
        return ' ' * nchars
    else:
        if len(s) <= nchars:
            return s + ' ' * (nchars - len(s))
        else:
            return s[0:nchars]


def write_colvars_traj(df, outputfile):
    columns = df.columns
    missing_step = False
    if 'step' not in columns:
        missing_step = True
    first_line = True
    with open(outputfile, 'w') as f_output:
        for i in range(0, len(df)):
            if first_line:
                # write labels
                f_output.write(f'# {"step":{it_width-2}} ')
                if not missing_step:
                    f_output.write(' ')
                    for label in columns:
                        if label != 'step':
                            if isinstance(df.iloc[i][label], float) or \
                               isinstance(df.iloc[i][label], int):
                                f_output.write(f' {wrap_string(label, cv_width)}')
                            else:
                                raise RuntimeError('Currently does not support non-numeric fields.')
                f_output.write('\n')
                first_line = False
            if missing_step:
                step = i
            else:
                step = df.iloc[i]['step']
            f_output.write(f'{int(step):{it_width}d} ')
            f_output.write(' ')
            for label in columns:
                if label != 'step':
                    f_output.write(f' {df.iloc[i][label]:{cv_width}.{cv_prec}e}')
            f_output.write('\n')


if __name__ == '__main__':
    from plot_colvars_traj import Colvars_traj
    traj = Colvars_traj(['experiment.traj'])
    write_colvars_traj(traj.as_pandas(), 'xxx')
    traj2 = Colvars_traj(['xxx'])
    write_colvars_traj(traj2.as_pandas(), 'xxx2')
