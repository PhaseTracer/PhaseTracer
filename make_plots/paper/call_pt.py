"""
Run program and parse output
"""

import subprocess


def getter(lines):
    entries = {}
    for line in lines:
        if "=" in line:
            words = line.split("=")
            try:
                key = words[0].strip()
                entries[key] = float(words[1].strip())
            except:
                pass
    return entries

def cstr(x):
    if x is True:
        return "1"
    elif x is False:
        return "0"
    else:
        return str(x)

def transpose_list_dict(list_dict):
    return {k: [d.get(k) for d in list_dict] for k in list_dict[0].keys()}

def call_pt(program_name, *args):
    arg_list = [cstr(a) for a in args]
    r = subprocess.run(["../../bin/{}".format(program_name)] + arg_list,
                        stdout=subprocess.PIPE, encoding='utf8')

    lines = r.stdout.split("\n")
    return getter(lines)
