from collections import OrderedDict
import os
import re

import numpy as np


_types = {"float": np.float32, "double": np.float64,
          "int32": np.int32, "int64": np.int64,
          "uchar": np.ubyte}


def read_nrrd(file):
    with open(file, "rb") as fh:
        hdr = read_header(fh)
        data = read_data(fh, hdr)
    return data, hdr


def read_header(file):

    it = iter(file)
    magic_line = next(it)
    if hasattr(magic_line, "decode"):
        need_decode = True
        magic_line = magic_line.decode("ascii", "ignore")  # type: ignore[assignment]
        if not magic_line.startswith("NRRD"):  # type: ignore[arg-type]
            raise NotImplementedError

    header = OrderedDict()

    for line in it:
        if need_decode:
            line = line.decode("ascii", "ignore")  # type: ignore[assignment]

        line = line.rstrip()
        if line.startswith("#"):  # type: ignore[arg-type]
            continue
        elif line == "":
            break

        key, value = re.split(r"[:=?]", line, 1)  # type: ignore[type-var]
        key, value = key.strip(), value.strip()  # type: ignore[attr-defined]

        value = _get_value_type(key, value)

        header[key] = value

    return header


def _get_value_type(key, value):

    if key in ["dimension", "space dimension"]:
        return int(value)
    elif key in ["endian", "encoding"]:
        return value
    elif key in ["type"]:
        return _types[value]
    elif key in ["sizes"]:
        return [int(x) for x in value.split()]
    else:
        pass
        # raise NotImplementedError


def read_data(file, header):

    size = np.array(header["sizes"])
    dtype = header["type"]

    total_data_points = size.prod(dtype=np.int64)
    dtype_size = np.dtype(dtype).itemsize
    file.seek(-1 * dtype_size * total_data_points, os.SEEK_END)
    data = np.fromfile(file, dtype=dtype, sep="")
    data = data.reshape((size))
    return data


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    folderName = "RSMCRT/data/absorb/"
    fileName = "absorb.nrrd"
    file = folderName + fileName
    grid, hdr = read_nrrd(file)
        
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot()
    
    depths = np.linspace(-2.0, 2.0, int(hdr["sizes"][0]))
    ymid = int(hdr["sizes"][1]/2)
    zmid = int(hdr["sizes"][2]/2)
    data = grid
    
    fluence = np.mean(np.mean(data, axis =2), axis = 1)
    
    ax1.plot(depths, fluence)
    ax1.set_xlabel("Depth (cm)")
    ax1.set_ylabel("Fluence (-)")
    ax1.set_xlim([depths[-14],1.6])
    ax1.set_ylim([0, np.max(fluence)*1.1])
    
    #""" Used for validate 2
    c1 = 5.76
    k1 = 1.00
    c2 = 1.31
    k2 = 10.2
    delta = 0.047
    norm = 0.115
    #"""
    
    """ Used for validate 3
    c1 = 6.27
    k1 = 1.00
    c2 = 1.18
    k2 = 14.4
    delta = 0.261
    norm = 0.0151
    #"""
    
    fittingFunction = norm * (c1* np.exp((depths-1.95)*k1/delta) - c2*np.exp((depths-1.95)*k2/delta))
    ax1.plot(depths, fittingFunction)
       
    plt.show()

