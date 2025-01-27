from collections import OrderedDict
import os
import re

import numpy as np
class read_nrrd_class:
    def __init__(self):
        self._types ={"float": np.float32, "double": np.float64,
            "int32": np.int32, "int64": np.int64,
            "uchar": np.ubyte}

    def read_nrrd(self, file):
        with open(file, "rb") as fh:
            hdr = self.read_header(fh)
            data = self.read_data(fh, hdr)
        return data, hdr

    def read_header(self, file):

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

            splitList = re.split(r"[:=?]", line)  # type: ignore[type-var]
            key = splitList[0]
            value = splitList[-1]
            key, value = key.strip(), value.strip()  # type: ignore[attr-defined]

            value = self._get_value_type(key, value)

            header[key] = value

        return header


    def _get_value_type(self, key, value):

        if key in ["dimension", "space dimension", "num_spheres", "numOptProp", "nphotons"]:
            return int(value)
        elif key in ["endian", "encoding", "annulus_type", "focus_type", "units", "grid_data", "real_size", "source", "experiment", "dector", "symmetryType"]:
            return value
        elif key in ["type"]:
            return self._types[value]
        elif key in ["sizes"]:
            return [int(x) for x in value.split()]
        elif key in ["radius", "focalLength", "rhi", "rlo", "sigma", "beam_size", "tau", 
                     "musb", "muab", "musc", "muac", "hgga", '"mua%   1"', '"mus%   1"', '"mur%   1"',
                     '"hgg%   1"', '"n%   1"', '"mua%   2"', '"mus%   2"', '"mur%   2"', '"hgg%   2"', 
                     '"n%   2"', '"mua%   3"', '"mus%   3"', '"mur%   3"', '"hgg%   3"', '"n%   3"', 
                     '"position%   1"', '"position%   2"', '"position%   3"', '"boundinglength%   1"',
                     '"boundinglength%   2"', '"boundinglength%   3"', '"BoxDimensions%   1"', 
                     '"BoxDimensions%   2"', '"BoxDimensions%   3"']:
            return float(value)
        else:
            print(f"Error not implemented {key}")
            print(repr(key))
            pass
            # raise NotImplementedError


    def read_data(self, file, header):

        size = np.array(header["sizes"])
        dtype = header["type"]

        total_data_points = size.prod(dtype=np.int64)
        dtype_size = np.dtype(dtype).itemsize
        file.seek(-1 * dtype_size * total_data_points, os.SEEK_END)
        data = np.fromfile(file, dtype=dtype, sep="")
        data = data.reshape((size))
        return data
    
