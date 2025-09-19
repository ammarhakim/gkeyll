import urllib.request
import urllib.error
import os
import numpy

adas_data_dir = "."

ADAS_FILES = {
    "H":  {"acd": "acd12_h.dat",  "scd": "scd12_h.dat",  "plt": "plt12_h.dat",  "prb": "prb12_h.dat", "zmax": 1},
    "He": {"acd": "acd96_he.dat", "scd": "scd96_he.dat", "plt": "plt96_he.dat", "prb": "prb96_he.dat", "zmax": 2},
    "Li": {"acd": "acd96_li.dat", "scd": "scd96_li.dat", "plt": "plt96_li.dat", "prb": "prb96_li.dat", "zmax": 3},
    "Be": {"acd": "acd96_be.dat", "scd": "scd96_be.dat", "plt": "plt96_be.dat", "prb": "prb96_be.dat", "zmax": 4},
    "B":  {"acd": "acd89_b.dat",  "scd": "scd89_b.dat",  "plt": "plt89_b.dat",  "prb": "prb89_b.dat",  "zmax": 5},
    "C":  {"acd": "acd96_c.dat",  "scd": "scd96_c.dat",  "plt": "plt96_c.dat",  "prb": "prb96_n.dat",  "zmax": 6},
    "N":  {"acd": "acd96_n.dat",  "scd": "scd96_n.dat",  "plt": "plt96_n.dat",  "prb": "prb96_n.dat",  "zmax": 7},
    "O":  {"acd": "acd96_o.dat",  "scd": "scd96_o.dat",  "plt": "plt96_o.dat",  "prb": "prb96_o.dat",  "zmax": 8},
    "Ar": {"acd": "acd89_ar.dat", "scd": "scd89_ar.dat", "plt": "plt89_ar.dat", "prb": "prb89_ar.dat", "zmax": 18}
}

def fetch_file_urllib(filename, loc):
    base_url = "https://open.adas.ac.uk/download/adf11"
    filename_mod = filename.split("_")[0] + "/" + filename
    full_url = f"{base_url}/{filename_mod}"
    try:
        with urllib.request.urlopen(full_url) as response:
            if response.status != 200:
                raise urllib.error.HTTPError(full_url, response.status, "Failed to download", response.headers, None)
            content_bytes = response.read()
        if len(content_bytes) < 1000:
            raise ValueError(f'Could not fetch a valid file for {filename} from ADAS! Response was too short.')
        with open(loc, "wb") as f:
            f.write(content_bytes)
    except Exception as e:
        print(f"Error downloading {filename}: {e}")
        raise

def get_adas_file_loc(filename):
    if filename == "none":
        return
    loc = os.path.join(adas_data_dir, filename)
    if os.path.exists(loc):
        return loc
    fetch_file_urllib(filename, loc)
    return loc

class adas_adf11:
    def __init__(self, filename):
        self.filepath = os.path.join(adas_data_dir, filename)
        self.filename = filename
        self.file_type = self.filename[:3]
        try:
            self.imp = self.filename.split("_")[1].split(".")[0]
        except Exception:
            self.imp = None
        self.load()
        self.getClist()
    def getClist(self):
        self.clogNe = self.logNe.tolist()
        self.clogT = self.logT.tolist()
        self.clogdata = self.logdata.ravel().tolist()
    def load(self):
        loc = get_adas_file_loc(self.filename)
        with open(loc) as f:
            header = f.readline()
            self.n_ion, self.n_ne, self.n_T = numpy.int_(header.split()[:3])
            f.readline()
            line = f.readline()
            if all([a.isdigit() for a in line.split()]):
                self.metastables = numpy.int_(line.split())
                f.readline()
                line = f.readline()
            else:
                self.metastables = numpy.ones(self.n_ion + 1, dtype=int)
            logNe = []
            while len(logNe) < self.n_ne:
                logNe += [float(n) for n in line.split()]
                line = f.readline()
            logT = []
            while len(logT) < self.n_T:
                logT += [float(t) for t in line.split()]
                line = f.readline()
            subheader = line
            ldata, self.Z, self.MGRD, self.MPRT = [], [], [], []
            ind = 0
            while True:
                ind += 1
                try:
                    iprt, igrd, typ, z = subheader.split("/")[1:5]
                    self.Z.append(int(z.split("=")[1]))
                    self.MGRD.append(int(igrd.split("=")[1]))
                    self.MPRT.append(int(iprt.split("=")[1]))
                except Exception:
                    self.Z.append(ind + 1)
                    self.MGRD.append(1)
                    self.MPRT.append(1)
                drcofd = []
                while len(drcofd) < self.n_ne * self.n_T:
                    line = f.readline()
                    drcofd += [float(L) for L in line.split()]
                ldata.append(numpy.array(drcofd).reshape(self.n_T, self.n_ne))
                subheader = f.readline().replace("-", " ")
                if len(subheader) == 0 or subheader.isspace() or subheader[0] == "C":
                    break
        self.logNe = numpy.array(logNe)
        self.logT = numpy.array(logT)
        self.logdata = numpy.array(ldata)
        self.meta_ind = list(zip(self.Z, self.MGRD, self.MPRT))
    def data_for_charge_state(self, charge_state):
        return self.logdata[int(charge_state), :, :]

def write_files(adas_dict):
    for name, props in adas_dict.items():
        zmax = props["zmax"]
        ioniz_np, recomb_np = [], []
        for zi in range(zmax):
            ioniz = adas_adf11(props["scd"])
            ioniz_dat = ioniz.logdata[zi, :, :] - 6.0
            ioniz_np.append(ioniz_dat.flatten())
            recomb = adas_adf11(props["acd"])
            recomb_dat = recomb.logdata[zi, :, :] - 6.0
            recomb_np.append(recomb_dat.flatten())
        ioniz_np = numpy.array(ioniz_np)
        recomb_np = numpy.array(recomb_np)
        ioniz_np.tofile(f"ioniz_{name.lower()}.npy")
        recomb_np.tofile(f"recomb_{name.lower()}.npy")
        ioniz.logT.tofile(f"logT_{name.lower()}.npy")
        ioniz.logNe.tofile(f"logN_{name.lower()}.npy")

def main():
    # Download all ADAS files for supported elements
    for element, files in ADAS_FILES.items():
        print(f"Retrieving files for element {element}")
        for key in ["acd", "scd", "plt", "prb"]:
            get_adas_file_loc(files[key])
    # Convert and write numpy arrays for all elements in the dictionary
    write_files(ADAS_FILES)

if __name__ == "__main__":
    main()