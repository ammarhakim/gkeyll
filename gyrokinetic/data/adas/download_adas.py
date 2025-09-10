import requests
import os
import shutil
import sys

adas_data_dir = "."

#if arg =atomic element symbol, attempt to download that set
#if arg =adf11 filename, attempt to download that file
#default arg: download default set of files
def main(arg="all"):
    files = adas_files_dict()
    if arg == "all":
        for x in files.keys():
            print(f"Retrieving files for element {x}")
            download_element(x,files)
    elif len(arg)<=2:
        if arg in files:
            download_element(arg,files)
        else:
            print(f"Atomic symbol {arg} not found")
            print("Supported elements:")
            for x in files.keys():
                print(x)
    else:
        print(f"Attempting to locate file {arg}")
        print(f"Will download if not in {adas_data_dir}")
        filename = get_adas_file_loc(arg)
        
#Return full path for specific file and download if necessary
def get_adas_file_loc(filename):
    #if adas_data_dir doesn't exist - create it
    
    if filename == "none":
        #Don't load a file
        return
    elif os.path.exists(adas_data_dir + os.sep + filename):
        #File is already on system
        return adas_data_dir + os.sep + filename
    elif os.path.exists(filename):
        #filename is full path
        return filename
    else:
        loc = adas_data_dir + os.sep + filename
        fetch_file(filename,loc)
        return loc

    
#Fetch file from open.adas.ac.uk    
def fetch_file(filename,loc):
    url = "https://open.adas.ac.uk/download/adf11/"
    
    filename_mod = filename.split("_")[0] + "/" + filename
    
    r = requests.get(url + "/" + filename_mod)
    
    if(len(r.text)) < 1000:
        raise ValueError(f'Could not fetch {filename} from ADAS!')
    
    with open(loc, "wb") as f:
        f.write(r.content)

def download_element(symbol,files):
    files_of_element = [files[symbol][k] for k in files[symbol]]
    for x in files_of_element:
        get_adas_file_loc(x)#only download new files
        
def adas_files_dict():
    files = {}
    files["H"] = {}
    files["He"] = {}
    files["Li"] = {}
    files["Be"] = {}
    files["B"] = {}
    files["C"] = {}
    files["N"] = {}
    files["O"] = {}
    files["Ar"] = {}
    files["H"]["acd"] = "acd12_h.dat"
    files["H"]["scd"] = "scd12_h.dat"
    files["H"]["plt"] = "plt12_h.dat"
    files["H"]["prb"] = "prb12_h.dat"
    files["He"]["acd"] = "acd96_he.dat"
    files["He"]["scd"] = "scd96_he.dat"
    files["He"]["plt"] = "plt96_he.dat"
    files["He"]["prb"] = "prb96_he.dat"
    files["Li"]["acd"] = "acd96_li.dat"
    files["Li"]["scd"] = "scd96_li.dat"
    files["Li"]["plt"] = "plt96_li.dat"
    files["Li"]["prb"] = "prb96_li.dat"
    files["Be"]["acd"] = "acd96_be.dat"
    files["Be"]["scd"] = "scd96_be.dat"
    files["Be"]["plt"] = "plt96_be.dat"
    files["Be"]["prb"] = "prb96_be.dat"
    files["B"]["acd"] = "acd89_b.dat"
    files["B"]["scd"] = "scd89_b.dat"
    files["B"]["plt"] = "plt89_b.dat"
    files["B"]["prb"] = "prb89_b.dat"
    files["C"]["acd"] = "acd96_c.dat"
    files["C"]["scd"] = "scd96_c.dat"
    files["C"]["plt"] = "plt96_c.dat"
    files["C"]["prb"] = "prb96_n.dat"
    files["N"]["acd"] = "acd96_n.dat"
    files["N"]["scd"] = "scd96_n.dat"
    files["N"]["plt"] = "plt96_n.dat"
    files["N"]["prb"] = "prb96_n.dat"
    files["O"]["acd"] = "acd96_o.dat"
    files["O"]["scd"] = "scd96_o.dat"
    files["O"]["plt"] = "plt96_o.dat"
    files["O"]["prb"] = "prb96_o.dat"
    files["Ar"]["acd"] = "acd89_ar.dat"
    files["Ar"]["scd"] = "scd89_ar.dat"
    files["Ar"]["plt"] = "plt89_ar.dat"
    files["Ar"]["prb"] = "prb89_ar.dat"
    return files

if __name__ == "__main__":
    if len(sys.argv)-1>0:
        main(sys.argv[1])
    else:
        main()
