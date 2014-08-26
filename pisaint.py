"""
PISA Interface Library
By Tom Newport, 2014

This library allows PDB ePISA Interface data to be downloaded, accessed and
searched. It may be imported as a python module or called via command line
arguments.

Run python pisaint.py help for command line help.
"""

import os
import xml.etree.ElementTree as ET
import urllib.request
import pandas
import numpy as np
import sys

class PisaConnection(object):
    """Class for downloading, caching and searching PISA files. PISA files are
    stored and cached in a folder."""
    def __init__(self, directory="PisaInterface_data"):
        """Constructor method. The cache directory can be changed using the
        directory keyword."""
        self.directory = directory
        if not os.path.exists(directory):
            os.makedirs(directory)
    def filename(self, pdb_id):
        """Returns the filename where the PISA result for the specified
        PDB ID should be stored."""
        return os.path.join(self.directory, self.sanitise_pdb(pdb_id)+
                            '.interfaces.pisa')
    def file_exists(self, pdb_id):
        """Checks if the PISA result for the specified PDB ID has been
        downloaded."""
        return os.path.isfile(self.filename(self.sanitise_pdb(pdb_id)))
    def get_file(self, pdb_id, force_download=False):
        """Returns the filename of the PISA result for the specified PDB ID
        after making sure it has been downloaded.
        If force_download is set to true, the file will be downloaded again.
        """
        filename = self.filename(pdb_id)
        if self.file_exists(pdb_id) and not force_download:
            return filename
        try:
            url = "http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?"
            url += self.sanitise_pdb(pdb_id)
            url_request = urllib.request.urlopen(url)
        except:
            return False
        with open(filename, "w+b") as pisa_file:
            pisa_file.write(url_request.read())
        return filename
    def sanitise_pdb(self, pdb_id):
        """Standardises the pdb code (converts to lower case)."""
        pdb_id = pdb_id.lower()
        return pdb_id
    def get_pisa(self, pdb_id):
        """Returns a parsed PISA file for the specified PDB Code. The file
        will be downloaded if it does not already exist."""
        return PisaOutput(self.get_file(pdb_id), self.sanitise_pdb(pdb_id))
    def get_pit_cache_name(self, pdb_id, css_min):
        """Returns the filename of the PISA interface table cache with a
        specified minimum CSS."""
        return os.path.join(self.directory, self.sanitise_pdb(pdb_id)+'.cache.'
                            +str(css_min)+'.interfaces.pisa')
    def pit_cache_exists(self, pdb_id, css_min):
        """Checks if a PISA interface table cache file exists."""
        return os.path.isfile(self.get_pit_cache_name(pdb_id, css_min))
    def get_pit(self, pdb_id, css_min=0.75):
        """Returns the PISA interface table for the specified PDB ID at the
        specified css_min."""
        if not self.pit_cache_exists(pdb_id, css_min):
            resulting_dataframe = self.get_pisa(pdb_id).by_residues(css_min)
            resulting_dataframe.to_csv(self.get_pit_cache_name(
                pdb_id, css_min))
        else:
            resulting_dataframe = pandas.DataFrame.from_csv(
                self.get_pit_cache_name(pdb_id, css_min))
        return resulting_dataframe
    def search(self, pdb_id, chain, start, end, css_min=0.75):
        """Returns a table summarising the interfaces between the specified
        amino acids and other chains or molecules.

        pdb_id - A PDB code.
        chain - The chain of the molecule specified by the PDB code.
        start - The first amino acid (inclusive).
        end - The last amino acid (inclusive).
        css_min - Minimum CSS for an interface to be considered.

        Returns a dataframe composed of:
            interface - Interface ID
            bsa_proportion - Proportion of the interface's buried surface area
            represented by specified amino acids
            to_chain - The chain that is bound to
            bsa_sum - The total buried surface area represented by specified
            amino acids.
            chain - should be identical to the "chain" argument.

            The total buried surface area of the interface is equal to
            bsa_sum/bsa_proportion

        """
        resulting_dataframe = self.get_pit(pdb_id, css_min=css_min)
        if len(resulting_dataframe) == 0:
            return pandas.DataFrame()
        resulting_dataframe = resulting_dataframe[
            resulting_dataframe["chain"] == chain]
        resulting_dataframe = resulting_dataframe[
            resulting_dataframe["seq_num"] >= start]
        resulting_dataframe = resulting_dataframe[
            resulting_dataframe["seq_num"] <= end]
        resulting_dataframe_sum = resulting_dataframe.groupby(
            "interface").sum()
        result = pandas.DataFrame()
        result["bsa_proportion"] = resulting_dataframe_sum["bsa_contrib"]
        result["bsa_sum"] = resulting_dataframe_sum["bsa"]
        resulting_dataframe_summary = resulting_dataframe.drop_duplicates(
            subset="interface")
        result = pandas.merge(result, resulting_dataframe_summary,
                              left_index=True, right_on="interface")
        return result[["interface", "bsa_proportion", "to_chain", "bsa_sum",
                       "chain"]]



class PisaOutput(object):
    """Class to parse the PISA Interface XML file returned from PDB ePISA"""
    def __init__(self, xmlfile, pdb_id):
        """Initialise the class using the filename of the xml file and the
        PDB code of the molecule of interest."""
        self.xmlfile = ET.parse(xmlfile)
        self.status = False
        self.pdb_id = pdb_id
        for pdb_entry in self.xmlfile.getroot().findall('pdb_entry'):
            if pdb_entry.find('pdb_code').text.lower() == self.pdb_id.lower():
                self.xmlroot = pdb_entry
                self.status = True
                break
        if not self.status:
            print("ERROR")
        self.interfaces = {}
        for interface in self.xmlroot.findall("interface"):
            interface_id = interface.find("id").text
            self.interfaces[interface_id] = PisaInterface(interface)
    def css_min(self, css_min):
        """Returns a dict of interfaces with a css greater than css_min"""
        return {key: value for (key, value) in self.interfaces.items()
                if value.css >= css_min}
    def by_residues(self, css_min, all_residues=False):
        """Returns a Pandas dataframe of all residues for all interfaces
        greater than css_min"""
        filtered = self.css_min(css_min)
        resulting_dataframes = []
        found = False
        for interface in filtered.values():
            resulting_dataframes.append(interface.molres_to_dataframe())
            if len(resulting_dataframes[-1]) > 0:
                found = True
        if not found:
            return pandas.DataFrame()
        resulting_dataframe = pandas.concat(resulting_dataframes)
        resulting_dataframe["pdb_code"] = self.pdb_id
        if all_residues == False:
            resulting_dataframe = resulting_dataframe[
                resulting_dataframe["bsa_contrib"] > 0]
        return resulting_dataframe

class PisaInterface(object):
    """Class to access individual PISA interfaces. Takes the XML root of the
    interface
    """
    def __init__(self, root):
        self.id = root.find("id").text
        self.type = root.find("type").text
        self.n_occ = root.find("n_occ").text
        self.int_area = root.find("int_area").text
        self.int_solv_en = root.find("int_solv_en").text
        self.pvalue = root.find("pvalue").text
        self.stab_en = root.find("stab_en").text
        self.css = float(root.find("css").text)
        self.overlap = root.find("overlap").text
        self.x_rel = root.find("x-rel").text
        self.fixed = root.find("fixed").text
        self.h_bonds = []
        self.salt_bridges = []
        self.ss_bonds = []
        self.cov_bonds = []
        for bond in root.find("h-bonds").findall("bond"):
            self.h_bonds.append(PisaBond(bond))
        for bond in root.find("salt-bridges").findall("bond"):
            self.salt_bridges.append(PisaBond(bond))
        for bond in root.find("ss-bonds").findall("bond"):
            self.ss_bonds.append(PisaBond(bond))
        for bond in root.find("cov-bonds").findall("bond"):
            self.cov_bonds.append(PisaBond(bond))
        self.molecules = {}
        for molecule in root.findall("molecule"):
            self.molecules[molecule.find("id").text] = PisaMolecule(molecule)
    def molres_to_dataframe(self):
        """Returns a dataframe of all the residues in all the molecules
        involved in the interface"""
        all_chains = [x.chain_id for x in self.molecules.values()]
        molecule_dataframes = []
        found = False
        for molecule in self.molecules.values():
            molecule_dataframe = molecule.residues_to_dataframe()
            if len(molecule_dataframe) > 0:
                found = True
            chain_id = list(molecule_dataframe["chain"])[0]
            chain_list = [x for x in all_chains]
            chain_list.remove(chain_id)
            molecule_dataframe["to_chain"] = ", ".join(chain_list)
            molecule_dataframes.append(molecule_dataframe)
        if not found:
            return pandas.DataFrame()
        resulting_dataframe = pandas.concat(molecule_dataframes)
        resulting_dataframe["interface"] = self.id
        return resulting_dataframe


class PisaBond(object):
    """Class for accessing ePISA bonds. Takes the XML root of the bond."""
    def __init__(self, root):
        """Constructor for the class"""
        self.chain_1 = root.find("chain-1").text
        self.res_1 = root.find("res-1").text
        self.seqnum_1 = root.find("seqnum-1").text
        self.inscode_1 = root.find("inscode-1").text
        self.atname_1 = root.find("atname-1").text
        self.chain_2 = root.find("chain-2").text
        self.res_2 = root.find("res-2").text
        self.seqnum_2 = root.find("seqnum-2").text
        self.inscode_2 = root.find("inscode-2").text
        self.atname_2 = root.find("atname-2").text

class PisaMolecule(object):
    """Class for accessing ePISA molecules. Takes the XML root of the
    molecule"""
    def __init__(self, root):
        """Constructor for the class"""
        self.id = root.find("id").text
        self.chain_id = root.find("chain_id").text
        self.mol_class = root.find("class").text
        self.symop_no = root.find("symop_no").text
        self.cell_i = np.array([float(x) for x in [root.find("cell_i").text,
                                                   root.find("cell_j").text,
                                                   root.find("cell_k").text]])
        self.matrix = np.zeros([3, 4])
        self.matrix[0, :] = [float(x) for x in [root.find("rxx").text,
                                                root.find("rxy").text,
                                                root.find("rxz").text,
                                                root.find("tx").text]]
        self.matrix[1, :] = [float(x) for x in [root.find("ryx").text,
                                                root.find("ryy").text,
                                                root.find("ryz").text,
                                                root.find("ty").text]]
        self.matrix[2, :] = [float(x) for x in [root.find("rzx").text,
                                                root.find("rzy").text,
                                                root.find("rzz").text,
                                                root.find("tz").text]]
        self.int_natoms = root.find("int_natoms").text
        self.int_nres = root.find("int_nres").text
        self.int_area = root.find("int_area").text
        self.int_solv_en = root.find("int_solv_en").text
        self.pvalue = root.find("pvalue").text
        self.residues = {}
        for residue in root.find("residues").findall("residue"):
            new_residue = PisaResidue(residue)
            self.residues[new_residue.ser_no] = new_residue
    def residues_to_dict(self):
        """Returns a dictionary of individual residues."""
        return {key:value.to_dict() for (key, value) in self.residues.items()}
    def residues_to_dataframe(self):
        """Returns a dataframe of individual residues."""
        resulting_dataframe = pandas.DataFrame.from_dict(
            list(self.residues_to_dict().values()))
        resulting_dataframe = resulting_dataframe.sort(["seq_num"],
                                                       ascending=[True])
        resulting_dataframe["chain"] = self.chain_id
        resulting_dataframe["mol_bsa"] = resulting_dataframe["bsa"].sum()
        resulting_dataframe["bsa_contrib"] = resulting_dataframe["bsa"]/\
        resulting_dataframe["mol_bsa"]
        return resulting_dataframe


class PisaResidue(object):
    """Class to access residues in the ePISA XML File. Takes the XML root
    of the residue"""
    def __init__(self, root):
        self.ser_no = root.find("ser_no").text
        self.name = root.find("name").text
        self.seq_num = int(root.find("seq_num").text)
        self.ins_code = root.find("ins_code").text
        self.bonds = root.find("bonds").text
        self.asa = float(root.find("asa").text)
        self.bsa = float(root.find("bsa").text)
        self.solv_en = root.find("solv_en").text
    def to_dict(self):
        """Returns the attributes of the residue as a dict"""
        return self.__dict__

def print_help():
    """Prints the command line documentation for the PISA Interface Search
    Library"""
    helpstring = '''
pisaint command line
====================
Tom Newport 2014

Any problems, questions or suggestions, email me (tom.newport@gmail.com)

PISA Interface Library (pisaint) may be called in one of 3 ways:

1: Simple - python pisaint.py [pdb_id] [chain] [start] [end]
   Example: >> python pisaint.py 2y26 A 10 300
            This will print all PISA interfaces present between
            residues 10 and 300 of chain A of pdb molecule 2y26

2: Batch -  python pisaint.py [pdb_id] [chain] [start] [end] [output]
            This will output to the specified file rather than print
            to the command window.

3: Advanced-python pisaint.py [pdb_id] [chain] [start] [end] [output] [cache]
            This allows the PISA data cache to be specified.


'''
    print(helpstring)

def pint_cmd():
    """Command line interface"""
    if len(sys.argv) not in [5, 6, 7]:
        print("\nError: Wrong number of arguments. \n")
        print_help()
    if len(sys.argv) == 5:
        pisa_main = PisaConnection()
        pisa_search = pisa_main.search(sys.argv[1], sys.argv[2],
                                       int(sys.argv[3]), int(sys.argv[4]))
        if len(pisa_search) > 0:
            print(pisa_search)
        else:
            print("No results found.")
    if len(sys.argv) == 6:
        pisa_main = PisaConnection()
        pisa_main.search(sys.argv[1], sys.argv[2], int(sys.argv[3]),
                         int(sys.argv[4])).to_csv(sys.argv[5])
    if len(sys.argv) == 7:
        pisa_main = PisaConnection(directory=sys.argv[6])
        pisa_main.search(sys.argv[1], sys.argv[2], int(sys.argv[3]),
                         int(sys.argv[4])).to_csv(sys.argv[5])

if __name__ == "__main__":
    pint_cmd()
