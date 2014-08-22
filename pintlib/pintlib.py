import os
import xml.etree.ElementTree as ET
import urllib.request
import pandas
import numpy as np
import sys

class PisaConnection:
	def __init__(self, directory = "pisa_interface_data"):
		self.directory = directory
		if not os.path.exists(directory):
			os.makedirs(directory)
	def filename(self, pdb_id):
		return os.path.join(self.directory,self.sanitise_pdb(pdb_id)+'.interfaces.pisa')
	def file_exists(self, pdb_id):
		return os.path.isfile(self.filename(self.sanitise_pdb(pdb_id)))
	def get_file(self,pdb_id,force_download = False):
		filename = self.filename(pdb_id)
		if self.file_exists(pdb_id) and not force_download:
			return filename
		try:
			url = "http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?"+self.sanitise_pdb(pdb_id)
			u = urllib.request.urlopen(url)
		except:
			return False
		with open(filename,"b+w") as f:
			f.write(u.read())
		return filename
	def sanitise_pdb(self, pdb_id):
		pdb_id = pdb_id.lower()
		return pdb_id
	def get_pisa(self, pdb_id):
		return Pisa_output(self.get_file(pdb_id), self.sanitise_pdb(pdb_id))
	def get_pit_cache_name(self, pdb_id, css_min):
		return os.path.join(self.directory,self.sanitise_pdb(pdb_id)+'.cache.'+str(css_min)+'.interfaces.pisa')
	def pit_cache_exists(self, pdb_id, css_min):
		return os.path.isfile(self.get_pit_cache_name(pdb_id, css_min))
	def get_pisa_interface_table(self, pdb_id,css_min = 0.75):
		if not self.pit_cache_exists(pdb_id, css_min):
			df = self.get_pisa(pdb_id).by_residues(css_min)
			df.to_csv(self.get_pit_cache_name(pdb_id, css_min))
		else:
			df = pandas.DataFrame.from_csv(self.get_pit_cache_name(pdb_id, css_min))
		return df
	def search(self, pdb_id, chain, start, end, css_min = 0.75):
		df = self.get_pisa_interface_table(pdb_id, css_min = css_min)
		df = df[df["chain"] == chain]
		df = df[df["seq_num"] >= start]
		df = df[df["seq_num"] <= end]
		df_sum = df.groupby("interface").sum()
		result = pandas.DataFrame()
		result["bsa_proportion"] = df_sum["bsa_contrib"]
		result["bsa_sum"] = df_sum["bsa"]
		df_summary = df.drop_duplicates(subset="interface")
		result = pandas.merge(result,df_summary, left_index = True, right_on = "interface")
		return result[["interface","bsa_proportion","to_chain","bsa_sum","chain"]]



class Pisa_output:
	def __init__(self, xmlfile, pdb_id):
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
			self.interfaces[interface_id] = Pisa_interface(interface)
	def css_min(self,css_min, chains = "", all_residues = False):
		return {key: value for (key, value) in self.interfaces.items() if value.css >= css_min}
	def by_residues(self,css_min, all_residues = False):
		filtered = self.css_min(css_min)
		dfs = []
		found = False
		for interface in filtered.values():
			dfs.append(interface.molres_to_dataframe())
			if len(dfs[-1]) > 0:
				found = True
		if not found:
			return pandas.DataFrame()
		df = pandas.concat(dfs)
		df["pdb_code"] = self.pdb_id
		if all_residues == False:
			df = df[df["bsa_contrib"] > 0]
		return df
			
class Pisa_interface:
	def __init__(self, root):
		self.id		  = root.find("id"		 ).text
		self.type		= root.find("type"	   ).text
		self.n_occ	   = root.find("n_occ"	  ).text
		self.int_area	= root.find("int_area"   ).text
		self.int_solv_en = root.find("int_solv_en").text
		self.pvalue	  = root.find("pvalue"	 ).text
		self.stab_en	 = root.find("stab_en"	).text
		self.css		 = float(root.find("css").text)
		self.overlap	 = root.find("overlap"	).text
		self.x_rel	   = root.find("x-rel"	  ).text
		self.fixed	   = root.find("fixed"	  ).text
		self.h_bonds = []
		self.salt_bridges = []
		self.ss_bonds = []
		self.cov_bonds = []
		for bond in root.find("h-bonds").findall("bond"):
			self.h_bonds.append(Pisa_bond(bond))
		for bond in root.find("salt-bridges").findall("bond"):
			self.salt_bridges.append(Pisa_bond(bond))
		for bond in root.find("ss-bonds").findall("bond"):
			self.ss_bonds.append(Pisa_bond(bond))
		for bond in root.find("cov-bonds").findall("bond"):
			self.cov_bonds.append(Pisa_bond(bond))
		self.molecules = {}
		for molecule in root.findall("molecule"):
			self.molecules[molecule.find("id").text] = Pisa_molecule(molecule)
	def molres_to_dataframe(self):
		all_chains = [x.chain_id for x in self.molecules.values()]
		d = []
		found = False;
		for molecule in self.molecules.values():
			mdf = molecule.residues_to_dataframe()
			if len(mdf) > 0:
				found = True
			chain_id = list(mdf["chain"])[0]
			chain_list = [x for x in all_chains]
			chain_list.remove(chain_id)
			mdf["to_chain"] = ",".join(chain_list)
			d.append(mdf)
		if not found:
			return pandas.DataFrame()
		df = pandas.concat(d)
		df["interface"] = self.id
		return df
		
			
class Pisa_bond:
	def __init__(self,root):
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
		
class Pisa_molecule:
	def __init__(self,root):
		self.id = root.find("id").text
		self.chain_id = root.find("chain_id").text
		self.mol_class = root.find("class").text
		self.symop_no = root.find("symop_no").text
		self.cell_i = np.array([float(x) for x in [root.find("cell_i").text,root.find("cell_j").text,root.find("cell_k").text]])
		self.matrix = np.zeros([3,4])
		self.matrix[0,:] = [float(x) for x in [root.find("rxx").text,root.find("rxy").text,root.find("rxz").text,root.find("tx").text]]
		self.matrix[1,:] = [float(x) for x in [root.find("ryx").text,root.find("ryy").text,root.find("ryz").text,root.find("ty").text]]
		self.matrix[2,:] = [float(x) for x in [root.find("rzx").text,root.find("rzy").text,root.find("rzz").text,root.find("tz").text]]
		self.int_natoms = root.find("int_natoms").text
		self.int_nres = root.find("int_nres").text
		self.int_area = root.find("int_area").text
		self.int_solv_en = root.find("int_solv_en").text
		self.pvalue = root.find("pvalue").text
		self.residues = {}
		for residue in root.find("residues").findall("residue"):
			new_residue = Pisa_residue(residue)
			self.residues[new_residue.ser_no] = new_residue
	def residues_to_dict(self):
		return {key:value.to_dict() for (key, value) in self.residues.items()}
	def residues_to_dataframe(self):
		df = pandas.DataFrame.from_dict(list(self.residues_to_dict().values()))
		df = df.sort(["seq_num"], ascending = [ True ])
		df["chain"] = self.chain_id
		df["mol_bsa"] = df["bsa"].sum()
		df["bsa_contrib"] = df["bsa"]/df["mol_bsa"]
		return df
			

class Pisa_residue:
	def __init__(self,root):
		self.ser_no = root.find("ser_no").text
		self.name = root.find("name").text
		self.seq_num = int(root.find("seq_num").text)
		self.ins_code = root.find("ins_code").text
		self.bonds = root.find("bonds").text
		self.asa = float(root.find("asa").text)
		self.bsa = float(root.find("bsa").text)
		self.solv_en = root.find("solv_en").text
	def to_dict(self):
		return self.__dict__

def print_help():
	help = '''
pintlib command line
====================
Tom Newport 2014

Any problems, questions or suggestions, email me (tom.newport@gmail.com)

PISA Interface Library (pintlib) may be called in one of 3 ways:

1: Simple - python pintlib.py [pdb_id] [chain] [start] [end]
   Example: >> python pintlib.py 2y26 A 10 300
            This will print all PISA interfaces present between
            residues 10 and 300 of chain A of pdb molecule 2y26

2: Batch -  python pintlib.py [pdb_id] [chain] [start] [end] [output]
            This will output to the specified file rather than print 
            to the command window.

3: Advanced-python pintlib.py [pdb_id] [chain] [start] [end] [output] [cache]
            This allows the PISA data cache to be specified.


'''
	print(help)

if __name__ == "__main__":
	if len(sys.argv) not in [5,6,7]:
		print("\nError: Wrong number of arguments. \n")
		print_help();
	if len(sys.argv) == 5:
		pisa_main = PisaConnection()
		print(pisa_main.search(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])))
	if len(sys.argv) == 6:
		pisa_main = PisaConnection()
		pisa_main.search(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])).to_csv(sys.argv[5])
	if len(sys.argv) == 7:
		pisa_main = PisaConnection(directory = sys.argv[6])
		pisa_main.search(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])).to_csv(sys.argv[5])

