import os
import re
import time

#Global variables
SequenceFileName = "sequence.fasta"
StatusFileName = "status.phyreworks"
ConfigFileName = "config.phyreworks"

Phyre2StatusURL = "http://www.sbg.bio.ic.ac.uk/phyre2/webscripts/jobmonitor.cgi?jobid="
Phyre2SubmitURL = ""
Phyre2ResultURL = ""

class PhyreWorks(object):
    """Class to submit several molecules to Phyre2 and track progress"""
    def __init__(self, email_address, maxlen=1200, overlap=600, 
                 directory="PhyreWorks_data"):
        """Constructor for the PhyreWorks class"""
        self.email_address = email_address
        self.maxlen = maxlen
        self.overlap = overlap
        self.directory = directory
        self.molecules = []
        if not os.path.exists(self.full_path()):
            os.makedirs(self.full_path())
    @classmethod
    def from_file(cls, directory="PhyreWorks_data"):
        with open(os.path.join(directory, ConfigFileName),"r") as config_file:
            lines = [x.split(":\t") for x in config_file.read().split("\n")]
        email = lines[0][1]
        maxlen = int(lines[1][1])
        overlap = int(lines[2][1])
        new_works = cls(email, maxlen=maxlen, overlap=overlap)
        for molecule in new_works.subfolders():
            molecule_file = os.path.join(new_works.full_path(),molecule)
            PhyreMolecule.from_file(new_works, molecule_file)
        return new_works
    def save(self):
        with open(os.path.join(self.directory, ConfigFileName),"w") as config_file:
            print(" Email Address:\t" + str(self.email_address),
                  file=config_file)
            print("Maximum Length:\t" + str(self.maxlen),
                  file=config_file)
            print("  Overlap Size:\t" + str(self.overlap),
                  file=config_file)
        for molecule in self.molecules:
            molecule.save()
    def full_path(self):
        """Returns the PhyreWorks root"""
        return self.directory
    def from_fasta(self, filename):
        """Loads new sequences from the specified FASTA file"""
        return None
    def subfolders(self):
        return [x[1] for x in os.walk(self.full_path())][0]

class PhyreMolecule(object):
    """Class to submit a single molecule to Phyre2 (which may be split into
    fragments)"""
    def __init__(self, parent, sequence, name, create=True):
        """Constructor for the PhyreMolecule class"""
        self.parent = parent
        self.parent.molecules.append(self)
        self.sequence = sequence
        self.name = name
        self.machinename = sanitise(name)
        self.fragments = []
        if not os.path.exists(self.full_path()):
            os.makedirs(self.full_path())
        if create:
            # Generate fragments based from full sequence:
            full_length = len(self.sequence)
            num_fragments = int(np.ceil((full_length-self.parent.overlap) \
                                  / (self.parent.maxlen-self.parent.overlap)))
            for fragment_number in range(0,num_fragments+1):
                start = int(np.round(1 + (full_length - self.parent.maxlen) \
                                       * (fragment_number/num_fragments)))
                end = start + self.parent.maxlen - 1
                PhyreFragment(self, start, end)
    @classmethod
    def from_file(cls, parent, folder):
        with open(os.path.join(folder, SequenceFileName), "r") as sequence_file:
            lines = sequence_file.read().split("\n")
            human_name = lines[0].split(">")[1]
            sequence = "".join(lines[1:])
            new_molecule = cls(parent, sequence, human_name, create=False)
            for fragment_name in new_molecule.subfolders():
                fragment_folder = os.path.join(new_molecule.full_path(),
                                               fragment_name)
                PhyreFragment.from_file(new_molecule,fragment_folder)
            return new_molecule
    def save(self):
        save_path = self.full_path()
        sequence_path = os.path.join(save_path, SequenceFileName)
        with open(sequence_path,"w") as sequence_file:
            print(">"+self.name, file=sequence_file)
            print(self.sequence, file=sequence_file)
        for fragment in self.fragments:
            fragment.save()
    def full_path(self):
        """Returns the PhyreWorks root"""
        return os.path.join(self.parent.full_path(), self.machinename)
    def subfolders(self):
        return [x[1] for x in os.walk(self.full_path())][0]
    def submit(self):
        """Submit all sequences to Phyre2"""
        return None
    def status(self):
        """Check the status of all sequences"""
        return None
    def download(self):
        """Download all newly completed sequences"""
        return None
    def process(self):
        """Process all newly completed sequences"""
        return None
        
import numpy as np
class PhyreFragment(object):
    """Single Phyre2 submission"""
    def __init__(self, parent, start, end):
        """Constructor for the Phyre2 Submission"""
        self.start = start
        self.end = end
        self.parent = parent
        self.phyre_id = ""
        self.submitted = 0
        self.lastchecked = 0
        self.completed = 0
        self.downloaded = 0
        self.parent.fragments.append(self)
        if not os.path.exists(self.full_path()):
            os.makedirs(self.full_path())
    @classmethod
    def from_file(cls, parent, folder):
        with open(os.path.join(folder, StatusFileName),"r") as status_file:
            lines = [line.split(":\t") 
                     for line in status_file.read().split("\n")]
        new_fragment = cls(parent,int(lines[0][1]), int(lines[1][1]))
        new_fragment.phyre_id = lines[2][1]
        new_fragment.submitted = float(lines[3][1])
        new_fragment.completed = float(lines[4][1])
        new_fragment.downloaded = float(lines[5][1])
        new_fragment.lastchecked = float(lines[6][1])
        return new_fragment
            #key, value = line.split("\n")[0].split(":\t")
    def full_path(self):
        """Returns the PhyreWorks root"""
        return os.path.join(self.parent.full_path(), "_".join(["fragment",
                            str(self.start), str(self.end)]))
    def save(self):
        status_filename = os.path.join(self.full_path(), StatusFileName)
        with open(status_filename,"w") as status_file:
            print("First amino acid:\t" + str(self.start),
                  file=status_file)
            print(" Last amino acid:\t" + str(self.end),
                  file=status_file)
            print("   Phyre2 job id:\t" + str(self.phyre_id),
                  file=status_file)
            print("   Job submitted:\t" + str(self.submitted),
                  file=status_file)
            print("   Job completed:\t" + str(self.completed),
                  file=status_file)
            print("  Job downloaded:\t" + str(self.downloaded),
                  file=status_file)
            print("Job last checked:\t" + str(self.lastchecked),
                  file=status_file)
                  
    def sequence(self):
        """Gets amino acid sequence for this fragment"""
        self.parent.sequence[self.start-1:self.end]
    def padded_sequence(self):
        """Gets the amino acid sequence for this fragment padded with spaces
        for an alignment"""
        return " " * (self.start - 1) + self.sequence()
    def submit(self):
        return False
    def status(self):
        return False
    def download(self):
        return False
    def advance(self):
        return False
    def process(self):
        return False

def sanitise(dangerous):
    safe = re.sub('[^0-9a-zA-Z]+', '*', dangerous)
    return safe
    
        
        
bigseq = """MSNINNKDSSTEWNCKEDVGCVPPRRQNLNMERLDNENEDSVPDFMKKTFYLAAAGEGKK
LREKHDESCDEFCDAWNRSLADYKDIFQGKDMWNDGKYGEAKNHIKNAFGDMNNRKTMLN
EIEKGIKDETFSRENGLDVCKSQCEERSRDDTEDQFLRFFAEWEEEFCDGLNKHEEQLKS
CTKDINCDIKCSNFKDWLETKKDEYDIQSRVFEKKYANDNKSKHLNYLKEGMNKCKVKNP
EMVFKSGFANVAECRNLNVEGAGNKNSNNLKDLDSNSDKDGIVSESYKATKKNGESIMDR
IPKSFNKLFGYFSGSQEEEQKENDVSHRNNYDNILVDKFHRSSLLDKLDDRMFFDELNRD
NIMEEVLSKIPEPIIREAPKYVPKKPVPPQHIPRGDNVPRNIDVNGSKDEYSPETESANS
KIKPTYEENDEDKSKISIETSEIDRDKEPFRISEEKKVVKEDVQELENIEYEFDETFDFF
DEDAKRNIDDIRKDIQSQIMKSVENYNSEKEEFKRNIETQLIESGDGMNAGNYSSALQDN
STEIPTMVLVPGVLTVFLLTIIWVLVYKHSLIDRVHGTGKTKKEEKMKELESEEFPKEKY
NIEDMEETEKENEIEKIMDENKETEQTEEGNTEEFVQEKELNQETLENEIIVDHIKEEED
TRNVKEQESVLEESPIEELPVENNIGKINEEVEEFILNEIPLEEQVPKELPQEEIEEIVV
EELPIDEHLSSEETTVTEEDTFKGQLINEEKPVEEKSVSEEIPVEEKSVSEEIPVEEKSV
SEEIPVEEKSVSEEIPVEEKNVSEEIPVEEKNVSEEIPVEEKNVSEEIPVEEENVSEEIP
EGGIAIEDVPVDEETVTEEITVDEKIYDKLPNEIETVNEEMPVEDETLTEQISSEHERVP
EEIIEEKPFTEGEETESLTDEIVEEGVVTDDIPEEQIITEKVQEEEEFVTGELSEEDIIN
EKVQEEDESVTEELPEEDIINEKVQEEEESAYQEIVQDGSVTKDVEYKELVNDDVRDKEN
FVIEEDPFKGQLINEGLPVEEEFVTKELPVKEESVFEELTEEDQSVTKEIPVEEHSVFKE
VDEIESVTDEIVEEEGSVNEEVEEEVSVSGEVDETEYVTEEVEETELVNEEVLKEEGSAS
EKVVKEEGSASEKIVEEEGSVTEETASEEIVEEEGSVTEESSSEEIVEDEGSVTEESSSE
EIVEEEGSVIEESASEVIVEDEVSASKEIVEDEASVTEEVVEEEGAVSDEVQVTESVEDE
IINQGIVDEVIVEQEASVTDEIVKEDESVIEEIAVEETVTEVVEETKPMDEEIVDQGSVV
ENVEEKKTMDEEIVDQGSVVENVEEKKTMDEEIGDQGSVAEKVEEEELVTEEVIEREGSI
NEDIVKEASITEEVEQIEPVTGKVEKIESITDEIKEQLVPEDIKEEQLDFEEIVAQRASV
NDEKVEVASITEEVEEEEKSVSEEVLEEEGSTTEKVVKGSSTEVVVEEQGSVTENLLEEE
SASQGIVEKEEFVDEEDSVKDQNVFEKEGSVTEQLVEEEKGLINEDNEKEELITEMSEEI
KSVNEEIEETDLSTEEIIKQQGLATYEFVEEEKSLTDKLLEEESVTKEVGETELSTQEVV
DEKVSVTEEVIEEEKSVSEEVLEEGSATEEVVEEGSCTEVVVEEGSCTEVVVEEGSCTEI
VVEKGSDTEIVVEEGSATEIVVEEGSATEVVVEEGSATEEVVEEGSVSEEMLEEEGSATE
EVVEEGLSSDNVQKSKGVIENVGEIYSVKTAKDESMNEKIPLEKSSFVDDESFKGQGPTD
NVSVEDVNSEDIINEHTPLEETKIEELPTEYITTADIHTKGETETKYNLIYEKINEEVEK
AKFQEEKITENIPVERESVTEDIVQEPSLAQEVEQKESDTNEIEETKLANEKIIPEVSVT
ENVVEKEGLDTEEVLEEDESITEEIVEEEVSSSEEIVEEEESSSEEIVEEEESSSEEIVE
EEESSSEEIVEEEESSSEEIVEEEESSSEEIIEEVSSTEEVLEEEGSVTEEIVEEEVSTT
EEVKDIGSVSEEVLEEEGFGTEEFVGQQGSVIEEIVETESSTEKVLEDVGSNVEEIVQEE
GPVAQEIVHEEVSTTEKHDEVDRSTTEEIVEKVGSVSEEIIVEEVSASEEIVEEGSVTEE
VVEEEKLINEVGETESVTEEIVQKEVSDAEEVLGQEGSMNEEILEKESIVEEIVGPEGSV
TEEIVDHGSFAEEVKEEELVTEEAVQYEGSVTEEIKEEESITENEAIEESAFAEIIEEKG
PNTDEIVKEEGLDTEEIVNEVSVTDEVIEEEKLVNEQIVGEERSVTEKPVEVERSATEDL
VEEEASVTEKVSVHEGSTTEQILDESVAEEIVEEEVSVDDKIIEEEVSVDEVVEEEGSVI
EEIVEEEESVPEEILEEELSGSEEVLEDEWVTDAFMGQEGSVIEEIEEIVDGEGSITEEI
VEDGSANEKIVEEEPSRVEEVLGKEGFVIEEIIEEGSVIEQVEDTKTVSEKSEESSAIEE
VKEVKEEESISEKIVEKEESVTEEIVRQEESTTEKIVKDVSPTEDFVEQTDSVTEKVIEQ
EGSNTEVAEDVEEKESASDEHEQEDVSVNAQVTYEKKSVTKEIVDEVSRTEEIVEENGSV
TEGVDETGSVTEEIIEEATVTEEVVEDGSVTEEVVEDGSVIQEVVEDGSVTEEIVQENGS
VTEEIVEEEGSVNEEVEEEVSVSGEVDETEYVTEEVEEEGSVVEEIVEEEGSVSGEVDET
EYVTEEVEEEGSVVEEIVEEEGSVSGEVDETEYVTEEVEEEGSVVEEIVEEEGSVSGEVD
ETEYVTEEVEEEGSVVEEIVEEEGSVVEEIVEEEGSVSEVVDETELVNDEIVEQAPFTEE
VEEQVSVNDEIIEDASVAEAVEESESITESVSQEEETEKGFVIEKVEETGAVTEEIVQDG
LITEEILEESESVNGEIINKESDAEEILETEFLTEEVVGQAGSTSEEIVEEEGSVTKEVE
EKESVTEELVDEGSVTEELVDEGSVTEEVVEQGGSIAQEIVEEESATEEIIRDETNVEEV
LEKEGSATEEIVQDGSGTNDFVGKQGSVIEEVVEEEISTTEEKLKEEASAIEEFVEEESI
REDVLEESLVTENVVGQQESVTEEIVDGEGSFTEDIVEEEESVTEEIVVDEESVTKEIVE
DEELVTEEIVEDEGSFTEEIVEDEGSFTEEVIEERSLIEEVEDTETVAEKEEGSVIKEII
DEKSLTEKIVEEEKSVTEEVEEKESVKEEVEEQRLVVEEEGSATEGIVEDRLATEGIVDD
ILVTEEIVEDGLATDEFVEQQGSIIEEVLDDEGSVTEEIVEEEGSPNEEIVEGVSVIEED
DNIEPVSEEIVEGSVTEEMIKEGLENEVILDEDSITEEALEKEGSVSEEIVEEMGSLTEE
IVDEERSTSEDMIEEGSASEEIIQEESQVEEVVEEVSVIDEIVEEDELDTKEVVEEIEFN
TEEVVEHKEEEGSVAEEIVQEEKEGSVNEEIIEEVGSITEEMVEQDVSDNEEIVEERSVI
EEAEENVWIEKEVEEEGLDNEEVIDEEDSVSEQAEEEVYINEEILKESSDVEDVKVENEL
MNEEVNEETQSVAENNEEDKELDNYVVEETESVTEEVVVDEVPNSKEVQEIESIIEEIVE
DGLTTDDLVGQQGSVIEEVVEEVGSDSEEIVEEASITEEVEKKESVTEDILVEESVTGDI
LVEGSVTEEVVGEEKLVSEEIVTEEGSVAQEIVEEDAPATEEIDEIESVTEEVVEEEGPV
DEEIVQEEGSVTEEIIQGESKVEEVVEEQGSENEEIFVEEVSASQEIVQNESGTEEILEK
VSASQEIVQDGSVTEQIIEEQKPVTEEVVNEEESITHEIIQEESHVEKVVQQGSVAEEVV
ENPVSVTEEIVEKEGSVTEDIGQEGYVAEEIVEEEEFDNEEILEEESVAEEFVEEGFDNE
EIFVEQISDSEIVKEENSVNEEVLEEEGSYTEEILEEEGSYTEEILEEEGSFTEEILEEE
GSYNEEILEEEGSYNEEILEEEGSYNEEILEEEGSYNEEILEEEGSATDYFVGQGSDNEE
IIEEGSATDYFVGQGSDNEEIIEEGSATDYFVGQGSIIEEVLEEEGLDTEKVFENEGSAT
EFEETESFTEVVEETESVNENILEETSINEVQKIESITEDIKEQLVPEEIKEEQLDSEEI
KEEQLDSEEIKEEQLVPEEIKEEQLDSEEIKEEQLDSEEIKEEQLDSEEIKEEQLDSEEI
KEEQLVPEEIKEEQLDSEEIKEEQLVPEEIKEEQLDSEEIKEEQLDSEEIKEEQLDSEEI
KEEQLDSEEIKEQKGSVNEEVVEEEGSVTEEIKEQEESVNEEVLEEVEETESIKEEIVEG
GIATQEIIEEESDTKEVVEEEVIDSEKLVDAGSVTGEVMPEEVSVTDEVVEEGSTTEEVL
EEQKSVNEEVVEDGLTIDDFVGLQGSTTEEVVEEDGSAIEKILEEETATEEIVEKQVSVT
EDIVEKEGSVNEEIIEEASVAEEIIQGGSFTEEIVGQEESATEEVIDEEGLISNEIEEEE
EKSVTEEMIEEVEEVSVDDEVEEVSVAEEIVEEELVDDEILPEELSATEDVIEEVRSVTD
EIVQEESVCEEILEQEVSASEEYVDDKSVTDDFVGHERSVIQDVENTESVTEEIAEVDKS
VIEEAVEKQGSVTEEKVQEGVSAIEEIEELESVTEEIAEEDKSVIEEAVEKQGSVTEEIV
EEEELDTEEVLEDKSVTGDVVEQEGSGKDESEAKESFTEEVDELKSVKEEDQETEYISRE
IEEESATEQHSEQELSINKEVVETESLTKDIEEEKSTTQEILEETQSVNEEIVEEERDTD
EVLKEKVSPSEEVIEEQASTTEEFVEERSSTDEIVEVEDLFTEEVKEREGSVTEEIVEEG
SDTGEIVEEEGSDTEEILEEGSFNEEIVEEEGSITEEILQGSVTEEILQGSVTEEFVGQQ
GSVIEEIVETESAIEERVEEESATEEVDERESVTEVVEEEVSSSDEVVEGSIEEVIENEG
SVTEEILEHEVSADENFVGQAVSVIEEVEGTESVTEEVVEETESVSEEIVEVSPTENVVQ
QTDSVIEEVVEQKEGSFNEEIDIRELGDDGVEEREKISTEEVVGQDKSATGDVEEVSSTE
DEEEVSSTEGLEEVSSTEGLEEVSNTEDVEEVSSTEDVEEGSVAENVKETKSITEEVSVE
EDIITDKVSVEQEVMAEASVEENILTEVPVEEEIMTEKLSVEDKALNEKIMSEEEIVIED
GNVHEVVPAEVSVTEEIPGVEETTNNESHVKGENVVNEVVVDDNSVNDEIQFDDDSSIEI
YTVDSKDVFHKEENYDSFREEVRSDENIHIFRKKNENFVKKIENEKSSVGHVPTEYTSKE
NIAEEVPSHIMFKENVTEEVPKEVKYEENTVEEVFEEVTSKENIIEEAPGDINSEENIIE
EVPEVVTSEEIVVQNDVNNINTNIDHMFDFDLDDILKIPEAQRKLNKLRTLVDVHLDVIE
KMQRDEWKKNKRDFLYICLKEINKLSEDTLMKFYGSDSNKSANDNDTVMVIKILKDKWGT
GRIVDTIANSLNKSYNSLYHNLYIEMEKDLLINKTETFNKWSKQHWNKLDNWKEEKWFKL
FKRDLKIDMRNTYESDQNDEENNEESMNELSDELIKNNTSDNMRNEQKTLEKNKEYSNLS
KDLGLIEKQKIIWKSWIVKNVNNIENWFDEMWFKNVVNELKEKNDVSNVLQENASEASVE
NYISDVEDSKDMTNSINDSEKSKIVASSSKTTNEEYIILSRKDLIYNIIVMVHMMVLDQF
KYDELKYAKKRFLNRSIDKFIKEKKIKDKETVLDYYIDDIIKRIVEDTSHVNNIKEKAID
HYTSHDWFRLLRQGNKAQNSISQEVDILTEKYKKLISEEDNKKENNNNNDDEVKPELKDN
IKCEENVSFNILRRSNKKDQLPLEEKKNKNGDLNNTLTESDGKKININEAIEESKNGNED
KNSGNMEKSRKPRDRRTERKEKDQDLRVQLLQDYENIFEQIQNMENQNKDKNKKETKNIL
KTSIDLDKNLLREYKKEKNIINQSEKEFSVNEN"""

bigseq = "".join( bigseq.split() )

c = PhyreWorks.from_file()

#- PhyreWorks_data
#    - molecule1.tracker
#        (name)               (sequence)            (start)  (end)   (Phyre2ID)      (Submitted)         (Completed)
#        molecule1_fragment_1 AFIURIGBWGBWFNSDB...   13       120     39t32h329t3   2014-06-13 15:24:04
#    - results
#        - 39t32h329t3
#        - u2tj05484tj
#
#
#maxlen = 500
#overlap = 100
#for target in targets:
#    length = len(sequences[target].seq)
#    if length > maxlen:
#        for n in range(0, length//(maxlen - overlap)):
#            start = n * (maxlen - overlap)
#            end = n*(maxlen-overlap) + maxlen
#            fragment = str(sequences[target].seq)[start:end]
#            name = 'fragment_'+str(start)+'_'+str(end)+'_'+target
#            print( name,metaPrDOS_call(name, 'molecules@inbox.com', fragment ))
#            for i in range(0,100):
#                sleep(6.1)
#                print("waiting..." + str(i) + '%')