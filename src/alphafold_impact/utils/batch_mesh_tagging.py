"""
This file contains functions for using the SKR web python API
for submitting batch jobs using the NLM Medical Text Indexer (MTI).

See here for installation steps: https://github.com/lhncbc/skr_web_python_api

Get UMLS API key from here: https://uts.nlm.nih.gov/uts/profile
"""
from pathlib import Path
from kedro.config import OmegaConfigLoader
from kedro.framework.project import settings
from skr_web_api import Submission

conf_loader = OmegaConfigLoader(conf_source=str(Path.cwd() / settings.CONF_SOURCE))
CREDS = conf_loader["credentials"]
EMAIL = CREDS["email"]
API_KEY = CREDS["umls_api_key"]


def skr_web_python_api_response_mti(input_text: str) -> None:
    """Submits a given text input to the NLM Medical Text Indexer (MTI)
    via the SKR Web Python API and prints the response."""
    inst = Submission(EMAIL, API_KEY)
    inst.init_mti_interactive(input_text)
    response = inst.submit()
    print(f"response status: {response.status_code}")
    print(f"content:\n {response.content.decode()}")


def skr_web_python_api_response_gb() -> None:
    """Submits generics batch request with the given text input
    via the SKR Web Python API and prints the response."""
    # This reads in the file as 3 parts but returns
    # [0]: *** ERROR *** ERROR *** ERROR ***
    # [1]: *** ERROR *** ERROR *** ERROR ***
    # [2]: *** ERROR *** ERROR *** ERROR ***
    inst = Submission(EMAIL, API_KEY)
    inst.init_generic_batch("MTI", "-D")
    inst.set_batch_file("src/alphafold_impact/utils/docs.medline.txt")
    response = inst.submit()
    print(f"response status: {response.status_code}")
    print(f"content:\n {response.content.decode()}")


# For testing different file formats:
eg0 = "A spinal tap was performed and oligoclonal bands were \
detected in the cerebrospinal fluid.\n"

eg1 = """direct fall
shoulder

rotator cuff injury

out Complaining

right arm pain

numbness

zoster scarring

carpal tunnel Syndrome

diabetic

lipids

elevated so

factor
"""

eg2 = """direct fall
shoulder
rotator cuff injury
out Complaining
right arm pain
numbness
zoster scarring
carpal tunnel Syndrome
diabetic
lipids
elevated so
factor
"""

eg3 = """97479605|Higher neonatal cerebral blood flow correlates with worse childhood neurologic outcome.
98018928|Bupivacaine inhibition of L-type calcium current in ventricular cardiomyocytes of hamster.
"""

eg4 = """
PMID- 12477674
OWN - NLM
STAT- MEDLINE
DCOM- 20030507
LR  - 20131121
IS  - 0832-610X (Print)
IS  - 0832-610X (Linking)
VI  - 49
IP  - 10
DP  - 2002 Dec
TI  - Oxygen in air (FiO2 0.4) improves gas exchange in young healthy patients during
      general anesthesia.
PG  - 1040-3
AB  - PURPOSE: One hundred percent O(2) is used routinely for preoxygenation and
      induction of anesthesia. The higher the O(2) concentration the faster is the
      development of atelectasis, an important cause of impaired pulmonary gas exchange
      during general anesthesia (GA). We evaluated the effect of ventilation with 0.4
      FiO(2) in air, 0.4 FiO(2) in N(2)O and 100% O(2) following intubation on the
      development of impaired gas exchange. METHODS: Twenty-seven patients aged 18-40
      yr, undergoing elective laparoscopic cholecystectomy were administered 100% O(2) 
      for preoxygenation (three minutes) and ventilation by mask (two minutes).
      Following intubation these patients were randomly divided into three groups of
      nine each and ventilated either with 0.4 FiO(2) in air, 0.4 FiO(2) in N(2)O or
      100% O(2). Arterial blood gases were obtained before preoxygenation and 30 min
      following intubation for PaO(2) analysis. Subsequently PaO(2)/FiO(2) ratios were 
      calculated. Results were analyzed with Student's t test and one-way ANOVA. P
      value of &lt; or = 0.05 was considered significant. RESULTS: Ventilation of the
      lungs with O(2) in air (FiO(2) 0.4) significantly improved the PaO(2)/FiO(2)
      ratio from baseline, while 0.4 FiO(2) in N(2)O or 100% O(2) worsened the ratio
      (558 +/- 47 vs 472 +/- 28, 365 +/- 34 vs 472 +/- 22 and 351 +/- 23 vs 477 +/- 28 
      respectively; P &lt; 0.05). CONCLUSION: Ventilation of lungs with O(2) in air
      (FiO(2) 0.4) improves gas exchange in young healthy patients during GA.
FAU - Agarwal, Anil
AU  - Agarwal A
AD  - Department of Anesthesia, Sanjay Gandhi Post Graduate Institute of Medical
      Sciences, Lucknow, India. aagarwal@sgpgi.ac.in
FAU - Singh, Prabhat K
AU  - Singh PK
FAU - Dhiraj, Sanjay
AU  - Dhiraj S
FAU - Pandey, Chandra M
AU  - Pandey CM
FAU - Singh, Uttam
AU  - Singh U
LA  - eng
PT  - Clinical Trial
PT  - Journal Article
PT  - Randomized Controlled Trial
PL  - United States
TA  - Can J Anaesth
JT  - Canadian journal of anaesthesia = Journal canadien d'anesthesie
JID - 8701709
RN  - S88TT14065 (Oxygen)
SB  - IM
MH  - Adult
MH  - Air/*analysis
MH  - *Anesthesia, General
MH  - Cholecystectomy, Laparoscopic
MH  - Female
MH  - Humans
MH  - Male
MH  - Oxygen/*analysis
MH  - Pulmonary Gas Exchange
MH  - *Respiration, Artificial
EDAT- 2002/12/13 04:00
MHDA- 2003/05/08 05:00
CRDT- 2002/12/13 04:00
PHST- 2002/12/13 04:00 [pubmed]
PHST- 2003/05/08 05:00 [medline]
PHST- 2002/12/13 04:00 [entrez]
AID - 10.1007/BF03017898 [doi]
PST - ppublish
SO  - Can J Anaesth. 2002 Dec;49(10):1040-3. doi: 10.1007/BF03017898.

PMID- 25972241
TI  - Complement Peptide C3a Promotes Astrocyte Survival in Response to Ischemic Stress.
AB  - Astrocytes are the most numerous cells in the central nervous system with a range of homeostatic and regulatory functions. Under normal conditions as well as after ischemia, astrocytes promote neuronal survival. We have previously reported that  the complement-derived peptide C3a stimulates neuronal differentiation of neural  progenitor cells and protects the immature brain tissue against hypoxic-ischemic  injury. Here, we studied the effects of C3a on the response of mouse cortical astrocytes to ischemia. We have found that chemical ischemia, induced by combined inhibition of oxidative phosphorylation and glycolysis, upregulates the expression of C3a receptor in cultured primary astrocytes. C3a treatment protected wild-type but not C3a receptor-deficient astrocytes from cell death induced by chemical ischemia or oxygen-glucose deprivation by reducing ERK signaling and caspase-3 activation. C3a attenuated ischemia-induced upregulation  of glial fibrillary acidic protein; however, the protective effects of C3a were not dependent on the presence of the astrocyte intermediate filament system. Pre-treatment of astrocytes with C3a during recovery abrogated the ischemia-induced neuroprotective phenotype of astrocytes. Jointly, these results  provide the first evidence that the complement peptide C3a modulates the response of astrocytes to ischemia and increases their ability to cope with ischemic stress.
PT  - Journal Article
PT  - Research Support, Non-U.S. Gov't
TA  - Mol Neurobiol
JT  - Molecular neurobiology
JID - 8900963

PMID- 27752159
TI  - [Trends in Gleason scores of Chinese prostate carcinoma from 1995 to 2014].
AB  - OBJECTIVE: To assess the changing trends in Gleason score (GS) of Chinese prostate carcinoma (PCa) from January 1995 to December 2014. METHODS: In the study, 875 patients admitted to hospital from January 1995 to December 2004 (1995-2004) and from January 2005 to December 2014 (2005-2014) were divided into  two groups. The mean levels and proportions of GS, primary and secondary grades were studied. The patients were divided into four groups according to age: <60, 60-69, 70-79 and >/=80 years. Types of specimen included needle biopsy (NB), transurethral resection of the prostate (TURP) and radical prostatectomy (RP). Histological types were made up by acinar carcinoma and other types (including atrophic, pseudohyperplastic, foam, signet ring cell and ductal carcinoma, and so on). The total prostate-specific antigen (tPSA) involved groups of <20.0 mug/L and >/=20.0 mug/L. We observed the mean levels and proportions of GS in age, types of specimen, histological types and total prostate-specific antigen in different periods, and used SPSS 17.0 software for statistical analysis. RESULTS: Compared with 1995-2004, the mean levels of GS, primary and secondary grades decreased 0.32 (P=0.003), 0.19 (P=0.001) and 0.12 (P=0.016) in 2005-2014, respectively. The proportions of </=6 in GS increased 10.9% (P=0.003), and >/=8 decreased 14.0% (P<0.001). The difference of GS 7 was not statistically significant. In the primary grade, the ratio of grades</=3 increased 12.8% (P=0.001), and grade 4 decreased 7.4% (P=0.037), grade 5 decreased 5.5% (P=0.007). The ratio of secondary grades </=3 increased 7.6% (P=0.037). The difference of grades 4 and 5 was not statistically significant. CONCLUSION: GS in Chinese patients with PCa showed a downward trend, which is one of the notable features in the past 20 years in China. The types of specimen and age are important factors in GS, while the histological types and tPSA have less impact on the GS.
PT  - Journal Article
PT  - English Abstract
TA  - Beijing Da Xue Xue Bao
JT  - Beijing da xue xue bao. Yi xue ban = Journal of Peking University. Health sciences
JID - 101125284
"""

# inp = eg4
# print(inp)
# skr_web_python_api_response_mti(inp)


skr_web_python_api_response_gb()
