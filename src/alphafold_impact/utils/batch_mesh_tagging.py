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

inp = eg3
print(inp)
skr_web_python_api_response_mti(inp)
