"""
This file contains functions for using the SKR web python API
for submitting batch jobs using the NLM Medical Text Indexer (MTI).

See here for installation steps: https://github.com/lhncbc/skr_web_python_api

Get UMLS API key from here: https://uts.nlm.nih.gov/uts/profile
"""

from skr_web_api import Submission


def skr_web_python_api_generic_batch(
    email: str,
    api_key: str,
    txt_file_name: str,
    cmd: str,
    cmdargs: str,
) -> str:
    """Submits NIH Generic Batch job using the SKR Web API

    Args:
        email (str): User's email address.
        api_key (str): UTS API key.
        txt_file_name (str): Name of the text file to submit to be tagged.
            Each record should be on a new line and be in the format {id}|{text}.
            See example here: https://lhncbc.nlm.nih.gov/ii/about/help/single_line_ID_sample.txt
        cmd (str): Command to be executed.
        cmdargs (str): Arguments for the command. A full list of args can be found:
            https://lhncbc.nlm.nih.gov/ii/tools/MTI/help_info.html
    """
    serviceurl = (
        "https://ii.nlm.nih.gov/cgi-bin/II/UTS_Required/API_batchValidationII.pl"
    )

    inst = Submission(email, api_key)
    inst.set_serviceurl(serviceurl)
    inst.init_generic_batch(cmd, cmdargs)
    inst.set_batch_file(txt_file_name)
    inst.form["SingLinePMID"] = True
    response = inst.submit()
    return response.content.decode().replace("NOT DONE LOOP\n", "")
