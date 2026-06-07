import platform
import shutil

MAJOR_ARC_START = 5747
MAJOR_ARC_END = 407

PHENOTYPE = {
    "8251":  1, "8472":  1, "8473":  1,
    "12705": -1, "14798": -1, "16223": 1,
}
SNP_NAMES = ["8251", "8472", "8473", "12705", "14798", "16223"]
REF_SUFFIX = {"8251": "G", "8472": "C", "8473": "T", "12705": "C", "14798": "T", "16223": "C"}
ALT_SUFFIX = {"8251": "A", "8472": "T", "8473": "C", "12705": "T", "14798": "C", "16223": "T"}

COLOR_POS_BAR = "#a8e6cf"
COLOR_NEG_BAR = "#ffb3b3"
COLOR_PHEN_POS = "#89CFF0"
COLOR_PHEN_NEG = "#FADADD"

if platform.system() == 'Windows':
    HAS_SYSTEM_SORT = False
else:
    HAS_SYSTEM_SORT = shutil.which('sort') is not None
