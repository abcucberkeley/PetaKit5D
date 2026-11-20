"""Setup of the detection benchmark (C02): the first frame of channel 1 of the point-detection demo
dataset, copied to the work folder as frame.tif (found by pattern, so no folder name is spelled
out here).
usage: python detection_c02_data.py <demo data folder> <work folder>
"""

import glob
import os
import shutil
import sys

DATA, WORK = sys.argv[1], sys.argv[2]
tifs = sorted(
    glob.glob(
        os.path.join(DATA, "PetaKit5D_point_detection_demo_dataset", "*", "*", "ch1", "*.tif")
    )
)
if not tifs:
    sys.exit("no channel-1 frame found in the point-detection demo dataset")
os.makedirs(WORK, exist_ok=True)
shutil.copyfile(tifs[0], os.path.join(WORK, "frame.tif"))
