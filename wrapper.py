import sys
import os
import shutil
import subprocess
import imageio
import numpy as np
import skimage
import skimage.color
from cytomine.models import Job
from biaflows import CLASS_OBJSEG, CLASS_SPTCNT, CLASS_PIXCLA, CLASS_TRETRC, CLASS_LOOTRC, CLASS_OBJDET, CLASS_PRTTRK, CLASS_OBJTRK
from biaflows.helpers import BiaflowsJob, prepare_data, upload_data, upload_metrics, get_discipline
import time
import shutil
# code for workflow:
import pandas as pd
from skimage import measure
from pyCellExpansion import CellExpansion


def main(argv):
    base_path = "{}".format(os.getenv("HOME"))  # Mandatory for Singularity
    with BiaflowsJob.from_cli(argv) as bj:
        # Change following to the actual problem class of the workflow
        problem_cls = get_discipline(bj, default=CLASS_SPTCNT)

        bj.job.update(status=Job.RUNNING, progress=0,
                      statusComment="Initialisation...")

        # 1. Prepare data for workflow
        in_imgs, gt_imgs, in_path, gt_path, out_path, tmp_path = prepare_data(
            problem_cls, bj, is_2d=True, **bj.flags)

        # MAKE SURE TMP PATH IS UNIQUE
        # tmp_path += f"/{int(time.time())}"  # timestamp
        # os.mkdir(tmp_path)  # setup tmp
        tmp_path = gt_path

        # Read parameters
        maxpixels = bj.parameters.max_pixels
        discardcellswithoutcytoplasm = bj.parameters.discard_cells_without_cytoplasm
        print(f"Max pixels: {maxpixels} |\
                Require cyto: {discardcellswithoutcytoplasm}")

        # 2. Run image analysis workflow
        bj.job.update(progress=25, statusComment="Launching workflow...")

        # Add here the code for running the analysis script
        #"--chan", "{:d}".format(nuc_channel)
        for bfimg in in_imgs:
            if "Nuclei" in bfimg.filename:
                print(f"CellExpand: {bfimg.__dict__}")
                # Read Nuclei labels
                fn = os.path.join(in_path, bfimg.filename)
                imCellsNucleiLabels = imageio.imread(fn)
                # Expand Cells
                (imCellsNucleiLabels,
                 imCellsCellLabels,
                 imCellsCytoplasmLabels) = CellExpansion(
                    imCellsNucleiLabels=imCellsNucleiLabels,
                    discardcellswithoutcytoplasm=discardcellswithoutcytoplasm,
                    maxpixels=maxpixels)
                # Write intermediate results
                imageio.imwrite(
                    os.path.join(tmp_path,
                                 bfimg.filename.replace("Nuclei", "Cells")),
                    imCellsCellLabels)
                imageio.imwrite(
                    os.path.join(tmp_path,
                                 bfimg.filename.replace("Nuclei", "Cytoplasm")),
                    imCellsCytoplasmLabels)
                print(f"Wrote 2 masks to {tmp_path}")

                # Read Granules
                imCellsGranulesLabels = imageio.imread(os.path.join(in_path, 
                    bfimg.filename.replace("Nuclei", "Granules")))
                # Calc overlap
                numCells = np.max(imCellsCellLabels)
                CellNumGranules = np.zeros([numCells, 1], dtype=np.int16)
                granulesStats = pd.DataFrame(
                    measure.regionprops_table(
                        imCellsGranulesLabels, properties=('centroid', )))
                granulesStatsnp = np.ndarray.astype(
                    np.round(granulesStats.to_numpy()), dtype=np.uint16)
                granulesStatsInCellLabel = imCellsCellLabels[
                    granulesStatsnp[:, 0], granulesStatsnp[:, 1]]
                for i in range(1, numCells+1):
                    CellNumGranules[i-1, 0] = np.count_nonzero(
                        granulesStatsInCellLabel == i)
                print(f"Counts: {CellNumGranules}")
                pd.DataFrame(CellNumGranules, columns=["EstimatedGranules"]).to_csv(
                    os.path.join(tmp_path,
                                 bfimg.filename_no_extension+"_granulecounts.csv")
                )
                print(f"Wrote CSV to {tmp_path}")

        # Copy to out folder when we're done
        for bimg in in_imgs:
            if "Nuclei" in bimg.filename:
                shutil.copy(
                    os.path.join(tmp_path, bimg.filename.replace("Nuclei",
                                                                 "Cells")),
                    out_path)
                shutil.copy(
                    os.path.join(tmp_path, bimg.filename.replace("Nuclei",
                                                                 "Cytoplasm")),
                    out_path)

                shutil.copy(
                    os.path.join(
                        tmp_path, bimg.filename_no_extension+"_granulecounts.csv"),
                    out_path)
                print(f"Copied tmp files to {out_path}")

        # 3. Upload data to BIAFLOWS
        upload_data(problem_cls, bj, in_imgs, out_path, **bj.flags, monitor_params={
            "start": 60, "end": 90, "period": 0.1,
            "prefix": "Extracting and uploading polygons from masks"})

        # 4. Compute and upload metrics
        bj.job.update(
            progress=90, statusComment="Computing and uploading metrics...")
        upload_metrics(problem_cls, bj, in_imgs, gt_path,
                       out_path, tmp_path, **bj.flags)

        # 5. Pipeline finished
        shutil.rmtree(tmp_path)  # cleanup tmp
        bj.job.update(progress=100, status=Job.TERMINATED,
                      status_comment="Finished.")


if __name__ == "__main__":
    main(sys.argv[1:])
