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


def main(argv):
    base_path = "{}".format(os.getenv("HOME")) # Mandatory for Singularity
    with BiaflowsJob.from_cli(argv) as bj:
        # Change following to the actual problem class of the workflow
        problem_cls = get_discipline(bj, default=CLASS_OBJSEG)
        
        bj.job.update(status=Job.RUNNING, progress=0, statusComment="Initialisation...")
        
        # 1. Prepare data for workflow
        in_imgs, gt_imgs, in_path, gt_path, out_path, tmp_path = prepare_data(problem_cls, bj, is_2d=True, **bj.flags)
        
        # MAKE SURE TMP PATH IS UNIQUE
        tmp_path += f"/{int(time.time())}"  # timestamp
        os.mkdir(tmp_path)  # setup tmp

        # Make sure all images have at least 224x224 dimensions
        # and that minshape / maxshape * minshape >= 224
        # 0 = Grayscale (if input RGB, convert to grayscale)
        # 1,2,3 = rgb channel
        nuc_channel = bj.parameters.nuc_channel
        resized = {}
        for bfimg in in_imgs:
            fn = os.path.join(in_path, bfimg.filename)
            img = imageio.imread(fn)
            if len(img.shape) > 2 and nuc_channel == 0:
                gray_rgb = False
                if np.array_equal(img[:,:,0],img[:,:,1]) and np.array_equal(img[:,:,0],img[:,:,2]):
                    gray_rgb = True
                img = skimage.color.rgb2gray(img) * 255
                img = img.astype(np.uint8)
                # Invert intensity if not grayscale img ie. expect the image
                # to be H&E stained image with dark nuclei
                if not gray_rgb:
                    img = np.invert(img)
            minshape = min(img.shape[:2])
            maxshape = max(img.shape[:2])
            if minshape != maxshape or minshape < 224:
                resized[bfimg.filename] = img.shape
                padshape = []
                for i in range(2):
                    if img.shape[i] < max(224,maxshape):
                        padshape.append((0,max(224,maxshape)-img.shape[i]))
                    else:
                        padshape.append((0,0))
                if len(img.shape) == 3:
                    padshape.append((0,0))
                img = np.pad(img, padshape, 'constant', constant_values=0)
            imageio.imwrite(os.path.join(tmp_path, bfimg.filename), img)

        # 2. Run image analysis workflow
        bj.job.update(progress=25, statusComment="Launching workflow...")

        # Add here the code for running the analysis script
        #"--chan", "{:d}".format(nuc_channel)
        prob_thresh = bj.parameters.prob_threshold
        diameter = bj.parameters.diameter
        cp_model = bj.parameters.cp_model
        use_gpu = bj.parameters.use_gpu
        print(f"Chosen model: {cp_model} | Channel {nuc_channel} | Diameter {diameter} | Cell prob threshold {prob_thresh} | GPU {use_gpu}")
        cmd = ["python", "-m", "cellpose", "--dir", tmp_path, "--pretrained_model", f"{cp_model}", "--save_tif", "--no_npy", "--chan", "{:d}".format(nuc_channel), "--diameter", "{:f}".format(diameter), "--cellprob_threshold", "{:f}".format(prob_thresh)]
        if use_gpu:
            print("Using GPU!")
            cmd.append("--use_gpu")
        status = subprocess.run(cmd)

        if status.returncode != 0:
            print("Running Cellpose failed, terminate")
            sys.exit(1)

        # Crop to original shape
        for bimg in in_imgs:
            shape = resized.get(bimg.filename, None)
            if shape:
                img = imageio.imread(os.path.join(tmp_path,bimg.filename_no_extension+"_cp_masks.tif"))
                img = img[0:shape[0], 0:shape[1]]
                imageio.imwrite(os.path.join(out_path,bimg.filename), img)
            else:
                shutil.copy(os.path.join(tmp_path,bimg.filename_no_extension+"_cp_masks.tif"), os.path.join(out_path,bimg.filename))
        
        # 3. Upload data to BIAFLOWS
        upload_data(problem_cls, bj, in_imgs, out_path, **bj.flags, monitor_params={
            "start": 60, "end": 90, "period": 0.1,
            "prefix": "Extracting and uploading polygons from masks"})
        
        # 4. Compute and upload metrics
        bj.job.update(progress=90, statusComment="Computing and uploading metrics...")
        upload_metrics(problem_cls, bj, in_imgs, gt_path, out_path, tmp_path, **bj.flags)

        # 5. Pipeline finished
        shutil.rmtree(tmp_path)  # cleanup tmp
        bj.job.update(progress=100, status=Job.TERMINATED, status_comment="Finished.")


if __name__ == "__main__":
    main(sys.argv[1:])
