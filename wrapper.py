import random
import sys
import os
import shutil
import imageio
import numpy as np
import skimage
from biaflows import CLASS_SPTCNT
from biaflows.helpers import BiaflowsJob, prepare_data, get_discipline
# code for workflow:
from pyCountOverlap import count_overlap
import pandas as pd


def main(argv):
    with BiaflowsJob.from_cli(argv) as bj:
        
        print("Initialisation...")

        # 1. Prepare data for workflow
        # 1a. Get folders
        
        in_imgs, gt_imgs, in_path, gt_path, out_path, tmp_path = prepare_data(
            get_discipline(bj, default=CLASS_SPTCNT), bj, is_2d=True, **bj.flags)
        print(in_path, in_imgs)
        print(os.listdir(in_path), os.listdir('/data/in'))
        in_imgs = os.listdir('/data/in')
        # 1b. MAKE SURE TMP PATH IS UNIQUE
        tmp_path = gt_path
        # 1c. Read parameters from commandline
        
        cmsuffix = bj.parameters.cell_mask_suffix
        amsuffix = bj.parameters.aggregate_mask_suffix
        c_counts = bj.parameters.column_name_counts
        c_cells = bj.parameters.column_name_cells
        
        print(f"Parameters: Cell suffix: {cmsuffix} |\
                Aggregate suffix: {amsuffix}")

        # 2. Run image analysis workflow
        print("Launching workflow...")

        # 2a. Add here the code for running the analysis script     
        def find_pairs(in_imgs, cmsuffix='_C', amsuffix='_A'):
            pairs = []
            base_names = set()
            skipped_base_names = set()
            
            for filename in in_imgs:
                try:
                    base_name, extension = filename.rsplit('.', 1)
                    print(base_name, extension)
                    if base_name.endswith(cmsuffix):
                        base_names.add(base_name[:-len(cmsuffix)])
                    print(base_names)
                except ValueError:
                    print(f"Error: Invalid filename format for {filename}. Skipping.")
                    continue
                
            for base_name in base_names:
                cm_filename = f'{base_name}{cmsuffix}.{extension}'
                am_filename = f'{base_name}{amsuffix}.{extension}'
                
                if cm_filename in in_imgs and am_filename in in_imgs:
                    pairs.append((cm_filename, am_filename))
                else:
                    skipped_base_names.add(base_name)
                    print(f"Error: Mismatched or missing pair for {base_name}. Skipping.")
            
            base_names -= skipped_base_names  # Using set difference
            
            return base_names, pairs
        
        print(f"IN: {in_imgs}. Suffixes: {cmsuffix}, {amsuffix}")
        base_names, mask_pairs = find_pairs(in_imgs, cmsuffix=cmsuffix, 
                                            amsuffix=amsuffix)
        print(f"Found {len(base_names)} pairs: {base_names}")
        
        # Create an empty DataFrame to store results
        result_df = pd.DataFrame()

        # Iterate over pairs and concatenate results
        for basename, (big_mask, small_mask) in zip(base_names, mask_pairs):
            big_mask = os.path.join(in_path, big_mask)
            small_mask = os.path.join(in_path, small_mask)
            df = count_overlap(big_mask, small_mask, columnName=c_counts)
            df['Basefilename'] = basename
            df = df.reset_index()
            df = df.rename(columns={"index": c_cells})
            result_df = pd.concat([result_df, df.set_index('Basefilename')])
            print(f"Counted {basename}: {result_df}")
            
        # Write the result_df to a CSV file
        result_df.to_csv(os.path.join(out_path, 'counts.csv'), index=True)

        print("CSV file written successfully.")

        # 3. Pipeline finished
        # 3a. cleanup tmp
        shutil.rmtree(tmp_path)  
        print("Finished.")


if __name__ == "__main__":
    main(sys.argv[1:])
