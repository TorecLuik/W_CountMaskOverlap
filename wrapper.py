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
import time


def main(argv):
    with BiaflowsJob.from_cli(argv) as bj:
        
        print("Initialisation...")

        # 1. Prepare data for workflow
        # 1a. Get folders
        
        in_imgs, gt_imgs, in_path, gt_path, out_path, tmp_path = prepare_data(
            get_discipline(bj, default=CLASS_SPTCNT), bj, is_2d=True, **bj.flags)
        # 1b. MAKE SURE TMP PATH IS UNIQUE
        try:
            tmp_path += f"/{int(time.time() * 1000)}"  # timestamp in ms
            os.mkdir(tmp_path)  # setup tmp
        except FileExistsError:
            tmp_path += f"/{int(time.time() * 10000)}"  # timestamp in ms
            os.mkdir(tmp_path)  # setup tmp
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
            '''
            Finds pairs of filenames based on provided suffixes.

            Args:
                in_imgs (list): List of input filenames.
                cmsuffix (str, optional): Suffix for the first pair element. Defaults to '_C'.
                amsuffix (str, optional): Suffix for the second pair element. Defaults to '_A'.

            Returns:
                tuple: A tuple containing two elements:
                    - set: Set of base names without suffixes, for which pairs were found.
                    - list: List of pairs (filename with cmsuffix, filename with amsuffix).

            Note:
                This function assumes that pairs will have the same extension in the input filenames.

                Prints if the input filename format is invalid or if pairs are mismatched.

            Examples:
                >>> in_imgs = ['file1_C.txt', 'file1_A.txt', 'file2_A.txt']
                >>> find_pairs(in_imgs)
                Error: Mismatched or missing pair for file2. Skipping.
                ({'file1'}, [('file1_C.txt', 'file1_A.txt')])
            '''
            pairs = []
            base_names = set()
            skipped_base_names = set()

            in_imgs = [img.filename for img in in_imgs]
            for filename in in_imgs:
                try:
                    if cmsuffix in filename:
                        base_name, extension = filename.rsplit(cmsuffix, 1)
                    elif amsuffix in filename:
                        base_name, extension = filename.rsplit(amsuffix, 1)
                    else:
                        raise ValueError(f"Invalid filename format: {filename}")
                    
                    base_names.add((base_name, extension))
                except ValueError:
                    print(f"Error: Invalid filename format for {filename} for suffixes {cmsuffix}/{amsuffix}. Skipping.")
                    continue

            for base_name, extension in base_names:
                cm_filename = f'{base_name}{cmsuffix}{extension}'
                am_filename = f'{base_name}{amsuffix}{extension}'
                if cm_filename in in_imgs and am_filename in in_imgs:
                    pairs.append((cm_filename, am_filename))
                else:
                    skipped_base_names.add(base_name)
                    print(f"Error: Mismatched or missing pair for {base_name}. Skipping.")

            base_names = {item[0] for item in base_names} # keep only the actual base_name
            base_names -= skipped_base_names  # Using set difference

            return base_names, pairs
                       
        base_names, mask_pairs = find_pairs(in_imgs, cmsuffix=cmsuffix, 
                                            amsuffix=amsuffix)
        print(f"Found {len(base_names)} pairs: {base_names}")
        
        # Create an empty DataFrame to store results
        result_df = pd.DataFrame()

        # Iterate over pairs and concatenate results
        for basename, (big_mask, small_mask) in zip(base_names, mask_pairs):
            print(basename, big_mask, small_mask)
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
