import argparse
def runinbatch(foldd,pid):
    import os
    import pathlib
    sss1  = '#SBATCH -p short\n#SBATCH --mem=30G\n#SBATCH --time 00:30:00\n'
    sss2  = '#SBATCH -p long\n#SBATCH --mem=50G\n#SBATCH --time 00:30:00\n'
    divice ='cpu'
    patha = pathlib.Path(foldd)
    flist  = patha.glob('*ch_0.tiff')
    idis = []
    for file in flist:
        if os.path.isfile(f'{str(file).split(".")[0]}.h5') and pathlib.Path(f'{str(file).split(".")[0]}.h5').stat().st_size > 10000:
            ned = os.popen(
                f"sbatch << EOF \n#!/bin/bash\n#SBATCH -p short\n#SBATCH --mem=5G\n#SBATCH --time 01:00:00\n#SBATCH --job-name={str(file).split('.')[0].split(os.sep)[-1]}_s2\n#SBATCH -o /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/surmout/{str(file).split('.')[0].split(os.sep)[-1]}_s2.out\npython /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/hpc_iman.py -i {file}\nEOF").read()
        else:
            id = os.popen(f"sbatch << EOF \n#!/bin/bash\n{sss1}#SBATCH --job-name={str(file).split('.')[0].split(os.sep)[-1]}_s1\n#SBATCH -o /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/surmout/{str(file).split('.')[0].split(os.sep)[-1]}_s1.out\npython /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/YeaZ-GUI-master/Launch_NN_command_line.py -i {str(file)} -m {str(file).split('.')[0]}.h5 --device {divice} --path_to_weights /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/YeaZ-GUI-master/unet/weights_budding_BF_multilab_0_1\nEOF").read()
            id = id.split(' ')[3]
            id2 = os.popen(f"sbatch << EOF \n#!/bin/bash\n{sss2}#SBATCH --dependency=afternotok:{id}\n#SBATCH --job-name={str(file).split('.')[0].split(os.sep)[-1]}_s1\n#SBATCH -o /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/surmout/{str(file).split('.')[0].split(os.sep)[-1]}_s1.out\npython /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/YeaZ-GUI-master/Launch_NN_command_line.py -i {str(file)} -m {str(file).split('.')[0]}.h5 --device {divice} --path_to_weights /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/YeaZ-GUI-master/unet/weights_budding_BF_multilab_0_1\nEOF").read()
            id2 = id2.split(' ')[3]
            ned = os.popen(
                    f"sbatch << EOF \n#!/bin/bash\n#SBATCH -p short\n#SBATCH --mem=5G\n#SBATCH --time 01:00:00\n#SBATCH  --dependency=afterok:{id}?{id2}\n#SBATCH --job-name={str(file).split('.')[0].split(os.sep)[-1]}_s3\n#SBATCH -o /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/surmout/{str(file).split('.')[0].split(os.sep)[-1]}_s3.out\npython /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/hpc_iman.py -i {file}\nEOF").read()
        idis.append(ned)
    listofiseis = ','.join(*idis)
    id3 = os.popen(f"sbatch << EOF \n#!/bin/bash\n#SBATCH -p short\n#SBATCH --mem=5G\n#SBATCH --time 01:00:00\n#SBATCH  --dependency=afterok:{listofiseis}\n#SBATCH --job-name=merge\n#SBATCH -o /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/surmout/merge{pid}.out\npython /data/gpfs-1/users/cohenn_c/work/microscopy/microscopy/comine_simple.py -i {foldd} -p {pid}\nEOF").read()
    return None

def nd2tiff_s1(filelocation,pid):
    import cv2 as cv
    import pandas as pd
    import pathlib
    import nd2, os
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    path = pathlib.Path(filelocation)
    basefolder = path.parent
    pathlib.Path(basefolder.joinpath(f'{pid}_tiff')).mkdir(parents=True, exist_ok=True)
    with nd2.ND2File(filelocation) as ndfile:
        microscopesetup = pd.DataFrame(ndfile.text_info['description'].split('\r\n'),columns=['setup'])
        microscopesetup.to_csv(basefolder.joinpath(pid+"_microscopesetup.txt"),sep='\t')
    image_array =  nd2.imread(filelocation, xarray=True, dask=True)
    number_of_channels = image_array.shape[1]
    number_of_pos = image_array.shape[0]
    logger.debug(f'Reading the ND file into an xarray with shape: {image_array.shape}')
    for fov in range(0,number_of_pos):
        for channelg in range(0,number_of_channels):
            logger.debug(f"processing {pid}_p{fov}_ch_{str(channelg)}.tiff")
            imagee  = image_array[fov,channelg,:,:].compute()
            logger.debug(f"Gathering metadata for {pid}_p{fov}_ch_{str(channelg)}.tiff")
            name = imagee.coords['P'].values.tolist()
            well = name.split('_')[0]
            pos = int(name.split('_')[1])+1
            if os.path.isfile(str(basefolder.joinpath(f'{pid}_tiff',f'{pid}_{well}_p{pos}_ch_{str(channelg)}.tiff'))):
              pass
            else:
              logger.debug(f"Saving {pid}_{well}_p{pos}_ch_{str(channelg)}.tiff")
              cv.imwrite(str(basefolder.joinpath(f'{pid}_tiff',f'{pid}_{well}_p{pos}_ch_{str(channelg)}.tiff')), imagee.values)
              logger.debug(f"Saved {pid}_{well}_p{pos}_ch_{str(channelg)}.tiff")
    
    runinbatch(str(basefolder.joinpath(f'{pid}_tiff')),pid)
    return None
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process nd2 file using candida pipeline')
    parser.add_argument('-i', '--input', type=str, help='input ND2 file')
    parser.add_argument('-p', '--pid', type=str, help='plate id')
    args = parser.parse_args()
    nd2tiff_s1(args.input, args.pid)