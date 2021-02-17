"""shared io related variables"""

from pixell import enmap
import os.path as op
import glob

# Input files
map_dir = "/projects/ACT/yilung/work/galactic_center"
filedb = {
    'f090': {
        'coadd':  op.join(map_dir, 'coadd', 'gc_f090_map.fits'),
        'planck': op.join(map_dir, 'planck', 'planck_npipe_100_split*_map.fits')
    },
    'f150': {
        'coadd':  op.join(map_dir, 'coadd', 'gc_f150_map.fits'),
        'planck': op.join(map_dir, 'planck', 'planck_npipe_143_split*_map.fits')        
    },
    'f220': {
        'coadd':  op.join(map_dir, 'coadd', 'gc_f220_map.fits'),
        'planck': op.join(map_dir, 'planck', 'planck_npipe_217_split*_map.fits')        
    },
}

def load_map(path, box=None):
    files = glob.glob(path)
    if len(files) == 1:
        imap = enmap.read_map(files[0])
    elif len(files) == 2:
        imap1 = enmap.read_map(files[0])
        imap2 = enmap.read_map(files[1])
        imap  = (imap1 + imap2) / 2
    else: raise ValueError("Unknown format")
    if box is not None: return imap.submap(box)
    else: return imap
