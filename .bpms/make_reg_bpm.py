from bpm_reg_class import *
import glob

regfiles = glob.glob('*.reg')

for reg_file in regfiles:
    region = reg_to_bpm(reg_file,'20140406T214040.1_GCPair-F1_odi_g.6183.fits')
    region.reg_mask_ota()
