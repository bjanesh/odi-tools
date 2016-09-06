import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pyraf import iraf


class reg_to_bpm(object):
    def __init__(self, reg_file, img_file):
        self.reg_file = reg_file
        self.img_file = img_file
        self.ota = 'OTA'+self.reg_file.strip('bpm_xy, .reg')+'.SCI'
        self.mask_name = 'bpm_ota'+self.ota+'.fits'
        try:
            hdu_list = fits.open(self.img_file)
            hdu_ota = hdu_list[self.ota]
            self.empty = np.zeros_like(hdu_ota.data)
            hdu_list.close()
        except KeyError:
            hdu_list = fits.open(self.img_file)
            hdu_ota = hdu_list['OTA33.SCI']
            self.empty = np.zeros_like(hdu_ota.data)
            hdu_list.close()


    def parse_box_reg(self):
        box_x = []
        box_y = []
        box_dx = []
        box_dy = []
        with open(self.reg_file) as reg:
            for line in reg:
                if line.startswith('box'):
                    reg_line = line[4:-1]
                    reg_info = reg_line.strip(')').split(',')
                    box_x.append(float(reg_info[0]))
                    box_y.append(float(reg_info[1]))
                    box_dx.append(float(reg_info[2]))
                    box_dy.append(float(reg_info[3]))
        self.box_x = box_x
        self.box_y = box_y
        self.box_dx = box_dx
        self.box_dy = box_dy

        return box_x,box_y,box_dx,box_dy

    def mask_reg(self,x,y,dx,dy):
        img_array = self.empty
        x1,x2 = (x-1) - 0.5*dx, x + 0.5*dx
        x1 = x1 - 1
        if x1 < 0:
            x1 = 0
        y1,y2 = (y-1) - 0.5*dy, y + 0.5*dy
        if y1 < 0:
            y1 = 0
        img_array[int(y1):int(y2),int(x1):int(x2)] = 1.0


        return img_array

    def reg_mask_ota(self):
        print self.reg_file
        box_x,box_y,box_dx,box_dy = self.parse_box_reg()
        for i,reg in enumerate(box_x):
            ota_reg_mask = self.mask_reg(self.box_x[i],
                                         self.box_y[i],
                                         self.box_dx[i],
                                         self.box_dy[i])

        hdu = fits.PrimaryHDU(ota_reg_mask.astype(float))
        mask_name = self.mask_name
        hdu.writeto(mask_name,clobber=True)
        iraf.unlearn(iraf.imutil.imcopy)
        iraf.imutil.imcopy.setParam('input',mask_name)
        iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
        iraf.imutil.imcopy.setParam('verbose','no')
        iraf.imutil.imcopy(mode='h')
        return ota_reg_mask
